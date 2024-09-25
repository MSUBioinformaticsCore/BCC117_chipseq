#' @param args[1] path to sample sheet
#' @param args[2] path to results dir
#' @param args[3] use grey list T/F
#' https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf

args <- commandArgs(TRUE)

# install -----------------------------------------------------------------

installed = installed.packages()[,"Package"]

# DiffBind

if(!"DiffBind" %in% installed) {
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("DiffBind")
  
}

# profileplyr

if(!"profileplyr" %in% installed) {
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("profileplyr")
  
}

# tidyverse

if(!"tidyverse" %in% installed) {
  
  install.packages("tidyverse")
  
}

# parallel

if(!"parallel" %in% installed) {
  
  install.packages("parallel")
  
}

library(DiffBind)
library(tidyverse)
library(parallel)

# set up ------------------------------------------------------------------

results.dir = args[2]
# results.dir = "/mnt/research/bioinformaticsCore/projects/grossmanl/20221220_ChIPSeq/results"
diffbind.dir = paste0(results.dir, "/diffbind")
if(!dir.exists(diffbind.dir)){dir.create(diffbind.dir)}

# Load sample sheet  -------------------------------------------------------

# samples = read.csv("/mnt/research/bioinformaticsCore/projects/grossmanl/20221220_ChIPSeq/data/DiffBind_samplesheet_consensus.csv")
samples = read.csv(args[1])

# DBA for each comparison -------------------------------------------------

sample_list = list()
dbaList = list()

conditions = unique(samples$Condition)
comparisons = combn(conditions, 2, simplify = FALSE)
comparisons = lapply(comparisons, function(x) {x[order(x)]})

for (i in 1:length(comparisons)){
  
  sample_list[[i]] = 
    samples %>%
    filter(Condition %in% comparisons[[i]])
  
  dbaList[[i]] = dba(sampleSheet=sample_list[[i]], bRemoveM=FALSE)

}

# count and normalize ----------------------------------------

# dbaOb <- dba.count(dbaOb)

dbaList = mclapply(dbaList, 
                   dba.count, 
                   mc.cores = detectCores()-1)

# dbaOb <- dba.normalize(dbaOb)

dbaList = mclapply(dbaList, 
                   dba.normalize, 
                   mc.cores = detectCores()-1)

# model design ------------------------------------------------------------

dbaList = mclapply(dbaList, 
                   dba.contrast, 
                   mc.cores = detectCores()-1)

# Differential binding ----------------------------------------------------

if(args[3]){
  
  print(paste0("use greylist =", args[3]))
  
  #dbaOb <- dba.analyze(dbaOb, bGreylist=TRUE, bBlacklist=FALSE)
  
  dbaList = mclapply(dbaList, 
                     dba.analyze, 
                     bGreylist=TRUE,
                     bBlacklist=FALSE,
                     mc.cores = detectCores()-1)
  
  file_tag="_greylisted"
  
} else {
  
  print(paste0("use greylist =", args[3]))
  
  dbaList = mclapply(dbaList, 
                     dba.analyze, 
                     bGreylist=FALSE,
                     bBlacklist=FALSE,
                     mc.cores = detectCores()-1)
  file_tag=""
  
}

# save DBA ----------------------------------------------------------------

save(dbaList, file = paste0(diffbind.dir, "/DBAlist", file_tag, ".Rdata"))

# output results ----------------------------------------------------------

res_list = list()

for (i in 1:length(dbaList)) {
  
  group1 = comparisons[[i]][1]
  group2 = comparisons[[i]][2]

  res_list[[i]] = dba.report(dbaList[[i]], 
                            contrast=1, 
                            th = .05, 
                            fold = 0, 
                            DataType = "DBA_DATA_FRAME")
  
  names(res_list[[i]])[names(res_list[[i]]) == 'p-value'] <- 'pvalue'
  
  res_list[[i]] =
    res_list[[i]] %>%
    mutate(Contrast = paste0(group1, "_vs_", group2)) %>%
    mutate(peakID=row_number()) %>%
    mutate(peakID=paste0(group1, "_vs_", group2, "_peak", peakID))
   
  write.csv(res_list[[i]], 
              file = paste0(diffbind.dir, 
                            "/", group1, "_vs_", group2,
                            file_tag, ".csv"))
  
  print(paste0(diffbind.dir, 
               "/", group1, "_vs_", group2,
               file_tag, ".csv saved"))
  
  bed = 
    res_list[[i]] %>%
    mutate(Strand = ".") %>%
    mutate(score = -10*log10(FDR)) %>%
    select(Chr,
           Start,
           End,
           peakID,
           score,
           Strand,
           Fold,
           pvalue,
           FDR)
  
  write.table(bed, 
            file = paste0(diffbind.dir, 
                          "/", group1, "_vs_", group2,
                          file_tag, ".bed"),
            sep="\t",
            quote=F,
            row.names=F,
            col.names=F)
  
  print(paste0(diffbind.dir, 
               "/", group1, "_vs_", group2,
               file_tag, ".bed saved"))
  }

# volcano plots -----------------------------------------------------------

for (i in 1:length(dbaList)){
  
  group1 = comparisons[[i]][1]
  group2 = comparisons[[i]][2]
  
  pdf(file = paste0(diffbind.dir, "/", group1, "_vs_", group2, file_tag, "_volcano_plot.pdf"))
  dba.plotVolcano(dbaList[[i]],
                  contrast=1, 
                  th = .05)
  dev.off()
  print(paste0(diffbind.dir, "/", group1, "_vs_", group2, file_tag, "_volcano_plot.pdf saved"))
}

# Profile plots -----------------------------------------------------------

library(profileplyr)

profiles = mclapply(dbaList, dba.plotProfile, sites=1, merge=FALSE)

for (i in 1:length(profiles)){
  
  group1 = comparisons[[i]][1]
  group2 = comparisons[[i]][2]
  
  labs = c(paste0(group1, "_rep1"),
           paste0(group1, "_rep2"),
           paste0(group1, "_rep3"),
           paste0(group2, "_rep1"),
           paste0(group2, "_rep2"),
           paste0(group2, "_rep3"))
  
  rownames(sampleData(profiles[[i]])) = labs
             
  pdf(file = paste0(diffbind.dir, 
                    "/",
                    group1, 
                    "_vs_", 
                    group2, 
                    file_tag, "_profile_plot.pdf"),
      width = 8)
  dba.plotProfile(profiles[[i]])
  dev.off()
  print(paste0(diffbind.dir, "/", group1, "_vs_", group2, file_tag, "_profile_plot.pdf saved"))
  
}

# session info ------------------------------------------------------------

sessionInfo()
