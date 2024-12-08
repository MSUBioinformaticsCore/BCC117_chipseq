#' @param args[1] path to sample sheet
#' @param args[2] path to results dir
#' @param args[3] baseline condition
#' @param args[4] greylist 
#' https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf

args <- commandArgs(TRUE)

# install ----------------------------------------------------------------------

installed = installed.packages()[,"Package"]

# DiffBind

if(!"DiffBind" %in% installed) {
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("DiffBind")
  
}

# profileplyr

if(!"profileplyr" %in% installed) {
  
  install_version("ggnewscale", version = "0.4.9")
  
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  
  BiocManager::install("profileplyr")
  
}

# tidyverse

if(!"tidyverse" %in% installed) {
  
  install.packages("tidyverse")
  
}


library(DiffBind)
library(tidyverse)

# set up -----------------------------------------------------------------------

results.dir = args[2]
# results.dir = "/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/merged_library"

diffbind.dir = paste0(results.dir, "/diffbind")
if(!dir.exists(diffbind.dir)){dir.create(diffbind.dir)}

baseline_condition = args[3]
# baseline_condition = "warm"

# Load sample sheet and greylist -----------------------------------------------

# sample_sheet = read.csv("/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/data/sample_sheets/DiffBind_sample_sheet_narrowPeak.csv")
sample_sheet = read.csv(args[1])

# custom_greylist = readRDS("/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/merged_library/greylist.Rds")
custom_greylist = readRDS(args[4])

# analyse ----------------------------------------------------------------------

dbaOb = dba(sampleSheet= sample_sheet, bRemoveM=FALSE) %>%
  dba.blacklist(blacklist = FALSE, greylist=custom_greylist$master) %>%
  dba.count() %>%
  dba.normalize() %>%
  dba.contrast(reorderMeta=list(Condition=baseline_condition)) %>%
  dba.analyze(bBlacklist = FALSE, bGreylist = FALSE)

# save DBA ---------------------------------------------------------------------

save(dbaOb, file = paste0(diffbind.dir, "/DBA.Rdata"))

# output results ----------------------------------------------------------

 for (i in 1:length(dbaOb$contrasts)) {
  
  group1 = dbaOb$contrasts[[i]]$name1
  group2 = dbaOb$contrasts[[i]]$name2
  
  res = dba.report(dbaOb, 
                   contrast=i, 
                   th = .05, 
                   fold = 0, 
                   DataType = "DBA_DATA_FRAME")
  
  names(res)[names(res) == 'p-value'] <- 'pvalue'
  
  res =
    res %>%
    mutate(Contrast = paste0(group1, "_vs_", group2)) %>%
    mutate(peakID=row_number()) %>%
    mutate(peakID=paste0(group1, "_vs_", group2, "_peak", peakID))
  
  write.csv(res, 
            file = paste0(diffbind.dir, 
                          "/", group1, "_vs_", group2, ".csv"))
  
  print(paste0(diffbind.dir, 
               "/", group1, "_vs_", group2, ".csv saved"))
  
  bed = 
    res %>%
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
                            "/", group1, "_vs_", group2, ".bed"),
              sep="\t",
              quote=F,
              row.names=F,
              col.names=F)
  
  print(paste0(diffbind.dir, 
               "/", group1, "_vs_", group2, ".bed saved"))
 }

# volcano plots -----------------------------------------------------------

for (i in 1:length(dbaOb$contrasts)){
  
  group1 = dbaOb$contrasts[[i]]$name1
  group2 = dbaOb$contrasts[[i]]$name2
  
  pdf(file = paste0(diffbind.dir, "/", group1, "_vs_", group2, "_volcano_plot.pdf"))
  dba.plotVolcano(dbaOb,
                  contrast=i, 
                  th = .05)
  dev.off()
  print(paste0(diffbind.dir, "/", group1, "_vs_", group2, "_volcano_plot.pdf saved"))
}

# Profile plots -----------------------------------------------------------
# https://content.cruk.cam.ac.uk/bioinformatics/software/DiffBind/plotProfileDemo.html

library(profileplyr)

profiles_list = list() 

for (i in 1:length(dbaOb$contrasts)){
  
  rep = dba.report(dbaOb, 
                   contrast=i, 
                   th = .05) 
  
  repList <- GRangesList(Gain=rep[rep$Fold>0,],Loss=rep[rep$Fold<0,])
  
  sample_names = names(which(dbaOb$contrasts[[i]]$group1 | dbaOb$contrasts[[i]]$group2 == TRUE))
  
  profiles_list[[i]] = dba.plotProfile(dbaOb, 
                                       samples = dbaOb$contrasts[[i]]$group1 | dbaOb$contrasts[[i]]$group2,
                                       sites=repList, 
                                       scores="Fold",
                                       merge=FALSE)
  
  rownames(sampleData(profiles_list[[i]])) = sample_names
  
  group1 = dbaOb$contrasts[[i]]$name1
  group2 = dbaOb$contrasts[[i]]$name2
  
  pdf(file = paste0(diffbind.dir, 
                    "/",
                    group1, 
                    "_vs_", 
                    group2, "_profile_plot.pdf"),
      width = 8)
  dba.plotProfile(profiles_list[[i]])
  dev.off()
  
  print(paste0(diffbind.dir, "/", group1, "_vs_", group2, "_profile_plot.pdf saved"))
  
}

save(profiles_list, file = paste0(diffbind.dir,"profiles.Rdata"))
     
# session info ------------------------------------------------------------
sessionInfo()
     