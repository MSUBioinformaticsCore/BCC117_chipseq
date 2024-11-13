#`getExclusionList`
#https://dozmorovlab.github.io/excluderanges/#cutrun-excludable-sets

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("AnnotationHub", update = FALSE) 
# BiocManager::install("GenomicRanges", update = FALSE)
# BiocManager::install("plyranges", update = FALSE)

#' @param args[1] outdir

args <- commandArgs(TRUE)

suppressMessages(library(GenomicRanges))
suppressMessages(library(AnnotationHub))
suppressMessages(library(tidyverse))

ah <- AnnotationHub()
query_data <- subset(ah, preparerclass == "excluderanges")

exclude <- query_data[["AH107341"]]

exclude <- 
  exclude %>% 
  sort() %>% 
  keepStandardChromosomes(pruning.mode = "tidy")

exclude = as.data.frame(exclude)
exclude$seqnames = gsub("chr", "", exclude$seqnames)

mt_pt = data.frame(seqnames = c("Mt", "Pt"),
                   start = c(1,1),
                   end = c(366924, 154478),
                   width = c(366924, 154478),
                   strand = c("*", "*"),
                   name = c("Mitochondria", "Chloroplast"))

exclude = rbind(exclude, mt_pt)

write.table(as.data.frame(exclude), 
            file = paste0(args[1], "/TAIR10.Klasfeld.arabidopsis_greenscreen_20inputs_MtPt.bed"),
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)