#`getExclusionList`
#https://dozmorovlab.github.io/excluderanges/#cutrun-excludable-sets

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("AnnotationHub", update = FALSE) 
# BiocManager::install("GenomicRanges", update = FALSE)
# BiocManager::install("plyranges", update = FALSE)

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

write.table(as.data.frame(exclude), 
            file = "/mnt/research/bioinformaticsCore/shared/Genomes/arabidopsis_thaliana/TAIR10.Klasfeld.arabidopsis_greenscreen_20inputs.bed",
            sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)