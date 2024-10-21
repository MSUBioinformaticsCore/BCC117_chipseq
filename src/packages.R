if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

bioc = c("org.At.tair.db", "topGO", "AnnotationHub", 
         "GenomicRanges", "plyranges", "profileplyr",
         "DiffBind")

install.packages("tidyverse")