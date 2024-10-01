# split diffbind anno files into gain/loss files
#' @param args[1] path to diffbind files
#' @param args[2] path to anno files
args <- commandArgs(TRUE)

library(tidyverse)

# set up ------------------------------------------------------------------

dbPath = args[1]
annoPath = args[2]

#dbPath = "/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/mergedLibrary/diffbind"
#annoPath = "/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/mergedLibrary/annotatePeaks"

dbFiles = list.files(dbPath, pattern = ".csv", full.names=T)
dbList = lapply(dbFiles, read.csv)

dbBedFiles = list.files(dbPath, pattern = ".bed", full.names=T)
dbBedList = lapply(dbBedFiles, read.delim, header = F)

annoFiles = list.files(annoPath, pattern= "_vs_", full.names=T)
annoFiles = annoFiles[grep("anno", basename(annoFiles))]
annoList = lapply(annoFiles, read.delim)

# rename PeakID cols ---------------------------------------------------

dbList = lapply(dbList, function(df){
  df = 
    df %>%
    dplyr::rename(PeakID = "peakID");
  df
})

annoList = lapply(annoList, function(df){
  df = 
    df %>%
    dplyr::rename(PeakID = colnames(df)[1]) 
})

# Add gene names to DB files and split by fc sign ----------------------

gene2peak = lapply(annoList, function(df){
  df = 
    df %>%
    select(PeakID, Entrez.ID, Gene.Name)
})

gainFiles = gsub(".bed.anno", "_gain.anno", annoFiles)
lossFiles = gsub(".bed.anno", "_loss.anno", annoFiles)

gainBedFiles = gsub(".bed", "_gain.bed", dbBedFiles)
lossBedFiles = gsub(".bed", "_loss.bed", dbBedFiles)

## change the peak names
for (i in 1:length(dbList)){
  
  dbList[[i]] = left_join(dbList[[i]], gene2peak[[i]])
  write.csv(dbList[[i]], row.names = FALSE, file = dbFiles[i])
  
  gainPeaks = 
    dbList[[i]] %>%
    filter(Fold > 0) %>%
    pull(PeakID)
  
  lossPeaks = 
    dbList[[i]] %>%
    filter(Fold < 0) %>%
    pull(PeakID)
  
  gainBed = 
    dbBedList[[i]] %>%
    filter(V4 %in% gainPeaks)
  
  write.table(gainBed, 
              row.names = FALSE,
              col.name = FALSE,
              sep = "\t",
              quote=F,
              file = gainBedFiles[i])
  
  lossBed = 
    dbBedList[[i]] %>%
    filter(V4 %in% lossPeaks)
  
  write.table(lossBed, 
              row.names = FALSE,
              col.name = FALSE,
              sep = "\t",
              quote=F,
              file = lossBedFiles[i])
  
  gainAnno = 
    annoList[[i]] %>%
    filter(PeakID %in% gainPeaks)
  
  write.table(gainAnno, 
              row.names = FALSE,
              sep = "\t",
              quote=F,
              file = gainFiles[i])
  
  lossAnno = 
    annoList[[i]] %>%
    filter(PeakID %in% lossPeaks)
  
  write.table(lossAnno, 
              row.names = FALSE,
              sep = "\t",
              quote=F,
              file = lossFiles[i])

}


