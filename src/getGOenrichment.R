#' @param args[1] path to ChIPSeqFUnctions.R
#' @param args[2] path to input dir
#' @param args[3] path to output dir
#' @param args[4] file extension
#' @param args[5] orgdb from bioconductor
#' @param args[6] path to Homer gene file

args <- commandArgs(TRUE)

# functions ---------------------------------------------------------------

processAnno = function(annoDF){
  annoDF = annoDF %>% dplyr::rename(PeakID = colnames(annoDF)[1])
  #annoDF = left_join(annoDF,gene_map)
  return(annoDF)
}

addOntologyCol = function(go, ont){
  go$Ontology = ont
  return(go)
}

# Set up ------------------------------------------------------------------

set.seed(1)
source(args[1])
# source("/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/src/ChIPSeqFunctions.R")
library(tidyverse)
library(topGO)
library(parallel)

orgdb = args[5]
# orgdb = "org.At.tair.db"
library(orgdb, character.only=TRUE)

indir = args[2]
# indir = "/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/mergedLibrary/annotatePeaks"
outdir = args[3]
# outdir = "/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/mergedLibrary"
enrichdir = paste0(outdir, "/geneset_enrichment")
if(!dir.exists(enrichdir)){dir.create(enrichdir)}
# genes = read.delim("/mnt/research/bioinformaticsCore/software/HOMER/data/accession/arabidopsis.description")
genes = read.delim(args[6], header=T)
# genes = genes %>% dplyr::rename(Entrez.ID = "GeneID")
# gene_map = genes %>% dplyr::select(Entrez.ID, homer_symbols, current_symbols)

# read in annotation files ------------------------------------------------

ext = args[4]
# ext = ".anno"
files = list.files(indir, pattern = ext, full.names = T)
annoList = mclapply(files, read.delim, mc.cores = detectCores()-1)
annoList = mclapply(annoList, processAnno, mc.cores = detectCores()-1)

# get enriched go ---------------------------------------------------------

goiList = mclapply(annoList, 
                   function(annoDF){goi = annoDF %>% pull(Gene.Name); 
                                    unique(goi)},
                   mc.cores = detectCores()-1)

# these are throwing errors with mclapply for some reason :-/
gobpList = lapply(goiList, 
                    findEnrichedGOBP,
                    background = genes$name,
                    org.db = orgdb,
                    id_type = "symbol",
                    min_size = 5,
                    max_size = 200)

print("findEnrichedGOBP finished")

gomfList = lapply(goiList, 
                  findEnrichedGOMF,
                  background = genes$name,
                  org.db = orgdb,
                  id_type = "symbol",
                  min_size = 5,
                  max_size = 200)

print("findEnrichedGOMF finished")

goccList = lapply(goiList, 
                  findEnrichedGOCC,
                  background = genes$name,
                  org.db = orgdb,
                  id_type = "symbol",
                  min_size = 5,
                  max_size = 200)

print("findEnrichedGOCC finished")

# add ontology columns

gobpList = mclapply(gobpList, 
                    addOntologyCol, 
                    "Biological.Process",
                    mc.cores = detectCores()-1)

print("addOntologyCol gobpList finished")

gomfList = mclapply(gomfList, 
                    addOntologyCol, 
                    "Molecular.Function",
                    mc.cores = detectCores()-1)

print("addOntologyCol gomfList finished")

goccList = mclapply(goccList, 
                    addOntologyCol, 
                    "Cellular.Compartment",
                    mc.cores = detectCores()-1)

print("addOntologyCol goccList finished")

# replace failures
empty =  data.frame(matrix(ncol = 8, nrow = 0))
colnames(empty) = c("GO.ID", 
                    "Term", 
                    "Annotated", 
                    "Significant", 
                    "Expected", 
                    "pval", 
                    "FDR", 
                    "OverlappingGenes")

findFails = function(go){
  if(class(go)!="data.frame"){
    go=empty
  }
  return(go)
}

gobpList = mclapply(gobpList, 
                    findFails, 
                    mc.cores = detectCores()-1)

print("findFails gobpList finished")

gomfList = mclapply(gomfList, 
                    findFails, 
                    mc.cores = detectCores()-1)

print("findFails gomfList finished")

goccList = mclapply(goccList, 
                    findFails, 
                    mc.cores = detectCores()-1)

print("findFails goccList finished")

# bind results from different ontologies
final = list()
for (i in 1:length(files)){
  print(paste0(i))
  final[[i]] = rbind(gobpList[[i]], gomfList[[i]], goccList[[i]])
  final[[i]] = final[[i]] %>% dplyr::filter(Significant > 0)
  file.name = basename(files[i])
  file.name = gsub(ext, ".GO.csv", file.name)
  write.csv(final[[i]], file = paste0(enrichdir, "/", file.name))

}

print("bind results from different ontologies finished")
names(final) = gsub(ext, ".GO", basename(files))
save(final, file = paste0(enrichdir, "/GOenrichment", ext, ".Rdata"))

print("all done!")
