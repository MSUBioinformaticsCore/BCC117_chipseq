#' @param args[1] path chrom sizes file
#' @param args[2] path to narrowPeak files 
#' @param args[3] pattern

args <- commandArgs(TRUE)
library(consensusSeekeR)
# Chromosomes information

chrom_file = args[1]
# chrom_file = "/mnt/research/bioinformaticsCore/shared/Genomes/arabidopsis_thaliana/ensembl_release60/sizes.genome"
chromsizes = read.delim(args[1])
# chromsizes = read.delim("/mnt/research/bioinformaticsCore/shared/Genomes/arabidopsis_thaliana/ensembl_release60/sizes.genome", header=F)

chrInfo <- Seqinfo(seqnames=chromsizes$V1,
                   seqlengths=chromsizes$V2, isCircular=rep(FALSE, nrow(chromsizes)),
                   genome=chrom_file)
narrow_peak_path = args[2]
# narrow_peak_path = "/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/merged_library/macs3/narrow_peak"
narrow_peak_files = list.files(narrow_peak_path, 
                               pattern = ".narrowPeak",
                               full.names = TRUE)
pattern = args[3]
# pattern = "cold"
narrow_peak_files = narrow_peak_files[grepl(pattern, narrow_peak_files)]

peak_names = basename(narrow_peak_files)

narrowPeak_list = list()
peak_list = list()
  
for(file in narrow_peak_files){
  
  REP = basename(file)
  result = readNarrowPeakFile(file, extractRegions = TRUE, extractPeaks = TRUE)
  narrowPeak_list[[REP]] = result$narrowPeak
  peak_list[[REP]] = result$peak
  
  names(narrowPeak_list[[REP]]) = rep(REP, length(narrowPeak_list[[REP]]))
  names(peak_list[[REP]]) = rep(REP, length(peak_list[[REP]]))

}

narrowPeak_grl = GRangesList(narrowPeak_list)
narrowPeak_unlist = unlist(narrowPeak_grl)
peak_grl = GRangesList(peak_list)
peak_unlist =  unlist(peak_grl)

results <- findConsensusPeakRegions(
  narrowPeaks = narrowPeak_unlist,
  peaks = peak_unlist,
  chrInfo = chrInfo,
  extendingSize = 25,
  expandToFitPeakRegion = TRUE,
  shrinkToFitPeakRegion = TRUE,
  minNbrExp = length(narrow_peak_files)-1,
  nbrThreads = 16)

library(GenomicRanges)
library(rtracklayer)
export.bed(format = "bed", 
           object=results$consensusRanges, 
           con = paste0(narrow_peak_path, "/", pattern, "_consensus.peaks.bed"))
