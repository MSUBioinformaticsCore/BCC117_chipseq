#' @param args[1] path chrom sizes file
#' @param args[2] path to narrowPeak files 
#' @param args[3] pattern

args <- commandArgs(TRUE) # Capture command line arguments
library(consensusSeekeR) # Load the consensusSeekeR library

# Chromosomes information
chrom_file = args[1] # Path to chromosome sizes file
# chrom_file = "/mnt/research/bioinformaticsCore/shared/Genomes/arabidopsis_thaliana/ensembl_release60/sizes.genome"
chromsizes = read.delim(args[1], header=F) # Read chromosome sizes file
# chromsizes = read.delim("/mnt/research/bioinformaticsCore/shared/Genomes/arabidopsis_thaliana/ensembl_release60/sizes.genome", header=F)

# Create chromosome information object
chrInfo <- Seqinfo(seqnames=chromsizes$V1,
                   seqlengths=chromsizes$V2, 
                   isCircular=rep(FALSE, nrow(chromsizes)),
                   genome=chrom_file)

narrow_peak_path = args[2] # Path to narrowPeak files
# narrow_peak_path = "/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/merged_library/macs3/narrow_peak"

# List all narrowPeak files in the directory
narrow_peak_files = list.files(narrow_peak_path, 
                               pattern = ".narrowPeak$",
                               full.names = TRUE)

pattern = args[3] # Pattern to filter narrowPeak files
# pattern = "cold"
narrow_peak_files = narrow_peak_files[grepl(pattern, narrow_peak_files)] # Filter files based on pattern

narrowPeak_list = list() # Initialize list for narrowPeak data
peak_list = list() # Initialize list for peak data

# Process each narrowPeak file
for(file in narrow_peak_files){
  
  REP = basename(file) # Get base name of the file
  result = readNarrowPeakFile(file, extractRegions = TRUE, extractPeaks = TRUE) # Read narrowPeak file
  narrowPeak_list[[REP]] = result$narrowPeak # Store narrowPeak regions
  peak_list[[REP]] = result$peak # Store peaks
  
  # Assign file names to narrowPeak and peak lists
  names(narrowPeak_list[[REP]]) = rep(REP, length(narrowPeak_list[[REP]]))
  names(peak_list[[REP]]) = rep(REP, length(peak_list[[REP]]))
}

# Create GRangesList objects
narrowPeak_grl = GRangesList(narrowPeak_list)
narrowPeak_unlist = unlist(narrowPeak_grl) # Unlist narrowPeak regions
peak_grl = GRangesList(peak_list)
peak_unlist =  unlist(peak_grl) # Unlist peaks

median_width = median(width(narrowPeak_unlist)) # Calculate median width of narrowPeak regions

# Find consensus peak regions
results <- findConsensusPeakRegions(
  narrowPeaks = narrowPeak_unlist,
  peaks = peak_unlist,
  extendingSize = round(median_width/2), # Extend size based on median width
  chrInfo = chrInfo,
  expandToFitPeakRegion = TRUE, # Expand regions to fit peaks
  shrinkToFitPeakRegion = TRUE, # Shrink regions to fit peaks
  minNbrExp = length(narrow_peak_files)-1, # Minimum number of experiments
  nbrThreads = 16) # Number of threads to use

# Assign consensus IDs to results
consensusIDS = as.data.frame(paste0(pattern, "_consensusPeak_", 1:length(results$consensusRanges)))
colnames(consensusIDS) = "consensusID"
mcols(results$consensusRanges) <- consensusIDS

library(GenomicRanges) # Load GenomicRanges library

# Find overlaps between narrowPeak regions and consensus ranges
hits = as.data.frame(findOverlaps(narrowPeak_unlist, results$consensusRanges))
consensus_hits = as.data.frame(results$consensusRanges[hits$subjectHits]) # Extract consensus ranges
narrowPeak_hits = as.data.frame(narrowPeak_unlist[hits$queryHits], row.names = NULL) # Extract narrowPeak regions

# Rename columns for clarity
colnames(consensus_hits)[1:5] = paste0("consensus_", colnames(consensus_hits)[1:5])
peak2consensus = cbind(narrowPeak_hits, consensus_hits) # Combine narrowPeak and consensus data

# Write detailed peak-to-consensus mapping to a CSV file
write.csv(peak2consensus, file = paste0(narrow_peak_path, "/", pattern, "_consensus.peaks_info.csv"))

library(tidyverse) # Load tidyverse for data manipulation

# Make the consensus score qValue of the constituent peak with the highest qValue
consensus_score = 
  peak2consensus %>%
  group_by(consensusID) %>%
  slice_max(qValue, n = 1) %>% # Select peak with highest qValue
  select(consensus_seqnames, consensus_start, 
         consensus_end, consensusID, qValue, strand) %>%
  mutate(strand = ".") # Set strand information

# Write consensus peaks to a BED file
write.table(consensus_score, 
            file = paste0(narrow_peak_path, "/", pattern, "_consensus.peaks.bed"),
            quote = FALSE,
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE)
