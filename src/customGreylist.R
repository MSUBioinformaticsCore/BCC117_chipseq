# https://www.bioconductor.org/packages/release/bioc/vignettes/GreyListChIP/inst/doc/GreyList-demo.pdf
# https://bioconductor.org/packages/release/bioc/vignettes/DiffBind/inst/doc/DiffBind.pdf
# https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html#looping-over-grangeslist-objects
#' @param args[1] Path to the DiffBind sample sheet (CSV file with sample details)
#' @param args[2] Path to the chromsizes file (tab-delimited file with chromosome sizes)
#' @param args[3] Output directory for the results

# Parse command-line arguments
args <- commandArgs(TRUE)

# Load the required library for GreyListChIP analysis
library(GreyListChIP)

# Set a random seed for reproducibility
set.seed(1)

# Read the sample sheet (provided as the first argument)
sample_sheet = read.csv(args[1])

# Read the chromsizes file (provided as the second argument)
chromsizes = read.delim(args[2])

# Set the output directory (provided as the third argument)
outdir = args[3]

# Example paths for testing purposes (commented out for production use)
# sample_sheet = read.csv("/path/to/DiffBind_sample_sheet_consensus.csv")
# chromsizes = "/path/to/sizes.genome"
# outdir = "/path/to/results_directory"

# Initialize an empty list to store greylist regions for each sample
gl_list = list()

# Construct names for the control input samples based on the sample sheet
inputs = paste0(sample_sheet$Condition, 
                "_REP",
                sample_sheet$Replicate,
                "_input")

# Loop through each row in the sample sheet to process control BAM files
for(i in 1:nrow(sample_sheet)){
  
  # Create a new GreyList object with the chromsizes file
  gl <- new("GreyList", karyoFile = chromsizes)
  
  # Count reads from the control BAM file specified in the sample sheet
  gl <- countReads(gl, sample_sheet$bamControl[i])
  
  # Calculate the threshold for greylist regions using multiple cores
  gl <- calcThreshold(gl, cores = 16)
  
  # Generate the greylist regions
  gl <- makeGreyList(gl)
  
  # Store the greylist regions in the list
  gl_list[[i]] = gl@regions
}

# Assign sample names to the greylist regions in the list
names(gl_list) = inputs

# Combine all greylist regions into a GRangesList object
grl = GRangesList(gl_list)

# Unlist all regions to get a single set of regions
all = unlist(grl)

# Reduce overlapping regions to create a master greylist
master = reduce(all)

# Create a final object containing the master greylist and individual control lists
final = list(master = master, controls = grl)

# Save the final greylist object to an RDS file in the output directory
saveRDS(final, file = paste0(outdir, "/greylist.Rds"))



