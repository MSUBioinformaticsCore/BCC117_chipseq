# ### DiffBind Sample sheet
# 
# This is a `.csv` file with the following columns:
#   
# 1. `SampleID`
# 2. `Condition`    
# 2. `Replicate`
# 4. `bamReads`: path to the IP `.bam` file
# 5. `bamControl`: path to the input `.bam` file
# 6. `Peaks`: path to the consensus peak `.bed` files
# 7. `PeakCaller`: bed
# 
# See `sample_sheets/DiffBind_samplesheet_consensus` for an example 

bam_path = "/mnt/ufs18/rs-013/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/mergedLibrary"
peak_path = "/mnt/ufs18/rs-013/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/results/bwa/mergedLibrary/macs2/narrowPeak"

ip = list.files(bam_path, 
                  pattern = ".+IP.+.bam$", 
                  full.names = T)

input = list.files(bam_path, 
                pattern = ".+input.+.bam$", 
                full.names = T)

consensus = list.files(peak_path, 
                       pattern = "_all.bed", 
                       full.names = T)

sheet = data.frame(bamReads = ip,
                   bamControl = input,
                   Peaks = c(rep(consensus[1],4),
                             rep(consensus[2],3)),
                   PeakCaller = "bed")
sheet = 
  sheet %>%
  mutate(SampleID=gsub(".mLb.clN.sorted.bam", 
                       "", basename(bamReads)),
         Condition = gsub("_.*", "", SampleID),
         Replicate = gsub(".*_IP", "", SampleID)) %>%
  select(SampleID, Condition, Replicate,
         bamReads, bamControl, Peaks, PeakCaller)
 
write.csv(sheet, file = "/mnt/ufs18/rs-013/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/data/sample_sheets/DiffBind_sample_sheet.csv", row.names = F) 
