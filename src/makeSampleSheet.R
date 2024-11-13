#' @param args[1] data dir
#' @param args[2] antibody

# make sample sheet for nfcore-chipseq
library(tidyverse)

#data_dir = "/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/data/fastq"
data_dir = args[1]
ab = args[2]

R1 = list.files(data_dir, 
                pattern = "R1_001.fastq.gz", 
                full.names = TRUE)

R2 = list.files(data_dir, 
                pattern = "R2_001.fastq.gz", 
                full.names = TRUE)

samp = data.frame(fastq_1 = R1,
                  fastq_2 = R2)

samp$sample = gsub("_S.*", "", basename(samp$fastq_1))

# replicates have to be sequential so rep 3 is actually rep 4
samp = 
  samp %>%
  filter(sample  != "warm_input3") %>%
  mutate(control = case_when(
    grepl("IP", sample) ~ gsub("IP", "input", sample),
    grepl("input", sample) ~ "")) %>%
  mutate(antibody = case_when(
    grepl("IP", sample) ~ ab,
    grepl("input", sample) ~ "")) %>%
  mutate(replicate = c(rep(c(1,2,3,4), 2), rep(c(1,2,3), 2)),
         control_replicate = case_when(
           grepl("input", sample) ~ "",
           grepl("IP", sample) ~ as.character(replicate))) %>%
  select(sample, fastq_1, fastq_2, replicate,
         antibody, control, control_replicate)

samp$sample = str_sub(samp$sample, start = 1, end = nchar(samp$sample)-1)
samp$control = str_sub(samp$control, start = 1, end = nchar(samp$control)-1)

sample_dir = paste0(dirname(data_dir), "/sample_sheets")
if(!dir.exists(sample_dir)){dir.create(sample_dir)}

write.table(samp, file = paste0(sample_dir, "/sample_sheet.csv"), 
            sep = ",", 
            quote=F,
            row.names = F)
