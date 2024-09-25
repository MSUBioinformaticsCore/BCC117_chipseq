# make sample sheet for nfcore-chipseq
library(tidyverse)

data_dir = "/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq/data"

R1 = list.files(data_dir, 
                pattern = "R1_001.fastq.gz", 
                full.names = TRUE)

R2 = list.files(data_dir, 
                pattern = "R2_001.fastq.gz", 
                full.names = TRUE)

samp = data.frame(fastq_1 = R1,
                  fastq_2 = R2)

samp$sample = gsub("_S.*", "", basename(samp$fastq_1))

samp = 
  samp %>%
  mutate(control = case_when(
    grepl("IP", sample) ~ gsub("IP", "input", sample),
    grepl("input", sample) ~ "")) %>%
  mutate(antibody = case_when(
    grepl("IP", sample) ~ "CAMTA3",
    grepl("input", sample) ~ "")) %>%
  select(sample, fastq_1, fastq_2, antibody, control) %>%
  filter(sample != "warm_input3")

write.table(samp, file = paste0(data_dir, "/sample_sheet.csv"), 
            sep = ",", 
            quote=F,
            row.names = F)
