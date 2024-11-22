# BCC117_chipseq
 Thomashow lab: ChIP-seq analysis of CAMTA3 at warm and long-term cold conditions

## Set up

### Set paths
```{bash, eval=FALSE}
## set paths
PIPES=/mnt/research/bioinformaticsCore/pipelines
PROJECT=/mnt/research/bioinformaticsCore/projects/thomashowm/BCC117_chipseq
GENOME=/mnt/research/bioinformaticsCore/shared/Genomes
SOFTWARE=/mnt/research/bioinformaticsCore/software
```

### Set project directories
```{bash, eval=FALSE}
mkdir $PROJECT/data
mkdir $PROJECT/results
mkdir $PROJECT/src
mkdir $PROJECT/run
```

### Make temp and run directories
```{bash, eval=FALSE}
mkdir $SCRATCH/thomashowm
mkdir $SCRATCH/thomashowm/BCC117_chipseq
```

### Download reference genome and annotation
```{bash, eval=FALSE}
mkdir $GENOME/arabidopsis_thaliana
mkdir $GENOME/arabidopsis_thaliana/ensembl_release60
cd $GENOME/arabidopsis_thaliana/ensembl_release60

wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz 

wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/gff3/arabidopsis_thaliana/Arabidopsis_thaliana.TAIR10.60.gff3.gz

mv $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.60.gff3.gz $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.60.gff.gz 

gunzip $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.60.gff.gz

module load cufflinks/2.2.1-linux-x86_64
gffread Arabidopsis_thaliana.TAIR10.60.gff -T -o Arabidopsis_thaliana.TAIR10.60.gtf
```

### Install R packages

Install the R packages required for this analysis:

 * [`org.At.tair.db`](https://bioconductor.org/packages/release/data/annotation/html/org.At.tair.db.html)    
 * [`topGO`](https://bioconductor.org/packages/release/bioc/html/topGO.html)   
 * [`AnnotationHub`](https://www.bioconductor.org/packages/release/bioc/html/AnnotationHub.html)   
 * [`GenomicRanges`](https://bioconductor.org/packages/release/bioc/html/GenomicRanges.html)   
 * [`plyranges`](https://www.bioconductor.org/packages/release/bioc/html/plyranges.html)   
 * [`profileplyr`](https://www.bioconductor.org/packages/release/bioc/html/profileplyr.html)   
 * [`DiffBind`](https://bioconductor.org/packages/release/bioc/html/DiffBind.html)   
 * [`tidyverse`](https://www.tidyverse.org)   
 
```{bash}
module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $PROJECT/src/packages.R
```
 
### Download exclusion regions

Exclusion lists identify areas of the genome that have low mappability rates or contain high artifactual signals. Reads mapping to these regions should be masked out before applying ChIP-seq peak-calling software such as MACS3. `getExclusionList.R` uses `excluderanges v 0.99.8` to download the exclusion list for TAIR10 defined by the [Greenscreen](https://academic.oup.com/plcell/article/34/12/4795/6705244) pipeline. Adds Mt and Pt chromosomes to blacklist.

**Script**: `getExclusionList.R`

**Arguments**

1. path to output directory

**Output**

`TAIR10.Klasfeld.arabidopsis_greenscreen_20inputs.bed` 

```{bash, eval=FALSE}
# find Pt and Mt sizes
cd $GENOME/arabidopsis_thaliana/ensembl_release60

gunzip -c Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz > Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
cut -f1,2 Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.fai > sizes.genome

module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $PROJECT/src/getExclusionList.R \
  $GENOME/arabidopsis_thaliana
```

### Make nf-core sample sheet

Sample sheet to specify the input to the nf-core/chipseq pipeline. See [here](https://nf-co.re/chipseq/2.0.0/docs/usage/) for more details. Script works for one antibody. Excludes `warm_input3`, `cold_IP3`, and `cold_input3`. 

**Script**: `makeSampleSheet.R` specific for this analysis   

**Arguments**

1. path to `data` directory containing the `.fastq.gz` files  
2. The ChIPed antibody    

**Output**: `data/sample_sheets/sample_sheet.csv` with the following columns

  * `sample`: Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (_).   
  * `fastq_1`: Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension “.fastq.gz” or “.fq.gz”.   
  * `fastq_2`: Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension “.fastq.gz” or “.fq.gz”.   
  * `antibody`:	Antibody name. This is required to segregate downstream analysis for different antibodies. Required when control is specified.
* `control`: Sample name for control (INPUT) sample.    

```{bash, eval=FALSE}
module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $PROJECT/src/makeSampleSheet.R \
  $PROJECT/data/fastq \
  CAMTA3 
```

## Install nfcore/chipseq pipeline

### Create and activate a new Conda environment with nf-core and nextflow 

```{bash}
module purge
module load Conda/3
which conda 

# if you haven't already
conda create \
  --name nf-core\
  python=3.8 \
  nf-core nextflow 

# activate the environment 
conda activate nf-core
```

### Download the most recent nf-core sigularity container 

Make a folder in the pipelines directory called `nfcore_chipseq`
```{bash}
PIPES=/mnt/research/bioinformaticsCore/pipelines
mkdir $PIPES/nfcore_chipseq
```

Download the most recent nf-core singularity container if you haven't already
```{bash}
cd $PIPES/nfcore_chipseq

nf-core download chipseq \
  -r 2.1.0 \
  -x none \
  -s singularity 
```

## Run the nf-core chipseq pipeline

```{bash}
nextflow run $PIPES/nfcore_chipseq/nf-core-chipseq_2.1.0/2_1_0 \
  --input $PROJECT/data/sample_sheets/sample_sheet.csv \
  --outdir $PROJECT/results \
  --fasta $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa \
  --gtf $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.60.gtf \
  --blacklist $GENOME/arabidopsis_thaliana/ensembl_release60/TAIR10.Klasfeld.arabidopsis_greenscreen_20inputs_MtPt.bed \
  --narrow_peak \
  --read_length 150 \
  -profile singularity \
  -c $PIPES/nfcore_chipseq/src/slurm.config.sh \
  -w $SCRATCH/thomashowm/BCC117_chipseq 
```

## Get chromosome size file

```{bash}
module purge
module load SAMtools/1.19.2-GCC-13.2.0

fasta=$GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa

samtools faidx $fasta

cut -f1,2 ${fasta}.fai > $GENOME/arabidopsis_thaliana/ensembl_release60/sizes.genome
```

## Consensus peaks

Uses [*consensusSeekeR*](https://www.bioconductor.org/packages/release/bioc/html/consensusSeekeR.html) to identify peaks shared between n-1 replicates of a condition. Within each peak, MACS3 reports the location highest signal intensity. *consensusSeekeR* calculates the median location of highest intensity among overlapping peaks. The width of the consensus peak is the merged width of the replicate peaks.

**Script**: `getConsensusPeaks.R`  

**Arguments**

1. path chrom sizes file    
2. path to narrowPeak files    
3. Pattern distinguishing which files belong to which group   

**Output**: `path/to/narrowPeak_files/<pattern>_consensus.peaks.bed` 

```{bash}
module purge
module load R-bundle-CRAN/2023.12-foss-2023a

# warm consensus
Rscript $PROJECT/src/getConsensusPeaks.R \
  $GENOME/arabidopsis_thaliana/ensembl_release60/sizes.genome \
  $PROJECT/results/bwa/merged_library/macs3/narrow_peak \
  warm 

# cold consensus
Rscript $PROJECT/src/getConsensusPeaks.R \
  $GENOME/arabidopsis_thaliana/ensembl_release60/sizes.genome \
  $PROJECT/results/bwa/merged_library/macs3/narrow_peak \
  cold 
```

## Find condition-specific consensus peaks 

```{bash}
module purge
module load BEDTools/2.31.0-GCC-12.3.0
cd $PROJECT/results/bwa/merged_library/macs3/narrow_peak

# cold only
bedtools intersect \
  -a cold_consensus.peaks.bed  \
  -b warm_consensus.peaks.bed \
  -v \
  > cold_only.consensus.bed 

# warm only  
bedtools intersect \
  -a warm_consensus.peaks.bed \
  -b cold_consensus.peaks.bed \
  -v \
  > warm_only.consensus.bed
  
```

## Differential binding analysis

We used [DiffBind v3.8.4](https://bioconductor.org/packages/release/bioc/html/DiffBind.html) to identify peaks differentially bound between each pair of conditions. Specifically, for each pair of conditions, a union peakset was derived that includes the consensus peaks from both conditions. For each replicate of the two conditions, the reads within the union peak set are counted and normalized by the total read depth. Then, the normalized counts from each condition are compared to see if there is a significant difference (FDR < .05) in the number of reads within each peak between conditions. 

This code outputs `.csv` files containing the locations, fold changes, and p-values for differentially bound peaks, volcano plots and binding profile heat maps for significantly gained and lost peaks.

### DiffBind Sample sheet

Make a `.csv` file with the following columns:

1. `SampleID`
1. `Condition`    
2. `Replicate`
3. `bamReads`: path to the IP `.bam` file
4. `bamControl`: path to the input `.bam` file
5. `Peaks`: path to the consensus peak `.bed` files
6. `PeakCaller`: bed

**Script**: `makeDiffBindSampleSheet.R`, specific to this project   

**Arguments** 

1. Path to the directory containing the `.bam` files   
2. Path to the directory containing the consensus peak `.bed` files   
3. Path to the data directory     

**Output**

`data/sample_sheets/DiffBind_sample_sheet_consensus.csv`
`data/sample_sheets/DiffBind_sample_sheet_narrowPeak.csv`

```{bash}
module purge
module load R-bundle-CRAN/2023.12-foss-2023a

# specific for this project
Rscript $PROJECT/src/makeDiffBindSampleSheet.R \
  $PROJECT/results/bwa/merged_library \
  $PROJECT/results/bwa/merged_library/macs3/narrow_peak \
  $PROJECT/data

```

### Make a custom Greylist

From the Diffbind tutorial:

> Greylists are specific to a ChIP-seq experiment, and are derived from the controls generated as part of the experiment. The idea is to analyze libraries that are not meant to show systematic enrichment (such as Inputs, in which no anti-body is introduced), and identify anomalous regions where a disproportionate degree of signal is present. These regions can then be excluded from subsequent analysis.
>  Application of greylists prevents identification of problematic genomic regions in the materials used in the experiment as being differentially bound. For example, these could include areas of high copy-number alterations in a cell line.
>  Prior to version 3.0, the default in DiffBind has been to simply subtract control reads from ChIP reads in order to dampen the magnitude of enrichment in anomalous regions. Greylists represent a more principled way of accomplishing this. If a greylist has been applied, the current default in DiffBind is to not subtract control reads.

### Make greylist

**Script:** `run_customGreylist.sh`

**Arguments**

1. Path to the DiffBind sample sheet 
2. Path to the sizes.genome file (tab-delimited file with chromosome sizes)   
3. Output directory for the results   
4. Path to customGreylist.R    

**Output**

`greylist.Rds`: An RDS file containing the greylist information.

```{bash}
cd $PROJECT/run

sbatch $PROJECT/src/run_customGreylist.sh \
  $PROJECT/data/sample_sheets/DiffBind_sample_sheet_consensus.csv \
  $GENOME/arabidopsis_thaliana/ensembl_release60/sizes.genome \
  $PROJECT/results/bwa/merged_library \
  $PROJECT/src/customGreylist.R 
```

### Run `runDiffBind.sh`

**Arguments**

1. Path to the DiffBind sample sheet    
2. Output directory for the results    
3. Baseline condition   
4. Path to `greylist.Rds`   
5. Path to `DiffBind.R`   

**Output**

* new `diffbind` directory to hold the results
* `diffbind/*.csv` with the DiffBind results
* `diffbind/*_volcano_plot.pdf` 
* `diffbind/*_profile_plot.pdf`      

```{bash, eval=F}

cd $PROJECT/run

sbatch $PROJECT/src/runDiffBind.sh \
  $PROJECT/data/sample_sheets/DiffBind_sample_sheet_consensus.csv \
  $PROJECT/results \
  "warm" \
  $PROJECT/results/bwa/merged_library/greylist.Rds \
  $PROJECT/src/DiffBind.R

```

## Peak annotation with HOMER

See HOMER's [installation instructions](http://homer.ucsd.edu/homer/introduction/install.html) for details. 
```{bash, eval=F}
module purge
module load SAMtools/1.19.2-GCC-13.2.0
module load BLAT/3.7-GCC-12.3.0

#Add HOMER to bashrc if you haven't
PATH=$PATH:/mnt/research/bioinformaticsCore/software/HOMER/.//bin/

# make an annotation dir
mkdir $PROJECT/results/annotatePeaks
```

Install the genome/annotation package for the genome of interest
```{bash,eval=FALSE}
#to see available organisms
perl /mnt/research/bioinformaticsCore/software/HOMER/configureHomer.pl -list
#to install organisms 
perl /mnt/research/bioinformaticsCore/software/HOMER/configureHomer.pl -install tair10
```

### Annotate peaks

**Script** `$PROJECT/src/annotatePeaks.sh`

**Output**

  * The `annotatePeaks/*.anno` annotation files generated by HOMER    
  * The `annotatePeaks/*.annStats` files listing the enrichment of peaks in types of genomic regions like promoters, exons, and introns.   

**Arguments**

  1. Path to the directory containing peak files in bed format   
  2. Path to output directory     
  3. Pattern in the file names you want to annotate   
  4. The genome/annotation package for your genome of interest
  5. Homer motif file if you have one   
  6. The gft file you want the annotations to come from     

```{bash,eval=FALSE}
cd $PROJECT/run

# consensus peaks
sbatch $PROJECT/src/annotatePeaks.sh \
  $PROJECT/results/bwa/merged_library/macs3/narrow_peak \
  $PROJECT/results/annotatePeaks \
  _consensus.peaks.bed \
  tair10 \
  NONE \
  $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.60.gtf 
  
# condition-specific peaks
sbatch $PROJECT/src/annotatePeaks.sh \
  $PROJECT/results/bwa/merged_library/macs3/narrow_peak \
  $PROJECT/results/annotatePeaks \
  _only.consensus.bed \
  tair10 \
  NONE \
  $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.60.gtf 
  
# diff peaks
sbatch $PROJECT/src/annotatePeaks.sh \
  $PROJECT/results/diffbind \
  $PROJECT/results/annotatePeaks \
  .bed \
  tair10 \
  NONE \
  $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.60.gtf

```

## Spilt DiffBind annation files into LOSS/GAIN

Make annotation and bed files containing only the gained or lost peaks from the DiffBind analyses. These will be used for gene set enrichment. 

**Script**: `splitDiffBind.R`

**Arguments**

1. path to diffbind files   
2. path to annotation files   

**Output**

 * `annotatePeaks/*_gain.anno` and `annotatePeaks/*_loss.anno` files for gene set enrichment  
 * `diffbind/*_gain.bed` and `diffbind/*_loss.bed` files for finding enriched motifs    
 * Gene names will be added to the `diffbind/*.csv` results file    

```{bash}
module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $PROJECT/src/splitDiffBind.R \
  $PROJECT/results/diffbind \
  $PROJECT/results/annotatePeaks
```

## Gene Ontology enrichment

`getGOenrichment.R` uses [`topGO v2.54.0`](https://bioconductor.org/packages/release/bioc/html/topGO.html) with Fisher's exact test and the gene ontology annotations from a species-specific bioconductor annotation package to check for enrichment of genes related to specific biological processes, molecular functions, and cellular compartments among genes near ChIP peaks. Only GO terms with 5 > n genes > 200 are included. GO terms containing several hundred genes are often overly broad and, therefore, uninformative. 

**Script**: `getGOenrichment.R`

**Arguments**

1. path to `ChIPSeqFunctions.R`   
2. input directory    
3. output directory   
4. extension of `.anno` files you want to query   
5. org.eg.db from bioconductor    
6. `*.descriptionfile` containing all genes used by `annotatePeaks.pl` (should be in `HOMER/data`)

**Output**

`geneset_enrichment/*GO.csv` files with the enrichment results

```{bash, eval=FALSE}

module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $PROJECT/src/getGOenrichment.R \
  $PROJECT/src/ChIPSeqFunctions.R \
  $PROJECT/results/annotatePeaks \
  $PROJECT/results \
  .anno \
  org.At.tair.db \
  $SOFTWARE/HOMER/data/accession/arabidopsis.description 
  
```

## Known and Novel Motif enrichment

Uses `meme` and `sea` to find *de novo* motifs and enrichment for known motifs.

**Script**: `findEnrichedMotifs.sh`

**Arguments**

1. Directory containing the peak files
2. Output directory
3. File extension of the peak files ie `.narrowPeak`
4. Path to motif database
5. Path to the reference genome

**Output**

* `motif_enrichment/*_denovo` directory for each peak set
  + `meme.html` - an HTML file that provides the results in an interactive, human-readable format    
  + `meme.txt` - a plain text file of the results for backwards compatibility with earlier versions of MEME   
  + `meme.xml` - an XML file that provides the results in a format designed for machine processing    
  + `logoN.png,.eps` - PNG and EPS images files containing sequence logos for each of the motifs found by MEME (where N is the motif number, logo1.png is shown above)   
* `motif_enrichment/*_known` directory for each peak set
  + `sea.html` - an HTML file that provides the results in an interactive, human-readable format   
  + `sea.tsv` - a TSV (tab-separated values) file that provides the results in a format suitable for parsing by scripts and viewing with Excel    
  + `sequences.tsv` - a TSV (tab-separated values) file that lists the true- and false-positive sequences identified by SEA    

```{bash, eval=F}
cd $PROJECT/run

MOTIFDB=/mnt/research/bioinformaticsCore/shared/motif_databases/20241006_MEME_motifs/motif_databases/ARABD/ArabidopsisDAPv1.meme

sbatch $PROJECT/src/findEnrichedMotifs.sh \
  $PROJECT/results/diffbind \
  $PROJECT/results \
  .bed \
  $MOTIFDB \
  $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
  
sbatch $PROJECT/src/findEnrichedMotifs.sh \
  $PROJECT/results/bwa/merged_library/macs3/narrow_peak \
  $PROJECT/results \
  _consensus.peaks.bed \
  $MOTIFDB \
  $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
  
sbatch $PROJECT/src/findEnrichedMotifs.sh \
  $PROJECT/results/bwa/merged_library/macs3/narrow_peak \
  $PROJECT/results \
  _only.consensus.bed \
  $MOTIFDB \
  $GENOME/arabidopsis_thaliana/ensembl_release60/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa
```

```

