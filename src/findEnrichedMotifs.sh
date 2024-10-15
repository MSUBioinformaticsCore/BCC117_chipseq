#!/bin/bash -login
#SBATCH --mem=64GB
#SBATCH --job-name=findEnrichedMotifs
#SBATCH --output=%x-%j.out
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --account=bioinformaticscore
#SBATCH --cpus-per-task=32

################################################################## 
# Print parameters
##################################################################

start=`date +%s`
echo "input dir:" $1
echo "output dir:" $2
echo "file extension:" $3
echo "motif database" $5
echo "genome" $6

indir=$1
outdir=$2
ext=$3
motif=$4
genome=$5

mkdir $outdir/motif_enrichment

################################################################## 
# bed to fasta
##################################################################

module purge 
module load BEDTools/2.31.0-GCC-12.3.0

cd $indir

for bed in *$ext
do

bedtools getfasta -fi $genome \
  -bed $bed \
  -name \
  > $outdir/motif_enrichment/$bed.fa

#remove duplicate sequences  
#https://stackoverflow.com/questions/61374573/how-can-i-eliminate-duplicated-sequences-in-fasta-file  
awk 'BEGIN {i = 1;} { if ($1 ~ /^>/) { tmp = h[i]; h[i] = $1; } else if (!a[$1]) { s[i] = $1; a[$1] = "1"; i++; } else { h[i] = tmp; } } END { for (j = 1; j < i; j++) { print h[j]; print s[j]; } }' < $outdir/motif_enrichment/$bed.fa > $outdir/motif_enrichment/$bed.nodup.fa

done

################################################################## 
# Look for enrichment of known motifs
##################################################################

module purge
module load MEME/5.5.4-gompi-2022b

cd $outdir/motif_enrichment

for fa in *${ext}.nodup.fa
do
sea --p $fa --m $motif --oc ${fa%.fa}_known_motifs
done

################################################################## 
# Look for de novo motifs
##################################################################

for fa in *${ext}.fa
do
meme $fa -nmotifs 10 -evt 0.01 -p 16 -oc ${fa%.fa}_denovo_motifs --use-hwthread-cpus
done

################################################################## 
# Finish
##################################################################

echo "Finished"
end=`date +%s`
runtime=$((end-start))
echo execution time was `expr $end - $start` 