#!/bin/bash -login
#SBATCH --mem=64GB
#SBATCH --job-name=annotatePeaks
#SBATCH --output=%x-%j.out
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --account=bioinformaticscore

################################################################## 
# Print parameters
##################################################################

echo "datadir:" $1
echo "outdir:" $2
echo "extension:" $3
echo "genome:" $4
echo "motif file:" $5
echo "gtf:" $6

################################################################## 
# Set up
##################################################################

start=`date +%s`
datadir=$1
outdir=$2
ext=$3
genome=$4
motif=$5
gtf=$6

################################################################## 
# annotate peaks
##################################################################

module purge
module load SAMtools/1.19.2-GCC-13.2.0

cd $datadir

for bed in *$ext
do

  awk '{gsub("\t.\t","\t+\t",$0); print;}' $bed > ${bed}.homer
  
  if [ $motif != "NONE" ]
  then
  
  annotatePeaks.pl ${bed}.homer \
    $genome \
    -m $motif \
    -annStats ${outdir}/${bed}.annStats \
    -gtf $gtf \
    > ${outdir}/${bed}.anno
    
  else
  
   annotatePeaks.pl ${bed}.homer \
    $genome \
    -annStats ${outdir}/${bed}.annStats \
    -gtf $gtf \
    > ${outdir}/${bed}.anno
    
  fi
done

rm *.homer
