#!/bin/sh -login
#SBATCH --mem=64GB
#SBATCH --job-name=GOenrichment
#SBATCH --output=%x-%j.out
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --account=bioinformaticscore
#SBATCH --cpus-per-task=16

################################################################## 
# Print parameters
##################################################################

start=`date +%s`
echo "ChIPSeqFunctions.R path:" $1
echo "input dir:" $2
echo "output dir:" $3
echo "file extension:" $4
echo "orgdb from bioconductor:" $5
echo "Homer gene file:" $6
echo "getGOenrichment.R path:" $7

################################################################## 
# Run getGOenrichment.R
##################################################################

module purge
module load GCC/8.3.0  OpenMPI/3.1.4  R/4.1.0

Rscript $7 $1 $2 $3 $4 $5 $6

################################################################## 
# Finish
##################################################################

echo "Finished"
end=`date +%s`
runtime=$((end-start))
echo execution time was `expr $end - $start` 