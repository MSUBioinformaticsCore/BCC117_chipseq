#!/bin/bash --login
#SBATCH --mem=256GB
#SBATCH --job-name=DiffBind
#SBATCH --output=%x-%j.out
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --account=bioinformaticscore
#SBATCH --cpus-per-task=16

################################################################## 
# Print parameters
##################################################################

start=`date +%s`
echo "sample sheet:" $1
echo "results dir:" $2
echo "DiffBind.R path:" $3

################################################################## 
# Run DiffBind.R
##################################################################

module purge
module load R-bundle-CRAN/2023.12-foss-2023a

Rscript $3 $1 $2

################################################################## 
# Finish
##################################################################

echo "Finished"
end=`date +%s`
runtime=$((end-start))
echo execution time was `expr $end - $start` 
