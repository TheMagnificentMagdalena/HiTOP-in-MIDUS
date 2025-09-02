#!/bin/bash -l        
#SBATCH --time=48:00:00
#SBATCH --ntasks=5
#SBATCH --mem=24g
#SBATCH --tmp=10g
#SBATCH --mail-type=BEGIN,END,FAIL  
#SBATCH --mail-user=march341@umn.edu
#SBATCH --job-name=imputation_job
#SBATCH --output=imputation_output_%j.txt

cd /users/2/march341/multiple_imputation

module load R/4.4.2-openblas-rocky8

Rscript Iterative_Efa_Velicer_Map.R
