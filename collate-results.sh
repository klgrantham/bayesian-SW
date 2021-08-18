#!/bin/env bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=4096
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --job-name=collate-results
#SBATCH --array=1-48
#SBATCH --output=./output/%A_%a.out

module load gcc/8.1.0
module load R/4.0.0-openblas

Rscript run-measures.R
