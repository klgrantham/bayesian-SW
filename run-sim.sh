#!/bin/env bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8192
#SBATCH --cpus-per-task=1
#SBATCH --time=5-00:00:00
#SBATCH --job-name=sim-study
#SBATCH --array=1-156
#SBATCH --output=./output/%A_%a.out

module load gcc/8.1.0
module load R/4.0.0-openblas

Rscript run-sim.R
