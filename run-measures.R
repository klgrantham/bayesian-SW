#!/bin/env Rscript
#
# Collates results files and calculates performance measures
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('performance-measures.R')

v <- Sys.getenv('SLURM_ARRAY_TASK_ID', NA)
vn <- as.integer(v)
row <- vn + 1
print(sprintf('We will run job #%d',row))
params <- read.csv('parameters.csv')
myparams = params[row,]
print('We will use the params:')
print(myparams)

reduce_all_results(2, myparams$clust_per_seq, myparams$periods,
                   myparams$subjects, myparams$WPICC,
                   myparams$CAC, myparams$theta)

collate_results(2, myparams$clust_per_seq, myparams$periods,
                myparams$subjects, myparams$WPICC,
                myparams$CAC, myparams$theta)

calculate_measures(2, myparams$clust_per_seq, myparams$periods,
                   myparams$subjects, myparams$WPICC,
                   myparams$CAC, myparams$theta)
