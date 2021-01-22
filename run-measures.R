#!/bin/env Rscript
#
# Collates results files and calculates performance measures
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

source('performance-measures.R')

v <- Sys.getenv('SLURM_ARRAY_TASK_ID', NA)
vn <- as.integer(v)
print(sprintf('We will run job #%d', vn))
params <- read.csv('parameters_batched.csv')
# Get unique specification with spec_id
myparams <- params[params$spec_id==vn,][1,]
print('We will use the params:')
print(myparams)

reduce_all_results(10, myparams$clust_per_seq, myparams$periods,
                   myparams$subjects, myparams$WPICC,
                   myparams$CAC, myparams$theta)

collate_results(10, myparams$clust_per_seq, myparams$periods,
                myparams$subjects, myparams$WPICC,
                myparams$CAC, myparams$theta)

calculate_measures(myparams$clust_per_seq, myparams$periods,
                   myparams$subjects, myparams$WPICC,
                   myparams$CAC, myparams$theta)
