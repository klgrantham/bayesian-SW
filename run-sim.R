#!/bin/env Rscript
#
# Extracts one set of input parameters from spreadsheet, then fits both models
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

v <- Sys.getenv('SLURM_ARRAY_TASK_ID', NA)
vn <- as.integer(v)
print(sprintf('We will run job #%d', vn))
params <- read.csv('parameters_batched.csv')
myparams <- params[params$task_id==vn,]
print('We will use the params:')
print(myparams)

start.time <- Sys.time()
source('sim.R')
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste('Source sim.R script took:', format(time.taken, digits=4)))

start.time <- Sys.time()
if (myparams$CAC==1.0) {
  stanmod <- stan_model('HH-ncp.stan')
} else {
  stanmod <- stan_model('HG-ncp.stan')
}
end.time <- Sys.time()
time.taken <- end.time - start.time
print(paste('Compile Stan took:', format(time.taken, digits=4)))

start.time.reps <- Sys.time()
for (i in myparams$start:myparams$end) {
  start.time <- Sys.time()
  fit_models(myparams, stanmod, i)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('Replicate runtime:', format(time.taken, digits=4)))
}
end.time.reps <- Sys.time()
time.taken.reps <- end.time.reps - start.time.reps
print(paste('Total runtime:', format(time.taken.reps, digits=4)))
