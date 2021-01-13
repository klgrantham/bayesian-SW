# Combine performance_measures files across trial configurations
# into a single data frame for formatting into a results table
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

library(tidyverse)

file.names <- dir('./performance_measures')

allbiasrows <- data.frame()
allRMSErows <- data.frame()
allcovrows <- data.frame()
for (i in 1:length(file.names)) {
  
  load(file.names[i])
  
  params <- results$params
  measures <- results$measures
  
  # Bias values
  biasvals <- measures %>%
    filter(measure %in% c('MCMC_mean_bias', 'REML_bias')) %>%
    mutate(method=sapply(strsplit(measure, '_'), '[', 1)) %>% 
    select(c(theta, WPICC, BPICC, method))
  biasrows <- cbind(params, biasvals)
  allbiasrows <- rbind(allbiasrows, biasrows)
  
  # RMSE values
  RMSEvals <- measures %>%
    filter(measure %in% c('MCMC_mean_mse', 'REML_mse')) %>%
    mutate(method=sapply(strsplit(measure, '_'), '[', 1)) %>%
    select(c(theta, WPICC, BPICC, method))
  RMSErows <- cbind(params, RMSEvals)
  allRMSErows <- rbind(allRMSErows, RMSErows)
  
  # Coverage values
  covvals <- measures %>%
    filter(measure %in% c('MCMC_coverage', 'REML_coverage')) %>%
    mutate(method=sapply(strsplit(measure, '_'), '[', 1)) %>%
    select(c(theta, method))
  covrows <- cbind(params, covvals)
  allcovrows <- rbind(allcovrows, covrows)
  
  # TODO: Include empirical and model-based SE
}

# Convert results to long format for formatting into a table
biasdf <- allbiasrows %>%
  pivot_longer(
    cols=c(theta, WPICC, BPICC),
    names_to='parameters',
    values_to='bias'
  )

RMSEdf <- allRMSErows %>%
  pivot_longer(
    cols=c(theta, WPICC, BPICC),
    names_to='parameters',
    values_to='RMSE'
  )
