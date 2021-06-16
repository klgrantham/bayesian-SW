# Combine performance_measures files across trial configurations
# into data frames for formatting into results tables
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

library(tidyverse)

file.names <- dir('./performance_measures')

allreprows <- data.frame()
allbiasrowsHH <- data.frame()
allMSErowsHH <- data.frame()
allcovrowsHH <- data.frame()
allintlenrowsHH <- data.frame()
allempSErowsHH <- data.frame()
allmodSErowsHH <- data.frame()
allpcterrmodSErowsHH <- data.frame()
allbiasrowsHG <- data.frame()
allMSErowsHG <- data.frame()
allcovrowsHG <- data.frame()
allintlenrowsHG <- data.frame()
allempSErowsHG <- data.frame()
allmodSErowsHG <- data.frame()
allpcterrmodSErowsHG <- data.frame()
for (i in 1:length(file.names)) {
  
  load(paste0('performance_measures/', file.names[i]))
  
  params <- results$params
  measures <- results$measures
  reps <- results$reps
  
  reprow <- cbind(params, reps)
  allreprows <- rbind(allreprows, reprow)
  
  if (params$r==1.0) {
    pars <- c('theta', 'WPICC')
    
    # Bias values
    biasvals <- measures %>%
      filter(measure %in% c('bias', 'MCSE_bias')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(all_of(pars), measure, Method))
    biasrows <- cbind(params, biasvals)
    allbiasrowsHH <- rbind(allbiasrowsHH, biasrows)
    
    # MSE values
    MSEvals <- measures %>%
      filter(measure %in% c('MSE', 'MCSE_MSE')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(all_of(pars), measure, Method))
    MSErows <- cbind(params, MSEvals)
    allMSErowsHH <- rbind(allMSErowsHH, MSErows)
    
    # Coverage values
    covvals <- measures %>%
      filter(measure %in% c('coverage', 'MCSE_coverage')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(theta, measure, Method))
    covrows <- cbind(params, covvals)
    allcovrowsHH <- rbind(allcovrowsHH, covrows)
    
    # Empirical SE
    empSEvals <- measures %>%
      filter(measure %in% c('empSE', 'MCSE_empSE')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(all_of(pars), measure, Method))
    empSErows <- cbind(params, empSEvals)
    allempSErowsHH <- rbind(allempSErowsHH, empSErows)
    
    # Average model-based SE
    modSEvals <- measures %>%
      filter(measure %in% c('avgmodSE', 'MCSE_avgmodSE')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(theta, measure, Method))
    modSErows <- cbind(params, modSEvals)
    allmodSErowsHH <- rbind(allmodSErowsHH, modSErows)
    
    # Relative % error in model-based SE
    pcterrmodSEvals <- measures %>%
      filter(measure %in% c('pcterrmodSE', 'MCSE_pcterrmodSE')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(theta, measure, Method))
    pcterrmodSErows <- cbind(params, pcterrmodSEvals)
    allpcterrmodSErowsHH <- rbind(allpcterrmodSErowsHH, pcterrmodSErows)
    
    # Average interval length
    intlenvals <- measures %>%
      filter(measure %in% c('avgintlength', 'MCSE_avgintlength')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(theta, measure, Method))
    intlenrows <- cbind(params, intlenvals)
    allintlenrowsHH <- rbind(allintlenrowsHH, intlenrows)
  } else {
    pars <- c('theta', 'WPICC', 'CAC', 'BPICC')
    
    # Bias values
    biasvals <- measures %>%
      filter(measure %in% c('bias', 'MCSE_bias')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(all_of(pars), measure, Method))
    biasrows <- cbind(params, biasvals)
    allbiasrowsHG <- rbind(allbiasrowsHG, biasrows)
    
    # MSE values
    MSEvals <- measures %>%
      filter(measure %in% c('MSE', 'MCSE_MSE')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(all_of(pars), measure, Method))
    MSErows <- cbind(params, MSEvals)
    allMSErowsHG <- rbind(allMSErowsHG, MSErows)
    
    # Coverage values
    covvals <- measures %>%
      filter(measure %in% c('coverage', 'MCSE_coverage')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(theta, measure, Method))
    covrows <- cbind(params, covvals)
    allcovrowsHG <- rbind(allcovrowsHG, covrows)
    
    # Empirical SE
    empSEvals <- measures %>%
      filter(measure %in% c('empSE', 'MCSE_empSE')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(all_of(pars), measure, Method))
    empSErows <- cbind(params, empSEvals)
    allempSErowsHG <- rbind(allempSErowsHG, empSErows)
    
    # Average model-based SE
    modSEvals <- measures %>%
      filter(measure %in% c('avgmodSE', 'MCSE_avgmodSE')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(theta, measure, Method))
    modSErows <- cbind(params, modSEvals)
    allmodSErowsHG <- rbind(allmodSErowsHG, modSErows)
    
    # Relative % error in model-based SE
    pcterrmodSEvals <- measures %>%
      filter(measure %in% c('pcterrmodSE', 'MCSE_pcterrmodSE')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(theta, measure, Method))
    pcterrmodSErows <- cbind(params, pcterrmodSEvals)
    allpcterrmodSErowsHG <- rbind(allpcterrmodSErowsHG, pcterrmodSErows)
    
    # Average interval length
    intlenvals <- measures %>%
      filter(measure %in% c('avgintlength', 'MCSE_avgintlength')) %>%
      mutate(Method=recode_factor(method, MCMC="Bayesian", REML="REML")) %>%
      select(c(theta, measure, Method))
    intlenrows <- cbind(params, intlenvals)
    allintlenrowsHG <- rbind(allintlenrowsHG, intlenrows)
  }
}

# Convert results to long format for formatting into a table
pars <- c('theta', 'WPICC')
biasdfHH <- allbiasrowsHH %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = bias, MCSE = MCSE_bias))

pars <- c('theta', 'WPICC', 'CAC', 'BPICC')
biasdfHG <- allbiasrowsHG %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = bias, MCSE = MCSE_bias))

pars <- c('theta', 'WPICC')
MSEdfHH <- allMSErowsHH %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = MSE, MCSE = MCSE_MSE))

pars <- c('theta', 'WPICC', 'CAC', 'BPICC')
MSEdfHG <- allMSErowsHG %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = MSE, MCSE = MCSE_MSE))

pars <- c('theta')
covdfHH <- allcovrowsHH %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = coverage, MCSE = MCSE_coverage))

pars <- c('theta')
covdfHG <- allcovrowsHG %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = coverage, MCSE = MCSE_coverage))

pars <- c('theta', 'WPICC')
empSEdfHH <- allempSErowsHH %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = empSE, MCSE = MCSE_empSE))

pars <- c('theta', 'WPICC', 'CAC', 'BPICC')
empSEdfHG <- allempSErowsHG %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = empSE, MCSE = MCSE_empSE))

pars <- c('theta')
modSEdfHH <- allmodSErowsHH %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = avgmodSE, MCSE = MCSE_avgmodSE))

pars <- c('theta')
modSEdfHG <- allmodSErowsHG %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = avgmodSE, MCSE = MCSE_avgmodSE))

pars <- c('theta')
pcterrmodSEdfHH <- allpcterrmodSErowsHH %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = pcterrmodSE, MCSE = MCSE_pcterrmodSE))

pars <- c('theta')
pcterrmodSEdfHG <- allpcterrmodSErowsHG %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = pcterrmodSE, MCSE = MCSE_pcterrmodSE))

pars <- c('theta')
intlendfHH <- allintlenrowsHH %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = avgintlength, MCSE = MCSE_avgintlength))

pars <- c('theta')
intlendfHG <- allintlenrowsHG %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  ) %>%
  pivot_wider(
    names_from='measure',
    values_from='value'
  ) %>%
  rename(c(value = avgintlength, MCSE = MCSE_avgintlength))
