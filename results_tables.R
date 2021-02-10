# Combine performance_measures files across trial configurations
# into a single data frame for formatting into a results table
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

library(tidyverse)
library(tables)
library(grid)
library(gridExtra)

file.names <- dir('./performance_measures')

allbiasrowsHH <- data.frame()
allMSErowsHH <- data.frame()
allcovrowsHH <- data.frame()
allempSErowsHH <- data.frame()
allmodSErowsHH <- data.frame()
allbiasrowsHG <- data.frame()
allMSErowsHG <- data.frame()
allcovrowsHG <- data.frame()
allempSErowsHG <- data.frame()
allmodSErowsHG <- data.frame()
for (i in 1:length(file.names)) {
  
  load(paste0('performance_measures/', file.names[i]))
  
  params <- results$params
  measures <- results$measures
  
  if (params$r==1.0) {
    pars <- c('theta', 'WPICC')
    parnames <- c("$\\theta$", "$\\rho_1$")
    
    # Bias values
    biasvals <- measures %>%
      filter(measure %in% c('bias', 'MCSE_bias')) %>%
      select(c(all_of(pars), measure, method))
    biasrows <- cbind(params, biasvals)
    allbiasrowsHH <- rbind(allbiasrowsHH, biasrows)
    
    # MSE values
    MSEvals <- measures %>%
      filter(measure %in% c('MSE', 'MCSE_MSE')) %>%
      select(c(all_of(pars), measure, method))
    MSErows <- cbind(params, MSEvals)
    allMSErowsHH <- rbind(allMSErowsHH, MSErows)
    
    # Coverage values
    covvals <- measures %>%
      filter(measure %in% c('coverage', 'MCSE_coverage')) %>%
      select(c(theta, measure, method))
    covrows <- cbind(params, covvals)
    allcovrowsHH <- rbind(allcovrowsHH, covrows)
    
    # Empirical SE
    empSEvals <- measures %>%
      filter(measure %in% c('empSE', 'MCSE_empSE')) %>%
      select(c(all_of(pars), measure, method))
    empSErows <- cbind(params, empSEvals)
    allempSErowsHH <- rbind(allempSErowsHH, empSErows)
    
    # Average model-based SE
    modSEvals <- measures %>%
      filter(measure %in% c('avgmodSE', 'MCSE_avgmodSE')) %>%
      select(c(all_of(pars), measure, method))
    modSErows <- cbind(params, modSEvals)
    allmodSErowsHH <- rbind(allmodSErowsHH, modSErows)
  } else {
    pars <- c('theta', 'WPICC', 'BPICC')
    parnames <- c("$\\theta$", "$\\rho_1$", "$\\rho_2$")
    
    # Bias values
    biasvals <- measures %>%
      filter(measure %in% c('bias', 'MCSE_bias')) %>%
      select(c(all_of(pars), measure, method))
    biasrows <- cbind(params, biasvals)
    allbiasrowsHG <- rbind(allbiasrowsHG, biasrows)
    
    # MSE values
    MSEvals <- measures %>%
      filter(measure %in% c('MSE', 'MCSE_MSE')) %>%
      select(c(all_of(pars), measure, method))
    MSErows <- cbind(params, MSEvals)
    allMSErowsHG <- rbind(allMSErowsHG, MSErows)
    
    # Coverage values
    covvals <- measures %>%
      filter(measure %in% c('coverage', 'MCSE_coverage')) %>%
      select(c(theta, measure, method))
    covrows <- cbind(params, covvals)
    allcovrowsHG <- rbind(allcovrowsHG, covrows)
    
    # Empirical SE
    empSEvals <- measures %>%
      filter(measure %in% c('empSE', 'MCSE_empSE')) %>%
      select(c(all_of(pars), measure, method))
    empSErows <- cbind(params, empSEvals)
    allempSErowsHG <- rbind(allempSErowsHG, empSErows)
    
    # Average model-based SE
    modSEvals <- measures %>%
      filter(measure %in% c('avgmodSE', 'MCSE_avgmodSE')) %>%
      select(c(all_of(pars), measure, method))
    modSErows <- cbind(params, modSEvals)
    allmodSErowsHG <- rbind(allmodSErowsHG, modSErows)
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

pars <- c('theta', 'WPICC', 'BPICC')
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

pars <- c('theta', 'WPICC', 'BPICC')
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

# TODO: duplicate remaining measures for HH and HG
covdf <- allcovrows %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  )

empSEdf <- allempSErows %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  )

modSEdf <- allmodSErows %>%
  pivot_longer(
    cols=all_of(pars),
    names_to='parameters',
    values_to='value'
  )

# Create formatted results table

# Create dummy results dataframe (for now)
dat <- expand.grid(
  S = c(1,2,5),
  Tp = c(5, 9),
  m = c(10, 100),
  rho1 = c(0.05, 0.1),
  r = 0.8,
  method = c('MCMC', 'REML'),
  measure = c('bias', 'MCSE_bias')
)
dat$theta <- rnorm(dim(dat)[1])
dat$WPICC <- rnorm(dim(dat)[1])
dat$BPICC <- rnorm(dim(dat)[1])

biasdf <- dat %>%
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

# Create Latex code for generating table
fmtmain <- function(x, digits){
  s <- format(x, digits=digits)
  s <- latexNumeric(s)
  s
}

fmtMCSE <- function(x, digits){
  s <- format(x, digits=digits)
  sf <- sprintf("$(%s)$", s)
  sf
}

create_Latex_table <- function(df, maindigits=1, MCSEdigits=1, pars, parnames) {
  toLatex(
    tabular(
      RowFactor(rho1, name="$\\rho_1$", spacing=1)
      * RowFactor(Tp, name="T", spacing=1)
      * RowFactor(m, spacing=1)
      * Factor(S)
      ~ Heading()*RowFactor(factor(parameters, levels=pars), levelnames=parnames)
      * Heading()*RowFactor(method)
      * Heading()*(Format(fmtmain(digits=maindigits))*value +
                   Format(fmtMCSE(digits=MCSEdigits))*MCSE)
      * Heading()*identity,
      data=df
    )
  )
}

create_Latex_table(biasdfHH, pars=c('theta', 'WPICC'), parnames=c("$\\theta$", "$\\rho_1$"))
create_Latex_table(biasdfHG, pars=c('theta', 'WPICC', 'BPICC'), parnames=c("$\\theta$", "$\\rho_1$", "$\\rho_2$"))

create_Latex_table(MSEdf, maindigits=2, MCSEdigits=2)
