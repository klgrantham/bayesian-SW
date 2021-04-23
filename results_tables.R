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


# Create formatted results table
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

create_Latex_table_multpars <- function(df, maindigits=1, MCSEdigits=1, pars, parnames) {
  toLatex(
    tabular(
      RowFactor(rho1, name="$\\rho_1$", spacing=1)
      * RowFactor(Tp, name="T", spacing=1)
      * RowFactor(m, spacing=1)
      * Factor(S)
      ~ Heading()*RowFactor(factor(parameters, levels=pars), levelnames=parnames)
      * Heading()*RowFactor(Method)
      * Heading()*(Format(fmtmain(digits=maindigits))*value +
                   Format(fmtMCSE(digits=MCSEdigits))*MCSE)
      * Heading()*identity,
      data=df
    )
  )
}

create_Latex_table <- function(df, maindigits=1, MCSEdigits=1) {
  toLatex(
    tabular(
      RowFactor(rho1, name="$\\rho_1$", spacing=1)
      * RowFactor(Tp, name="T", spacing=1)
      * RowFactor(m, spacing=1)
      * Factor(S)
      ~ Factor(r)
      * Heading()*RowFactor(Method)
      * Heading()*(Format(fmtmain(digits=maindigits))*value +
                   Format(fmtMCSE(digits=MCSEdigits))*MCSE)
      * Heading()*identity,
      data=df
    )
  )
}

biasdftheta <- rbind(biasdfHH %>% filter(parameters=='theta'),
                     biasdfHG %>% filter(parameters=='theta'))
create_Latex_table(biasdftheta, 2, 1)

biasdftheta_scaled <- biasdftheta %>% 
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(biasdftheta_scaled, 1, 1)

MSEdftheta <- rbind(MSEdfHH %>% filter(parameters=='theta'),
                    MSEdfHG %>% filter(parameters=='theta'))
create_Latex_table(MSEdftheta, 2, 1)

MSEdftheta_scaled <- MSEdftheta %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(MSEdftheta_scaled, 1, 1)

covdftheta <- rbind(covdfHH, covdfHG)
create_Latex_table(covdftheta, 2, 1)

pcterrmodSEdftheta <- rbind(pcterrmodSEdfHH, pcterrmodSEdfHG)
create_Latex_table(pcterrmodSEdftheta, 1, 2)

biasdfWPICC <- rbind(biasdfHH %>% filter(parameters=='WPICC'),
                     biasdfHG %>% filter(parameters=='WPICC'))
create_Latex_table(biasdfWPICC, 2, 1)

biasdfWPICC_scaled <- biasdfWPICC %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(biasdfWPICC_scaled, 1, 1)

MSEdfWPICC <- rbind(MSEdfHH %>% filter(parameters=='WPICC'),
                    MSEdfHG %>% filter(parameters=='WPICC'))
create_Latex_table(MSEdfWPICC, 2, 1)

MSEdfWPICC_scaled <- MSEdfWPICC %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(MSEdfWPICC_scaled, 1, 1)

biasdfCAC <- biasdfHG %>% filter(parameters=='CAC')
create_Latex_table(biasdfCAC, 2, 1)

biasdfCAC_scaled <- biasdfCAC %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(biasdfCAC_scaled, 1, 1)

MSEdfCAC <- MSEdfHG %>% filter(parameters=='CAC')
create_Latex_table(MSEdfCAC, 2, 1)

MSEdfCAC_scaled <- MSEdfCAC %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(MSEdfCAC_scaled, 1, 1)


create_Latex_table_multpars(biasdfHH, maindigits=2, MCSEdigits=1,
                   pars=c('theta', 'WPICC'), parnames=c("$\\theta$", "$\\rho_1$"))
create_Latex_table_multpars(biasdfHG, maindigits=2, MCSEdigits=1,
                   pars=c('theta', 'WPICC', 'CAC', 'BPICC'),
                   parnames=c("$\\theta$", "$\\rho_1$", "r", "$\\rho_2$"))
create_Latex_table_multpars(MSEdfHH, maindigits=2, MCSEdigits=1,
                   pars=c('theta', 'WPICC'), parnames=c("$\\theta$", "$\\rho_1$"))

# Calculate success rates
allreprows <- allreprows %>%
  mutate(Bayesian=MCMCreps/1e3*100, REML=REMLreps/1e3*100)

toLatex(
  tabular(
    Factor(r, name="$r$")
    * Factor(rho1, name="$\\rho_1$")
    * Factor(Tp, name="T")
    * Factor(m)
    * Factor(S)
    ~ (Format(fmtmain(digits=3))*Bayesian +
      Format(fmtmain(digits=3))*REML)
    * Heading()*identity,
    data=allreprows
  )
)

# Widen by putting r values in different columns
toLatex(
  tabular(
    Factor(rho1, name="$\\rho_1$")
    * Factor(Tp, name="T")
    * Factor(m)
    * Factor(S)
    ~ Factor(r, name="$r$")
    * (Format(fmtmain(digits=3))*Bayesian +
       Format(fmtmain(digits=3))*REML)
    * Heading()*identity,
    data=allreprows
  )
)