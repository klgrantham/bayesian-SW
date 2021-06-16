# Generate results tables
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

library(tidyverse)
library(tables)
library(grid)
library(gridExtra)

source('process-results.R')

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

intlendftheta <- rbind(intlendfHH, intlendfHG)
toLatex(
  tabular(
    RowFactor(rho1, name="$\\rho_1$", spacing=1)
    * RowFactor(Tp, name="T", spacing=1)
    * RowFactor(m, spacing=1)
    * Factor(S)
    ~ Factor(r)
    * Heading()*RowFactor(Method)
    * Heading()*(Format(fmtmain(digits=1))*value)
    * Heading()*identity,
    data=intlendftheta
  )
)

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