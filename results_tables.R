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

# Table 3: Calculate percentage of valid replicates
allreprows <- allreprows %>%
  mutate(Bayesian=MCMCreps/1e3*100, REML=REMLreps/1e3*100)
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

# Table 4: Bias for theta
biasdftheta <- rbind(biasdfHH %>% filter(parameters=='theta'),
                     biasdfHG %>% filter(parameters=='theta'))
biasdftheta_scaled <- biasdftheta %>% 
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(biasdftheta_scaled, 1, 1)

# Table 5: MSE for theta
MSEdftheta <- rbind(MSEdfHH %>% filter(parameters=='theta'),
                    MSEdfHG %>% filter(parameters=='theta'))
MSEdftheta_scaled <- MSEdftheta %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(MSEdftheta_scaled, 1, 1)

# Table 6: Confidence/credible interval coverage for theta
covdftheta <- rbind(covdfHH, covdfHG)
create_Latex_table(covdftheta, 3, 1)

# Table 7: Relative % error in model-based SE for theta
pcterrmodSEdftheta <- rbind(pcterrmodSEdfHH, pcterrmodSEdfHG)
create_Latex_table(pcterrmodSEdftheta, 1, 2)

# Table 8: Confidence/credible interval width for theta
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

# Table 9: Bias for WPICC
biasdfWPICC <- rbind(biasdfHH %>% filter(parameters=='WPICC'),
                     biasdfHG %>% filter(parameters=='WPICC'))
biasdfWPICC_scaled <- biasdfWPICC %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(biasdfWPICC_scaled, 1, 1)

# Table 10: MSE for WPICC
MSEdfWPICC <- rbind(MSEdfHH %>% filter(parameters=='WPICC'),
                    MSEdfHG %>% filter(parameters=='WPICC'))
MSEdfWPICC_scaled <- MSEdfWPICC %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(MSEdfWPICC_scaled, 1, 1)

# Table 11: Bias for CAC
biasdfCAC <- biasdfHG %>% filter(parameters=='CAC')
biasdfCAC_scaled <- biasdfCAC %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(biasdfCAC_scaled, 1, 1)

# Table 12: MSE for CAC
MSEdfCAC <- MSEdfHG %>% filter(parameters=='CAC')
MSEdfCAC_scaled <- MSEdfCAC %>%
  rename(value_unscaled=value, MCSE_unscaled=MCSE) %>%
  mutate(value=value_unscaled*1e4, MCSE=MCSE_unscaled*1e4)
create_Latex_table(MSEdfCAC_scaled, 1, 1)
