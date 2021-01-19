# Combine performance_measures files across trial configurations
# into a single data frame for formatting into a results table
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

library(tidyverse)
library(tables)
library(gridExtra)

file.names <- dir('./performance_measures')

allbiasrows <- data.frame()
allMSErows <- data.frame()
allcovrows <- data.frame()
allempSErows <- data.frame()
allmodSErows <- data.frame()
for (i in 1:length(file.names)) {
  
  load(paste0('performance_measures/', file.names[i]))
  
  params <- results$params
  measures <- results$measures
  
  if (params$r==1.0) {
    pars <- c('theta', 'WPICC')
    parnames <- c("$\\theta$", "$\\rho_1$")    
  } else {
    pars <- c('theta', 'WPICC', 'BPICC')
    parnames <- c("$\\theta$", "$\\rho_1$", "$\\rho_2$")
  }
  
  # Bias values
  biasvals <- measures %>%
    filter(measure %in% c('bias', 'MCSE_bias')) %>%
    select(c(all_of(pars), measure, method))
  biasrows <- cbind(params, biasvals)
  allbiasrows <- rbind(allbiasrows, biasrows)
  # May need to format values with desired number of digits here
  # as Format() in tables package doesn't seem to work
  
  # MSE values
  MSEvals <- measures %>%
    filter(measure %in% c('MSE', 'MCSE_MSE')) %>%
    select(c(all_of(pars), measure, method))
  MSErows <- cbind(params, MSEvals)
  allMSErows <- rbind(allMSErows, MSErows)
  
  # Coverage values
  covvals <- measures %>%
    filter(measure %in% c('coverage', 'MCSE_coverage')) %>%
    select(c(theta, measure, method))
  covrows <- cbind(params, covvals)
  allcovrows <- rbind(allcovrows, covrows)
  
  # Empirical SE
  empSEvals <- measures %>%
    filter(measure %in% c('empSE', 'MCSE_empSE')) %>%
    select(c(all_of(pars), measure, method))
  empSErows <- cbind(params, empSEvals)
  allempSErows <- rbind(allempSErows, empSErows)
  
  # Average model-based SE
  modSEvals <- measures %>%
    filter(measure %in% c('avgmodSE', 'MCSE_avgmodSE')) %>%
    select(c(all_of(pars), measure, method))
  modSErows <- cbind(params, modSEvals)
  allmodSErows <- rbind(allmodSErows, modSErows)
}

# Convert results to long format for formatting into a table
biasdf <- allbiasrows %>%
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

MSEdf <- allMSErows %>%
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

create_Latex_table <- function(df, maindigits=1, MCSEdigits=1) {
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

create_Latex_table(biasdf)

create_Latex_table(MSEdf, maindigits=2, MCSEdigits=2)


# Plot contents of results table
Tp.labs <- c("T = 5", "T = 9")
names(Tp.labs) <- c("5", "9")

m.labs <- c("m = 10", "m = 100")
names(m.labs) <- c("10", "100")

plot_measure <- function(measure) {
  ggplot(aes(x=S, y=measure, group=method)) +
    geom_line(aes(color=method)) +
    facet_grid(
      m ~ Tp,
      labeller = labeller(Tp = Tp.labs, m = m.labs)
    ) +
    theme(legend.position="none") +
    theme_bw()
}

p1 <- biasdf %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta') %>%
    ggplot(aes(x=S, y=value, group=method)) +
      geom_line(aes(color=method), show.legend=FALSE) +
      geom_point(aes(color=method), show.legend=FALSE) +
      facet_grid(
        m ~ Tp,
        labeller = labeller(Tp = Tp.labs, m = m.labs)
      ) +
      theme_bw()  +
      theme(
        strip.background = element_rect(
          color="white", fill="white", linetype="solid"
        )
      )
p2 <- biasdf %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta') %>%
  ggplot(aes(x=S, y=value, group=method)) +
  geom_point(aes(color=method), show.legend=FALSE) +
  geom_line(aes(color=method), show.legend=FALSE) +
  facet_grid(
    m ~ Tp,
    labeller = labeller(Tp = Tp.labs, m = m.labs)
  ) +
  theme_bw()  +
  theme(
    strip.background = element_rect(
      color="white", fill="white", linetype="solid"
    )
  )
p3 <- biasdf %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta') %>%
  ggplot(aes(x=S, y=value, group=method)) +
  geom_line(aes(color=method), show.legend=FALSE) +
  geom_point(aes(color=method), show.legend=FALSE) +
  facet_grid(
    m ~ Tp,
    labeller = labeller(Tp = Tp.labs, m = m.labs)
  ) +
  theme_bw()  +
  theme(
    strip.background = element_rect(
      color="white", fill="white", linetype="solid"
    )
  )
p4 <- biasdf %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta') %>%
  ggplot(aes(x=S, y=value, group=method)) +
  geom_line(aes(color=method), show.legend=FALSE) +
  geom_point(aes(color=method), show.legend=FALSE) +
  facet_grid(
    m ~ Tp,
    labeller = labeller(Tp = Tp.labs, m = m.labs)
  ) +
  theme_bw()  +
  theme(
    strip.background = element_rect(
      color="white", fill="white", linetype="solid"
    )
  )
#grid.arrange(p1, p2, p3, p4, nrow=2, ncol=2)
grid.arrange(p1, p2, nrow=2)
# TODO: include legend at bottom of arranged plot