# Combine performance_measures files across trial configurations
# into a single data frame for formatting into a results table
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

library(tidyverse)
library(tables)
library(knitr)
library(kableExtra)
library(gridExtra)

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
  # May need to format values with desired number of digits here
  # as Format() in tables package doesn't seem to work
  
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

# Create formatted results table

# Create dummy results dataframe (for now)
dat <- expand.grid(
  S = c(1,2,5),
  Tp = c(5, 9),
  m = c(10, 100),
  rho1 = c(0.05, 0.1),
  r = 0.8,
  method = c('MCMC', 'REML')
)
dat$theta <- rnorm(dim(dat)[1])
dat$WPICC <- rnorm(dim(dat)[1])
dat$BPICC <- rnorm(dim(dat)[1])

biasdf <- dat %>%
  pivot_longer(
    cols=c(theta, WPICC, BPICC),
    names_to='parameters',
    values_to='bias'
  )

# Create Latex code for generating table
toLatex(
  tabular(
    RowFactor(rho1, spacing=1)
    * RowFactor(Tp, spacing=1)
    * RowFactor(m, spacing=1)
    * Factor(S)
    ~ Heading()*RowFactor(parameters)
    * Heading()*RowFactor(method)
    * Heading()*bias
    * Heading()*identity,
    data=biasdf
  )
)

# kable and kableExtra don't seem to be able to generate
# Latex code and can't seem to handle nested columns and rows as well
kbl(biasdf, booktabs=T) %>%
  kable_styling() %>%
  collapse_rows(columns=1:4, valign="top")


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
    ggplot(aes(x=S, y=bias, group=method)) +
      geom_line(aes(color=method), show.legend=F) +
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
  ggplot(aes(x=S, y=bias, group=method)) +
  geom_line(aes(color=method), show.legend=F) +
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
  ggplot(aes(x=S, y=bias, group=method)) +
  geom_line(aes(color=method), show.legend=F) +
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
  ggplot(aes(x=S, y=bias, group=method)) +
  geom_line(aes(color=method), show.legend=F) +
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