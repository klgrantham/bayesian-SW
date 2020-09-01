# Plot prior and marginal posterior densities against true parameter values
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

rm(list=ls())
load('mcmc_results.Rda')
library(tidyverse)

plot_dens_true <- function(df, dfvals, dftrue, dftruevals, fillcolor){
  orderedparams <- c(
    "rho", "theta", "sig_sq_subject", "sig_sq_cluster",
    paste0('C',1:12)
  )
  df %>%
    mutate(
      parameters = ordered(parameters, levels=orderedparams)
    ) %>%
    ggplot(aes(x=dfvals)) +
      geom_density(aes(x=dfvals), fill=fillcolor, alpha=0.5) +
      geom_vline(aes(xintercept=dftruevals), color="steelblue", size=1, data=dftrue) +
      facet_wrap(~parameters, scales="free") +
      theme_bw()
}

sig2c_implied <- function(sig2e, rho){
  # Calculates implied cluster variance using draws from
  # subject variance and within-cluster correlation
  return(sig2e * rho/(1 - rho))
}

sig_sq_cluster <- sig2c_implied(sig_sq_subject, rho)
Cnames <- paste('C[', 1:clusters, ']', sep='')
pars <- c('rho', 'theta', 'sig_sq_subject', 'sig_sq_cluster', Cnames)
draws <- as.data.frame(fit, pars=pars)
parnames <- gsub("[][]", "", pars)
colnames(draws) <- parnames
gatherdraws <- draws %>% gather("parameters", "samplevals")

true_vals <- c(rho, theta, sig_sq_subject, sig_sq_cluster, C)
actual <- data.frame(parameters=parnames, truevals=true_vals)

rho_prior <- runif(2e4, 0, 1)
theta_prior <- runif(2e4, -1e3, 1e3)
sig2e_prior <- exp(runif(2e4, -5, 5))
sig2c_prior <- sig2c_implied(sig2e_prior, rho_prior)
C_prior <- rnorm(2e4, 0, sqrt(sig_sq_cluster))
Cmat <- matrix(rep(C_prior, times=clusters), nrow=length(C_prior), ncol=clusters)
Cdf <- data.frame(Cmat)
colnames(Cdf) <- gsub("[][]", "", Cnames)
priors <- data.frame(rho=rho_prior, theta=theta_prior,
                     sig_sq_subject=sig2e_prior, sig_sq_cluster=sig2c_prior,
                     Cdf)
gatherpriors <- priors %>% gather("parameters", "priorvals")

mcmc_recover_hist(draws, true_vals)

# Plot prior distributions with true parameter values
plot_dens_true(gatherpriors, gatherpriors$priorvals, actual, actual$truevals, "darkorange") +
  labs(
    title='Priors and true values',
    x='Value'
  )

# Plot marginal posterior distributions with true parameter values
plot_dens_true(gatherdraws, gatherdraws$samplevals, actual, actual$truevals, "darkseagreen") +
  labs(
    title='Posteriors and true values',
    x='Value'
  )

