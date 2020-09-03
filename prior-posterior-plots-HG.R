# Plot prior and marginal posterior densities against true parameter values
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

rm(list=ls())
load('mcmc_results.Rda')
library(tidyverse)
library(invgamma)

params <- c(
  "theta", "WPICC", "CAC", "sig_sq_subject", "sig_sq_cluster", "sig_sq_cp",
  paste0('C',1:12)
)

plot_dens_true <- function(df, dfvals, dftrue, dftruevals, fillcolor){
  # Plots a density (prior or posterior) against the true parameter value
  # for all parameters
  
  orderedparams <- params
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

plot_dens_true_all <- function(df, x, prior, post, dftrue, dftruevals){
  # Plots prior and posterior densities against the true parameter value
  # for all parameters, with x-axis range adjusted to show 99.8% of the
  # posterior probability mass
  
  orderedparams <- params
  df %>%
    mutate(
      parameters = ordered(parameters, levels=orderedparams)
    ) %>%
    ggplot(aes(x=x, y=post)) +
    geom_area(fill="darkseagreen", alpha=0.5) +
    geom_area(aes(x=x, y=prior), fill="darkorange", alpha=0.5) +
    geom_vline(aes(xintercept=dftruevals), color="steelblue", size=1, data=dftrue) +
    facet_wrap(~parameters, scales="free") +
    theme_bw()
}

sig2c_implied <- function(sig2e, WPICC, CAC){
  # Calculates implied cluster variance using draws from
  # subject variance and within-cluster correlation components
  
  return(CAC * sig2e * WPICC / (1 - WPICC))
}

sig2cp_implied <- function(sig2e, WPICC, sig2c){
  # Calculates implied between-cluster-period variance using
  # draws from subject and cluster variances and within-period ICC

  return(sig2e * WPICC / (1 - WPICC) - sig2c)
}

# Prepares posterior draws for plotting
Cnames <- paste('C[', 1:clusters, ']', sep='')
pars <- c('theta', 'WPICC', 'CAC', 'sig_sq_subject', 'sig_sq_cluster',
          'sig_sq_cp', Cnames)
draws <- as.data.frame(fit, pars=pars)
parnames <- gsub("[][]", "", pars)
colnames(draws) <- parnames
longdraws <- draws %>% 
  pivot_longer(
    everything(), names_to='parameters', values_to='samplevals'
  )

# Prepares true parameter values for plotting
true_vals <- c(theta, WPICC, CAC, sig_sq_subject, sig_sq_cluster, sig_sq_cp, C)
actual <- data.frame(parameters=parnames, truevals=true_vals)

# Redefines priors (based on associated Stan file) and prepares for plotting
theta_prior <- runif(2e4, -1e4, 1e4)
WPICC_prior <- rbeta(2e4, 1.2, 5)
CAC_prior <- rbeta(2e4, 5, 2)
sig2e_prior <- rinvgamma(2e4, 2, 2)
sig2c_prior <- sig2c_implied(sig2e_prior, WPICC_prior, CAC_prior)
sig2cp_prior <- sig2cp_implied(sig2e_prior, WPICC_prior, sig2c_prior)
C_prior <- rnorm(2e4, 0, sqrt(sig_sq_cluster))
Cmat <- matrix(rep(C_prior, times=clusters), nrow=length(C_prior), ncol=clusters)
Cdf <- data.frame(Cmat)
colnames(Cdf) <- gsub("[][]", "", Cnames)
priors <- data.frame(theta=theta_prior, WPICC=WPICC_prior, CAC=CAC_prior,
                     sig_sq_subject=sig2e_prior, sig_sq_cluster=sig2c_prior,
                     sig_sq_cp=sig2cp_prior, Cdf)
longpriors <- priors %>%
  pivot_longer(
    everything(), names_to='parameters', values_to='priorvals'
  )

# Benchmark plot of posteriors against true values
mcmc_recover_hist(draws, true_vals)

# Plot prior distributions with true parameter values
plot_dens_true(longpriors, longpriors$priorvals, actual, actual$truevals, "darkorange") +
  labs(
    title='Priors and true values',
    x='Value'
  )

# Plot marginal posterior distributions with true parameter values
plot_dens_true(longdraws, longdraws$samplevals, actual, actual$truevals, "darkseagreen") +
  labs(
    title='Posteriors and true values',
    x='Value'
  )

# Calculates densities for the range of values that covers 99.8% of the
# posterior probability mass (for diffuse priors, this will trim the priors
# so we can view the densities simultaneously without them being on very
# different scales)
vals <- NULL
for (par in colnames(draws)){
  edges <- quantile(draws[[par]], c(0.001, 0.999))
  post_dens <- density(draws[[par]], n=1000, from=edges[1], to=edges[2])
  prior_dens <- density(priors[[par]], n=1000, from=edges[1], to=edges[2])
  parvals <- tibble(
    parameters=par,
    x=post_dens$x,
    prior=prior_dens$y,
    post=post_dens$y
  )
  vals <- bind_rows(vals, parvals)
}

# Plot prior and posterior distributions with true parameter values
plot_dens_true_all(vals, vals$x, vals$prior, vals$post, actual, actual$truevals) +
  labs(
    title='Priors, posteriors and true values',
    x='Value'
  )
