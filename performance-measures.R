# Calculate performance measures for models with MCMC and REML estimation
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(rstan)
library(lme4)
library(tidyverse)
library(parameters)

collate_results <- function(Nrep, clust_per_seq, periods, subjects, WPICC, CAC, theta) {
  # Collate results across (valid) replicates and compute performance measures
  
  MCMC_means <- data.frame()
  MCMC_medians <- data.frame()
  MCMC_sds <- data.frame()
  MCMC_025 <- data.frame()
  MCMC_975 <- data.frame()
  MCMC_powvals <- data.frame()
  REML_ests <- data.frame()
  REML_stderrs <- data.frame()
  REML_stderrKRs <- data.frame()
  REML_CI_KR_lowers <- data.frame()
  REML_CI_KR_uppers <- data.frame()
  
  WPICC100 <- WPICC*100
  CAC100 <- CAC*100
  config <- paste0(
    '_S', clust_per_seq,
    '_T', periods,
    '_m', subjects,
    '_WPICC', WPICC100,
    '_CAC', CAC100,
    '_theta', theta
  )
  
  nonzerovar <- 0
  for (n in 1:Nrep) {
    
    infile=paste0(
      'reduced_results/',
      'reduced_results',
      config,
      sprintf('_seed%04d', n),
      '.Rda'
    )
    load(file=infile)
    # Contains: res$MCMC_res, res$REML_res, res$truevals, res$div, res$zerovar

    REML <- res$REML_res
    MCMC <- res$MCMC_res
    
    # Count REML results with all nonzero variance components
    if (isFALSE(res$zerovar)) {
      nonzerovar <- nonzerovar + 1
    }
    # Include nth REML results in results containers
    REML_ests <- rbind(REML_ests, REML$est)
    REML_stderrs <- rbind(REML_stderrs, REML$stderr)
    REML_stderrKRs <- rbind(REML_stderrKRs, REML$stderrKR)
    REML_CI_KR_lowers <- rbind(REML_CI_KR_lowers, REML$CI_KR_lower)
    REML_CI_KR_uppers <- rbind(REML_CI_KR_uppers, REML$CI_KR_upper)
    
    # Skip nth MCMC results if any divergent transitions
    if (res$div > 0) {
      next
    }
    # Include nth MCMC results in results containers 
    MCMC_means <- rbind(MCMC_means, MCMC$post_mean)
    MCMC_medians <- rbind(MCMC_medians, MCMC$post_median)
    MCMC_sds <- rbind(MCMC_sds, MCMC$post_sd)
    MCMC_025 <- rbind(MCMC_025, MCMC$post_025)
    MCMC_975 <- rbind(MCMC_975, MCMC$post_975)
    MCMC_powvals <- rbind(MCMC_powvals, MCMC$theta_greaterthanC)
  }
  rnames_MCMC <- rownames(MCMC)
  rnames_REML <- rownames(REML)
  colnames(REML_ests) <- rnames_REML
  colnames(REML_stderrs) <- rnames_REML
  colnames(REML_stderrKRs) <- rnames_REML
  colnames(REML_CI_KR_lowers) <- rnames_REML
  colnames(REML_CI_KR_uppers) <- rnames_REML
  colnames(MCMC_means) <- rnames_MCMC
  colnames(MCMC_medians) <- rnames_MCMC
  colnames(MCMC_sds) <- rnames_MCMC
  colnames(MCMC_025) <- rnames_MCMC
  colnames(MCMC_975) <- rnames_MCMC
  colnames(MCMC_powvals) <- rnames_MCMC
  
  true_params <- res$truevals
  
  # Calculate performance measures across all valid replicates
  
  # Number of valid MCMC replicates (those with divergent transitions excluded)
  MCMCreps <- dim(MCMC_means)[1]
  
  # Number of "valid" REML replicates
  REMLreps <- nonzerovar
  
  # Calculate bias
  MCMC_mean_bias <- bias(MCMC_means, true_params, 'MCMC_mean_bias')
  MCMC_median_bias <- bias(MCMC_medians, true_params, 'MCMC_median_bias')
  REML_bias <- bias(REML_ests, true_params, 'REML_bias')
  
  # Calculate MSE
  MCMC_mean_mse <- mse(MCMC_means, true_params, 'MCMC_mean_mse')
  MCMC_median_mse <- mse(MCMC_medians, true_params, 'MCMC_median_mse')
  REML_mse <- mse(REML_ests, true_params, 'REML_mse')
  
  # Calculate Monte Carlo standard deviation
  MCMC_mean_MCsd <- MC.sd(MCMC_means, 'MCMC_mean_MCsd')
  REML_MCsd <- MC.sd(REML_ests, 'REML_MCsd')
  
  # Calculate average standard errors
  MCMC_avg_SD <- as.data.frame(t(colMeans(MCMC_sds)))
  MCMC_avg_SD$measure <- 'MCMC_avg_SD'
  REML_avg_KR_SE <- as.data.frame(t(colMeans(REML_stderrKRs)))
  REML_avg_KR_SE$measure <- 'REML_avg_KR_SE'
  REML_avg_SE <- as.data.frame(t(colMeans(REML_stderrs)))
  REML_avg_SE$measure <- 'REML_avg_SE'
  
  # Calculate interval coverage
  true_params_MCMCrep <- true_params[rep(1, MCMCreps),] # true_params is a row vector
  MCMC_coverage <- coverage(MCMC_025, MCMC_975, true_params_MCMCrep, 'MCMC_coverage')
  true_params_REMLrep <- true_params[rep(1, REMLreps),]
  REML_KR_coverage <- coverage(REML_CI_KR_lowers, REML_CI_KR_uppers,
                               true_params_REMLrep, 'REML_coverage')

  # Calculate 'power'
  MCMC_theta_pow <- sum(MCMC_powvals$theta)/length(MCMC_powvals$theta)
  MCMC_pow <- c(MCMC_theta_pow, rep(NA, length(rnames_MCMC)-1))
  MCMC_pow_df <- as.data.frame(t(MCMC_pow))
  colnames(MCMC_pow_df) <- rnames_MCMC
  MCMC_pow_df$measure <- 'MCMC_power'

  # Include BPICC
  true_params$measure <- 'true_val'
  
  measures <- rbind(
    true_params,
    MCMC_mean_bias,
    MCMC_median_bias,
    REML_bias,
    MCMC_avg_SD,
    REML_avg_KR_SE,
    REML_avg_SE,
    MCMC_mean_mse,
    REML_mse,
    MCMC_mean_MCsd,
    REML_MCsd,
    MCMC_coverage,
    REML_KR_coverage,
    MCMC_pow_df
  )
  
  reps <- data.frame(
    Nreps=Nrep,
    MCMCreps=MCMCreps,
    REMLreps=REMLreps
  )
  
  params <- data.frame(
    S=clust_per_seq,
    Tp=periods,
    m=subjects,
    rho1=WPICC,
    r=CAC
  )

  results <- list(
    measures = measures,
    reps = reps,
    params = params
  )
  
  dir.create('performance_measures')

  save(
    'results',
    file=paste0(
      'performance_measures/',
      'performance_measures_',
      config,
      '.Rda')
  )
}

reduce_all_results <- function(Nrep, clust_per_seq, periods, subjects, WPICC, CAC, theta) {
  # Reduce all results files for a particular trial and parameter configuration
  
  for (n in 1:Nrep) {
    reduce_nth_results(n, clust_per_seq, periods, subjects, WPICC, CAC, theta)
  }
}

# Retrieve results for nth replicate
reduce_nth_results <- function(n, clust_per_seq, periods, subjects, WPICC, CAC, theta) {
  # Get estimates of interest from nth result file and save 'reduced_results' file
  
  WPICC100 <- WPICC*100
  CAC100 <- CAC*100
  config <- paste0(
    '_S', clust_per_seq,
    '_T', periods,
    '_m', subjects,
    '_WPICC', WPICC100,
    '_CAC', CAC100,
    '_theta', theta,
    sprintf('_seed%04d', n)
  )
  
  load(file=paste0('results/MCMC', config, '.Rda'))
  # Contains: fit (MCMC), parvals (true values)
  load(file=paste0('results/REML', config, '.Rda'))
  # Contains: remlfit (REML), parvals (true values)
  
  # Get number of divergent transitions
  div <- 0
  chains <- dim(fit)[2] # Extract number of chains from fit
  for (ch in 1:chains) {
    divergent <- get_sampler_params(fit, inc_warmup=FALSE)[[ch]][,'divergent__']
    div <- div + sum(divergent)
  }

  # Extract inference
  
  # Extract estimates from REML model
  REML_theta_est <- coef(summary(remlfit))['treat','Estimate']
  # Point estimate only: fixef(remlfit)['treat']
  REML_theta_stderr <- coef(summary(remlfit))['treat','Std. Error']
  
  # Calculate 95% confidence interval for theta using Kenward Roger correction
  confint <- ci_kenward(remlfit, ci=0.95)
  theta_CI_KR_low <- confint[confint$Parameter=='treat','CI_low']
  theta_CI_KR_high <- confint[confint$Parameter=='treat','CI_high']
  
  # Get KR-adjusted standard error for theta
  adj_SE <- se_kenward(remlfit)
  theta_SE_KR <- adj_SE[adj_SE$Parameter=='treat','SE']

  # Combine post-warmup posterior draws across chains
  if (CAC==1.0) {
    varcomps <- as.data.frame(VarCorr(remlfit))
    REML_sig_sq_c <- varcomps[1, 'vcov']
    REML_sig_sq_e <- varcomps[2, 'vcov']
    REML_WPICC <- REML_sig_sq_c / (REML_sig_sq_c + REML_sig_sq_e)
    
    est <- c(REML_theta_est, REML_WPICC, REML_sig_sq_e, REML_sig_sq_c)
    
    # Flag if cluster variance estimated as 0
    if (REML_sig_sq_c==0) {
      zerovar <- TRUE
    } else {
      zerovar <- FALSE
    }
    
    # Extract post-warmup MCMC draws for relevant parameters
    pars <- c('theta', 'WPICC', 'sig_sq_subject', 'sig_sq_cluster')
    post <- as.data.frame(fit, pars=pars)
  } else {
    varcomps <- as.data.frame(VarCorr(remlfit))
    REML_sig_sq_cp <- varcomps[1, 'vcov']
    REML_sig_sq_c <- varcomps[2, 'vcov']
    REML_sig_sq_e <- varcomps[3, 'vcov']
    sum_c_cp <- (REML_sig_sq_c + REML_sig_sq_cp)
    REML_WPICC <- sum_c_cp / (sum_c_cp + REML_sig_sq_e)
    REML_CAC <- REML_sig_sq_c / sum_c_cp
    if (is.na(REML_CAC)) { # TODO: Decide how best to handle this
      REML_CAC <- 0
    }
    REML_BPICC <- REML_WPICC * REML_CAC
    
    est=c(REML_theta_est, REML_WPICC, REML_CAC, REML_sig_sq_e,
          REML_sig_sq_c, REML_sig_sq_cp, REML_BPICC)
    
    # Flag if cluster variance estimated as 0
    if (REML_sig_sq_c==0 | REML_sig_sq_cp==0) {
      zerovar <- TRUE
    } else {
      zerovar <- FALSE
    }
    
    # Extract post-warmup MCMC draws for relevant parameters
    pars <- c('theta', 'WPICC', 'CAC', 'sig_sq_subject', 'sig_sq_cluster', 'sig_sq_cp')
    post <- as.data.frame(fit, pars=pars)
    
    # Derive BPICC posterior from posterior draws
    post$BPICC <- post$WPICC * post$CAC
    pars <- c('theta', 'WPICC', 'CAC', 'sig_sq_subject', 'sig_sq_cluster', 'sig_sq_cp', 'BPICC')
  }
  
  REML_res <- data.frame(
    est=est,
    stderr=c(REML_theta_stderr, rep(NA, (length(pars)-1))),
    stderrKR=c(theta_SE_KR, rep(NA, length(pars)-1)),
    CI_KR_lower=c(theta_CI_KR_low, rep(NA, length(pars)-1)),
    CI_KR_upper=c(theta_CI_KR_high, rep(NA, length(pars)-1)),
    row.names=pars
  )

  # Posterior means
  MCMC_means <- apply(post, 2, mean)

  # Posterior medians
  MCMC_medians <- apply(post, 2, quantile, probs=0.5)
  
  # Standard deviation of marginal posterior distributions
  MCMC_sds <- apply(post, 2, sd)
  
  # 95% credible intervals from 2.5% and 97.5% percentiles of posterior draws
  MCMC_credints <- apply(post, 2, quantile, probs=c(0.025, 0.975))
  
  # Posterior probability: P(theta > 0)
  # Calculate proportion of posterior draws for theta that are greater than 0
  prop <- sum(post$theta > 0)/length(post$theta)
  # Is this probability greater than cutoff C=0.975?
  MCMC_greaterthanC <- (prop > 0.975)
  
  MCMC_res <- data.frame(
    post_mean=MCMC_means,
    post_median=MCMC_medians,
    post_sd=MCMC_sds,
    post_025=MCMC_credints['2.5%',],
    post_975=MCMC_credints['97.5%',],
    theta_greaterthanC=c(MCMC_greaterthanC, rep(NA, (length(pars)-1)))
  )
  
  # Retrieve true parameter values
  parvals <- c(parvals, BPICC=WPICC*CAC)
  truevals <- as.data.frame(parvals[pars])
  
  res <- list(
    MCMC_res = MCMC_res,
    REML_res = REML_res,
    truevals = truevals,
    div = div,
    zerovar = zerovar
  )
  
  dir.create('reduced_results')
  
  save(
    'res',
    file=paste0(
      'reduced_results/',
      'reduced_results',
      config,
      '.Rda'
    )
  )
}

# Compute bias
bias <- function(sim_estimates, true_params, measure_name) {
  # Difference between average estimate across replications and true value
  bias_sim_estimates_df <- colMeans(sim_estimates) - true_params
  bias_sim_estimates_df$measure <- measure_name
  return(bias_sim_estimates_df)
}

# Compute MSE
mse <- function(sim_estimates, true_params, measure_name) {
  sim_estimates_centered <- sweep(as.matrix(sim_estimates), 2, as.matrix(true_params), '-')
  mse_sim_estimates <- apply(sim_estimates_centered^2, 2, mean)
  mse_sim_estimates_df <- as.data.frame(t(mse_sim_estimates)) # TODO: Find more elegant way
  mse_sim_estimates_df$measure <- measure_name
  return(mse_sim_estimates_df)
}

# Compute Monte Carlo standard deviation
MC.sd <- function(sim_estimates, measure_name) {
  MC_sd_sim_estimates <- apply(as.matrix(sim_estimates), 2, sd)
  MC_sd_sim_estimates_df <- as.data.frame(t(MC_sd_sim_estimates))
  MC_sd_sim_estimates_df$measure <- measure_name
  return(MC_sd_sim_estimates_df)
}

# Do confidence/credible intervals include true parameter value?
coverage <- function(sim_lower, sim_upper, true_params_rep, measure_name) {
  Nrep <- dim(sim_lower)[1]
  covered <- (true_params_rep >= sim_lower) & (true_params_rep <= sim_upper)
  covprop <- colSums(covered)/Nrep
  covprop_df <- as.data.frame(t(covprop))
  covprop_df$measure <- measure_name
  return(covprop_df)
}
