# Calculate performance measures for models with MCMC and REML estimation
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(rstan)
library(lme4)
library(tidyverse)
library(parameters)

calculate_measures <- function(clust_per_seq, periods, subjects, WPICC, CAC, theta) {
  # Calculate performance measures with estimates datasets and save results
  
  # Load estimates datasets
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
  
  infile=paste0(
    'estimates/',
    'estimates',
    config,
    '.Rda'
  )
  load(file=infile)
  # Contains: est$REML_ests, est$REML_stderrs, est$stderrKRs,
  #           est$REML_CI_KR_lowers, est$REML_CI_KR_uppers,
  #           est$MCMC_means, est$MCMC_medians, est$MCMC_sds,
  #           est$MCMC_025, est$MCMC_975, est$MCMC_powvals,
  #           est$true_params, est$nonzerovar
  
  # Calculate performance measures across all valid replicates
  
  # Number of valid MCMC estimates (those with divergent transitions excluded)
  MCMCreps <- dim(est$MCMC_means)[1]
  
  # Number of REML estimates
  REMLreps <- dim(est$REML_ests)[1]

  # Calculate bias
  MCMC_mean_bias <- bias(est$MCMC_means, est$true_params, 'bias', 'MCMC')
  REML_bias <- bias(est$REML_ests, est$true_params, 'bias', 'REML')
  # Calculate MCSE of bias estimates
  MCMC_mean_MCSE_bias <- MCSE_bias(est$MCMC_means, 'MCSE_bias', 'MCMC')
  REML_MCSE_bias <- MCSE_bias(est$REML_ests, 'MCSE_bias', 'REML')
  
  # Calculate MSE
  MCMC_mean_MSE <- MSE(est$MCMC_means, est$true_params, 'MSE', 'MCMC')
  REML_MSE <- MSE(est$REML_ests, est$true_params, 'MSE', 'REML')
  # Calculate MCSE of MSE estimates
  MCMC_mean_MCSE_MSE <- MCSE_MSE(est$MCMC_means, est$true_params,
                                 MCMC_mean_MSE, 'MCSE_MSE', 'MCMC')
  REML_MCSE_MSE <- MCSE_MSE(est$REML_ests, est$true_params, REML_MSE,
                            'MCSE_MSE', 'REML')
  
  # Calculate interval coverage
  true_params_MCMCrep <- est$true_params[rep(1, MCMCreps),] # true_params is a row vector
  MCMC_coverage <- coverage(est$MCMC_025, est$MCMC_975,
                            true_params_MCMCrep, 'coverage', 'MCMC')
  true_params_REMLrep <- est$true_params[rep(1, REMLreps),]
  REML_KR_coverage <- coverage(est$REML_CI_KR_lowers, est$REML_CI_KR_uppers,
                               true_params_REMLrep, 'coverage', 'REML')
  # Calculate MCSE of coverage estimates
  MCMC_MCSE_coverage <- MCSE_coverage(MCMCreps, MCMC_coverage,
                                      'MCSE_coverage', 'MCMC')
  REML_KR_MCSE_coverage <- MCSE_coverage(REMLreps, REML_KR_coverage,
                                         'MCSE_coverage', 'REML')
  
  # Calculate empirical standard error
  MCMC_mean_empSE <- empSE(est$MCMC_means, 'empSE', 'MCMC')
  REML_empSE <- empSE(est$REML_ests, 'empSE', 'REML')
  # Calculate MCSE of empirical SE
  MCMC_mean_MCSE_empSE <- MCSE_empSE(MCMCreps, MCMC_mean_empSE, 'MCSE_empSE', 'MCMC')
  REML_MCSE_empSE <- MCSE_empSE(REMLreps, REML_empSE, 'MCSE_empSE', 'REML')
  
  # Calculate average model standard error
  # Root-mean of squared model SEs
  MCMC_avgmodSE <- avgmodSE(est$MCMC_sds, 'avgmodSE', 'MCMC')
  REML_avgKRmodSE <- avgmodSE(est$REML_stderrKRs, 'avgmodSE', 'REML')
  # Calculate MCSE of average model SE
  MCMC_MCSE_avgmodSE <- MCSE_avgmodSE(est$MCMC_sds, MCMC_avgmodSE,
                                      'MCSE_avgmodSE', 'MCMC')
  REML_MCSE_avgKRmodSE <- MCSE_avgmodSE(est$REML_stderrKRs, REML_avgKRmodSE,
                                        'MCSE_avgmodSE', 'REML')

  # Calculate 'power'
  MCMC_pow <- colMeans(est$MCMC_powvals)
  MCMC_pow_df <- as.data.frame(t(MCMC_pow))
  MCMC_pow_df$measure <- 'power'
  MCMC_pow_df$method <- 'MCMC'
  
  # Label true parameter values
  est$true_params$measure <- 'true_val'
  est$true_params$method <- 'Both'
  
  measures <- rbind(
    est$true_params,
    MCMC_mean_bias,
    REML_bias,
    MCMC_mean_MCSE_bias,
    REML_MCSE_bias,
    MCMC_mean_MSE,
    REML_MSE,
    MCMC_mean_MCSE_MSE,
    REML_MCSE_MSE,
    MCMC_coverage,
    REML_KR_coverage,
    MCMC_MCSE_coverage,
    REML_KR_MCSE_coverage,
    MCMC_mean_empSE,
    REML_empSE,
    MCMC_mean_MCSE_empSE,
    REML_MCSE_empSE,
    MCMC_avgmodSE,
    REML_avgKRmodSE,
    MCMC_MCSE_avgmodSE,
    REML_MCSE_avgKRmodSE,
    MCMC_pow_df
  )
  
  reps <- data.frame(
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
      'performancemeasures',
      config,
      '.Rda')
  )
}

collate_results <- function(Nrep, clust_per_seq, periods, subjects, WPICC, CAC, theta) {
  # Collate results across (valid) replicates and save estimates
  
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
      'reducedresults',
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
  
  # Save estimates datasets
  est <- list(
    REML_ests = REML_ests,
    REML_stderrs = REML_stderrs,
    REML_stderrKRs = REML_stderrKRs,
    REML_CI_KR_lowers = REML_CI_KR_lowers,
    REML_CI_KR_uppers = REML_CI_KR_uppers,
    MCMC_means = MCMC_means,
    MCMC_medians = MCMC_medians,
    MCMC_sds = MCMC_sds,
    MCMC_025 = MCMC_025,
    MCMC_975 = MCMC_975,
    MCMC_powvals = MCMC_powvals,
    true_params = true_params,
    nonzerovar = nonzerovar
  )
  
  dir.create('estimates')

  save(
    'est',
    file=paste0(
      'estimates/',
      'estimates',
      config,
      '.Rda'
    )
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
  # Contains: fit (MCMC), parvals (true values), diagnostics
  load(file=paste0('results/REML', config, '.Rda'))
  # Contains: remlfit (REML), parvals (true values)
  
  # Get number of divergent transitions
  div <- diagnostics$ndiv_trans

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
    pars <- c('theta', 'WPICC', 'CAC', 'BPICC', 'sig_sq_subject', 'sig_sq_cluster', 'sig_sq_cp')
    post <- as.data.frame(fit, pars=pars)
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
  truevals <- parvals[pars]
  
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
      'reducedresults',
      config,
      '.Rda'
    )
  )
}

# Bias
bias <- function(sim_estimates, true_params, measure_name, method_name) {
  # Difference between average estimate across replications and true value
  bias_df <- colMeans(sim_estimates) - true_params
  bias_df$measure <- measure_name
  bias_df$method <- method_name
  return(bias_df)
}

# Monte Carlo standard error (MCSE) of bias estimate
MCSE_bias <- function(sim_estimates, measure_name, method_name) {
  nsim <- dim(sim_estimates)[1]
  var_est <- apply(as.matrix(sim_estimates), 2, var)
  MCSE_bias <- sqrt((1/nsim)*var_est)
  MCSE_bias_df <- as.data.frame(t(MCSE_bias))
  MCSE_bias_df$measure <- measure_name
  MCSE_bias_df$method <- method_name
  return(MCSE_bias_df)
}

# Mean squared error (MSE)
MSE <- function(sim_estimates, true_params, measure_name, method_name) {
  diff_est <- sweep(as.matrix(sim_estimates), 2, as.matrix(true_params), '-')
  mean_sq_diff <- apply(diff_est^2, 2, mean)
  MSE_df <- as.data.frame(t(mean_sq_diff))
  MSE_df$measure <- measure_name
  MSE_df$method <- method_name
  return(MSE_df)
}

# MCSE of MSE estimate
MCSE_MSE <- function(sim_estimates, true_params, MSE_ests, measure_name, method_name) {
  nsim <- dim(sim_estimates)[1]
  diff_est_true <- sweep(as.matrix(sim_estimates), 2, as.matrix(true_params), '-')
  sq_diff <- diff_est_true^2
  MSE_ests <- select(MSE_ests, -c('measure', 'method'))
  diff_sqdiff_MSE <- sweep(as.matrix(sq_diff), 2, as.matrix(MSE_ests), '-')
  mean_sq_diff <- apply(diff_sqdiff_MSE^2, 2, mean)
  MCSE_MSE <- sqrt((1/(nsim-1))*mean_sq_diff)
  MCSE_MSE_df <- as.data.frame(t(MCSE_MSE))
  MCSE_MSE_df$measure <- measure_name
  MCSE_MSE_df$method <- method_name
  return(MCSE_MSE_df)
}

# Interval coverage
# Do confidence/credible intervals include true parameter value?
coverage <- function(sim_lower, sim_upper, true_params_rep, measure_name, method_name) {
  covered <- (true_params_rep >= sim_lower) & (true_params_rep <= sim_upper)
  covprop <- colMeans(covered)
  covprop_df <- as.data.frame(t(covprop))
  covprop_df$measure <- measure_name
  covprop_df$method <- method_name
  return(covprop_df)
}

# MCSE of coverage
MCSE_coverage <- function(nsim, cov_ests, measure_name, method_name) {
  coverage_ests <- select(cov_ests, -c('measure', 'method'))
  MCSE_cov <- sqrt((1/nsim)*(as.matrix(coverage_ests) * (1 - as.matrix(coverage_ests))))
  MCSE_cov_df <- as.data.frame(MCSE_cov)
  MCSE_cov_df$measure <- measure_name
  MCSE_cov_df$method <- method_name
  return(MCSE_cov_df)
}

# Empirical standard error
empSE <- function(sim_estimates, measure_name, method_name) {
  empSE_est <- apply(as.matrix(sim_estimates), 2, sd)
  empSE_df <- as.data.frame(t(empSE_est))
  empSE_df$measure <- measure_name
  empSE_df$method <- method_name
  return(empSE_df)
}

# MCSE of empirical SE
MCSE_empSE <- function(nsim, empSE_ests, measure_name, method_name) {
  empSE_ests <- select(empSE_ests, -c('measure', 'method'))
  MCSE_empSE <- as.matrix(empSE_ests)/(sqrt(2*(nsim-1)))
  MCSE_empSE_df <- as.data.frame(MCSE_empSE)
  MCSE_empSE_df$measure <- measure_name
  MCSE_empSE_df$method <- method_name
  return(MCSE_empSE_df)
}

# Average model standard error
avgmodSE <- function(sim_estimates, measure_name, method_name) {
  avgmodSE <- sqrt(colMeans(sim_estimates^2))
  avgmodSE_df <- as.data.frame(t(avgmodSE))
  avgmodSE_df$measure <- measure_name
  avgmodSE_df$method <- method_name
  return(avgmodSE_df)
}

# MCSE of average model SE
MCSE_avgmodSE <- function(sim_estimates, avgmodSE_ests, measure_name, method_name) {
  nsim <- dim(sim_estimates)[1]
  avgmodSE_ests <- select(avgmodSE_ests, -c('measure', 'method'))
  
  sq_avgmodSE_ests <- as.matrix(avgmodSE_ests)^2
  sq_ests <- as.matrix(sim_estimates)^2
  varvar <- apply(sq_ests, 2, var)
  MCSE_avgmodSE <- sqrt(varvar/(4*nsim*(sq_avgmodSE_ests)))
  MCSE_avgmodSE_df <- as.data.frame(MCSE_avgmodSE)
  MCSE_avgmodSE_df$measure <- measure_name
  MCSE_avgmodSE_df$method <- method_name
  return(MCSE_avgmodSE_df)
}
