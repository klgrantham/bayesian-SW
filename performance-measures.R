# Calculate performance measures for models with Bayesian and REML estimation
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(rstan)
library(lme4)
library(tidyverse)

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
  # Contains: est$REML_ests, est$REML_stderrs, est$REML_stderrKRs,
  #           est$REML_CI_lowers, est$REML_CI_uppers,
  #           est$REML_CI_KR_lowers, est$REML_CI_KR_uppers,
  #           est$MCMC_means, est$MCMC_medians, est$MCMC_sds,
  #           est$MCMC_025, est$MCMC_975, est$MCMC_powvals,
  #           est$true_params, est$nonzerovar
  
  # Calculate performance measures across all valid replicates
  
  # Drop MCMC replicates with divergent transitions
  MCMC_means <- est$MCMC_means %>% filter(div==0) %>% select(-div)
  MCMC_medians <- est$MCMC_medians %>% filter(div==0) %>% select(-div)
  MCMC_sds <- est$MCMC_sds %>% filter(div==0) %>% select(-div)
  MCMC_025 <- est$MCMC_025 %>% filter(div==0) %>% select(-div)
  MCMC_975 <- est$MCMC_975 %>% filter(div==0) %>% select(-div)
  MCMC_powvals <- est$MCMC_powvals %>% filter(div==0) %>% select(-div)
  
  # Number of valid MCMC estimates (those with divergent transitions excluded)
  MCMCreps <- dim(MCMC_means)[1]
  
  # Number of REML estimates
  REMLreps <- dim(est$REML_ests)[1]
  
  # Drop REML replicates where CAC is NA from calculation of performance
  # measures for CAC only; retain replicates for other parameters
  if(CAC != 1.0) {
    REML_ests_CAC <- est$REML_ests %>%
      filter(!is.na(est$REML_ests$CAC)) %>%
      select(CAC, BPICC)
    REML_ests_exclCAC <- est$REML_ests %>%
      select(-c(CAC, BPICC))
    true_params_CAC <- est$true_params %>%
      select(c(CAC, BPICC))
    true_params_exclCAC <- est$true_params %>%
      select(-c(CAC, BPICC))
    col_order <- colnames(est$REML_ests)
    
    # Number of valid REML estimates for CAC
    REMLreps <- dim(REML_ests_CAC)[1]
  
    # Calculate bias
    REML_bias_CAC <- bias(REML_ests_CAC, true_params_CAC, 'bias', 'REML')
    REML_bias_exclCAC <- bias(REML_ests_exclCAC, true_params_exclCAC, 'bias', 'REML')
    REML_bias_join <- left_join(REML_bias_exclCAC, REML_bias_CAC, by=c('measure', 'method'))
    REML_bias <- REML_bias_join[, c(col_order, 'measure', 'method')]
    # Calculate MCSE of bias estimates
    REML_MCSE_bias_CAC <- MCSE_bias(REML_ests_CAC, 'MCSE_bias', 'REML')
    REML_MCSE_bias_exclCAC <- MCSE_bias(REML_ests_exclCAC, 'MCSE_bias', 'REML')
    REML_MCSE_bias_join <- left_join(REML_MCSE_bias_exclCAC, REML_MCSE_bias_CAC,
                                     by=c('measure', 'method'))
    REML_MCSE_bias <- REML_MCSE_bias_join[, c(col_order, 'measure', 'method')]
    
    # Calculate MSE
    REML_MSE_CAC <- MSE(REML_ests_CAC, true_params_CAC, 'MSE', 'REML')
    REML_MSE_exclCAC <- MSE(REML_ests_exclCAC, true_params_exclCAC, 'MSE', 'REML')
    REML_MSE_join <- left_join(REML_MSE_exclCAC, REML_MSE_CAC, by=c('measure', 'method'))
    REML_MSE <- REML_MSE_join[, c(col_order, 'measure', 'method')]
    # Calculate MCSE of bias estimates
    REML_MCSE_MSE_CAC <- MCSE_MSE(REML_ests_CAC, true_params_CAC, REML_MSE_CAC,
                                  'MCSE_MSE', 'REML')
    REML_MCSE_MSE_exclCAC <- MCSE_MSE(REML_ests_exclCAC, true_params_exclCAC,
                                      REML_MSE_exclCAC, 'MCSE_MSE', 'REML')
    REML_MCSE_MSE_join <- left_join(REML_MCSE_MSE_exclCAC, REML_MCSE_MSE_CAC,
                                     by=c('measure', 'method'))
    REML_MCSE_MSE <- REML_MCSE_MSE_join[, c(col_order, 'measure', 'method')]
  } else {
    # Number of valid REML estimates
    REMLreps <- dim(est$REML_ests)[1]
    
    # Calculate bias
    REML_bias <- bias(est$REML_ests, est$true_params, 'bias', 'REML')
    # Calculate MCSE of bias estimates
    REML_MCSE_bias <- MCSE_bias(est$REML_ests, 'MCSE_bias', 'REML')
    
    # Calculate MSE
    REML_MSE <- MSE(est$REML_ests, est$true_params, 'MSE', 'REML')
    # Calculate MCSE of MSE estimates
    REML_MCSE_MSE <- MCSE_MSE(est$REML_ests, est$true_params, REML_MSE,
                              'MCSE_MSE', 'REML')
  }
  
  # Calculate bias
  MCMC_median_bias <- bias(MCMC_medians, est$true_params, 'bias', 'MCMC')
  # Calculate MCSE of bias estimates
  MCMC_median_MCSE_bias <- MCSE_bias(MCMC_medians, 'MCSE_bias', 'MCMC')

  # Calculate MSE
  MCMC_median_MSE <- MSE(MCMC_medians, est$true_params, 'MSE', 'MCMC')
  # Calculate MCSE of MSE estimates
  MCMC_median_MCSE_MSE <- MCSE_MSE(MCMC_medians, est$true_params,
                                   MCMC_median_MSE, 'MCSE_MSE', 'MCMC')

  # Calculate interval coverage
  true_params_MCMCrep <- est$true_params[rep(1, dim(MCMC_025)[1]),] # true_params is a row vector
  MCMC_coverage <- coverage(MCMC_025, MCMC_975,
                            true_params_MCMCrep, 'coverage', 'MCMC')
  true_params_REMLrep <- est$true_params[rep(1, dim(est$REML_CI_lowers)[1]),]
  REML_coverage <- coverage(est$REML_CI_lowers, est$REML_CI_uppers,
                            true_params_REMLrep, 'coverage', 'REML')
  REML_KR_coverage <- coverage(est$REML_CI_KR_lowers, est$REML_CI_KR_uppers,
                               true_params_REMLrep, 'coverage', 'REML (KR)')
  # Calculate MCSE of coverage estimates
  MCMC_MCSE_coverage <- MCSE_coverage(MCMCreps, MCMC_coverage,
                                      'MCSE_coverage', 'MCMC')
  REML_MCSE_coverage <- MCSE_coverage(REMLreps, REML_coverage,
                                      'MCSE_coverage', 'REML')
  REML_KR_MCSE_coverage <- MCSE_coverage(REMLreps, REML_KR_coverage,
                                         'MCSE_coverage', 'REML (KR)')
  
  # Calculate empirical standard error
  MCMC_median_empSE <- empSE(MCMC_medians, 'empSE', 'MCMC')
  REML_empSE <- empSE(est$REML_ests, 'empSE', 'REML')
  # Calculate MCSE of empirical SE
  MCMC_median_MCSE_empSE <- MCSE_empSE(MCMCreps, MCMC_median_empSE, 'MCSE_empSE', 'MCMC')
  REML_MCSE_empSE <- MCSE_empSE(REMLreps, REML_empSE, 'MCSE_empSE', 'REML')
  
  # Calculate average model standard error
  # Root-mean of squared model SEs
  MCMC_avgmodSE <- avgmodSE(MCMC_sds, 'avgmodSE', 'MCMC')
  REML_avgmodSE <- avgmodSE(est$REML_stderrs, 'avgmodSE', 'REML')
  REML_avgKRmodSE <- avgmodSE(est$REML_stderrKRs, 'avgmodSE', 'REML (KR)')
  # Calculate MCSE of average model SE
  MCMC_MCSE_avgmodSE <- MCSE_avgmodSE(MCMC_sds, MCMC_avgmodSE,
                                      'MCSE_avgmodSE', 'MCMC')
  REML_MCSE_avgmodSE <- MCSE_avgmodSE(est$REML_stderrs, REML_avgmodSE,
                                      'MCSE_avgmodSE', 'REML')
  REML_MCSE_avgKRmodSE <- MCSE_avgmodSE(est$REML_stderrKRs, REML_avgKRmodSE,
                                        'MCSE_avgmodSE', 'REML (KR)')

  # Calculate relative % error in model standard error
  MCMC_pcterr_modSE <- pcterr_modSE(MCMC_avgmodSE, MCMC_median_empSE, 'pcterrmodSE', 'MCMC')
  REML_pcterr_modSE <- pcterr_modSE(REML_avgmodSE, REML_empSE, 'pcterrmodSE', 'REML')
  REML_KR_pcterr_modSE <- pcterr_modSE(REML_avgKRmodSE, REML_empSE, 'pcterrmodSE', 'REML (KR)')
  # Calculate MCSE of relative % error in model SE
  MCMC_MCSE_pcterr_modSE <- MCSE_pcterr_modSE(MCMC_sds, MCMC_avgmodSE, MCMC_median_empSE,
                                              'MCSE_pcterrmodSE', 'MCMC')
  REML_MCSE_pcterr_modSE <- MCSE_pcterr_modSE(est$REML_stderrs, REML_avgmodSE, REML_empSE,
                                              'MCSE_pcterrmodSE', 'REML')
  REML_KR_MCSE_pcterr_modSE <- MCSE_pcterr_modSE(est$REML_stderrKRs, REML_avgKRmodSE, REML_empSE,
                                                 'MCSE_pcterrmodSE', 'REML (KR)')
  
  # Calculate average interval length
  MCMC_avgintlen <- avgintlength(MCMC_025, MCMC_975,
                                 'avgintlength', 'MCMC')
  REML_avgintlen <- avgintlength(est$REML_CI_lowers, est$REML_CI_uppers,
                                 'avgintlength', 'REML')
  REML_KR_avgintlen <- avgintlength(est$REML_CI_KR_lowers, est$REML_CI_KR_uppers,
                                    'avgintlength', 'REML (KR)')
  # MCSE of average interval length (placeholders, not calculated)
  MCMC_MCSE_avgintlen <- MCSE_avgintlength(MCMC_avgintlen,
                                           'MCSE_avgintlength', 'MCMC')
  REML_MCSE_avgintlen <- MCSE_avgintlength(REML_avgintlen,
                                           'MCSE_avgintlength', 'REML')
  REML_KR_MCSE_avgintlen <- MCSE_avgintlength(REML_KR_avgintlen,
                                              'MCSE_avgintlength', 'REML (KR)')
  
  # Calculate 'power'
  MCMC_pow <- colMeans(MCMC_powvals)
  MCMC_pow_df <- as.data.frame(t(MCMC_pow))
  MCMC_pow_df$measure <- 'power'
  MCMC_pow_df$method <- 'MCMC'
  
  # Label true parameter values
  est$true_params$measure <- 'true_val'
  est$true_params$method <- 'Both'
  
  measures <- rbind(
    est$true_params,
    MCMC_median_bias,
    REML_bias,
    MCMC_median_MCSE_bias,
    REML_MCSE_bias,
    MCMC_median_MSE,
    REML_MSE,
    MCMC_median_MCSE_MSE,
    REML_MCSE_MSE,
    MCMC_coverage,
    REML_coverage,
    REML_KR_coverage,
    MCMC_MCSE_coverage,
    REML_MCSE_coverage,
    REML_KR_MCSE_coverage,
    MCMC_median_empSE,
    REML_empSE,
    MCMC_median_MCSE_empSE,
    REML_MCSE_empSE,
    MCMC_avgmodSE,
    REML_avgmodSE,
    REML_avgKRmodSE,
    MCMC_MCSE_avgmodSE,
    REML_MCSE_avgmodSE,
    REML_MCSE_avgKRmodSE,
    MCMC_pcterr_modSE,
    REML_pcterr_modSE,
    REML_KR_pcterr_modSE,
    MCMC_MCSE_pcterr_modSE,
    REML_MCSE_pcterr_modSE,
    REML_KR_MCSE_pcterr_modSE,
    MCMC_avgintlen,
    REML_avgintlen,
    REML_KR_avgintlen,
    MCMC_MCSE_avgintlen,
    REML_MCSE_avgintlen,
    REML_KR_MCSE_avgintlen,
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
    r=CAC,
    trt=theta
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
  # Collate results across all replicates and save estimates
  
  dir.create('estimates')
  
  MCMC_means <- data.frame()
  MCMC_medians <- data.frame()
  MCMC_sds <- data.frame()
  MCMC_025 <- data.frame()
  MCMC_975 <- data.frame()
  MCMC_powvals <- data.frame()
  REML_ests <- data.frame()
  REML_stderrs <- data.frame()
  REML_CI_lowers <- data.frame()
  REML_CI_uppers <- data.frame()
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
  
  zerovarreps <- 0
  divreps <- 0
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
    
    # Count REML replicates with at least one variance estimate of 0
    if (isTRUE(res$zerovar)) {
      zerovarreps <- zerovarreps + 1
    }
    
    # Combine REML results across all replicates (including those with invalid CAC)
    # Include nth REML results in results containers
    REML_ests <- rbind(REML_ests, REML$est)
    REML_stderrs <- rbind(REML_stderrs, REML$stderr)
    REML_CI_lowers <- rbind(REML_CI_lowers, REML$CI_lower)
    REML_CI_uppers <- rbind(REML_CI_uppers, REML$CI_upper)
    REML_stderrKRs <- rbind(REML_stderrKRs, REML$stderrKR)
    REML_CI_KR_lowers <- rbind(REML_CI_KR_lowers, REML$CI_KR_lower)
    REML_CI_KR_uppers <- rbind(REML_CI_KR_uppers, REML$CI_KR_upper)
    
    # Count MCMC replicates with at least one divergent transition
    if (res$div > 0) {
      divreps <- divreps + 1
      div <- 1
    } else{
      div <- 0
    }
    
    # Include nth MCMC results in results containers 
    MCMC_means <- rbind(MCMC_means, c(MCMC$post_mean, div))
    MCMC_medians <- rbind(MCMC_medians, c(MCMC$post_median, div))
    MCMC_sds <- rbind(MCMC_sds, c(MCMC$post_sd, div))
    MCMC_025 <- rbind(MCMC_025, c(MCMC$post_025, div))
    MCMC_975 <- rbind(MCMC_975, c(MCMC$post_975, div))
    MCMC_powvals <- rbind(MCMC_powvals, c(MCMC$theta_greaterthanC, div))
  }
  parnames_MCMC <- c(MCMC$parameter, 'div')
  parnames_REML <- REML$parameter
  colnames(REML_ests) <- parnames_REML
  colnames(REML_stderrs) <- parnames_REML
  colnames(REML_CI_lowers) <- parnames_REML
  colnames(REML_CI_uppers) <- parnames_REML
  colnames(REML_stderrKRs) <- parnames_REML
  colnames(REML_CI_KR_lowers) <- parnames_REML
  colnames(REML_CI_KR_uppers) <- parnames_REML
  colnames(MCMC_means) <- parnames_MCMC
  colnames(MCMC_medians) <- parnames_MCMC
  colnames(MCMC_sds) <- parnames_MCMC
  colnames(MCMC_025) <- parnames_MCMC
  colnames(MCMC_975) <- parnames_MCMC
  colnames(MCMC_powvals) <- parnames_MCMC
  
  true_params <- res$truevals
  
  # Save estimates datasets
  est <- list(
    REML_ests = REML_ests,
    REML_stderrs = REML_stderrs,
    REML_CI_lowers = REML_CI_lowers,
    REML_CI_uppers = REML_CI_uppers,
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
    zerovarreps = zerovarreps,
    divreps = divreps
  )
  
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
  
  dir.create('reduced_results')
  
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
  # Contains: draws (MCMC posterior draws), parvals (true values), diagnostics
  load(file=paste0('results/REML', config, '.Rda'))
  # Contains: sumreml (REML fit summary), parvals (true values), adj_ddf (KR-adjusted ddf), adj_SE (KR-adjusted SE)
  
  # Get number of divergent transitions
  div <- diagnostics$ndiv_trans

  # Extract inference
  
  # Extract estimates from REML model
  REML_theta_est <- coef(sumreml)['treat','Estimate']
  REML_theta_stderr <- coef(sumreml)['treat','Std. Error']
  
  # Get 95% confidence interval for theta using Kenward Roger correction
  
  # Get t test statistic with KR-adjusted ddf
  alpha <- (1 + 0.95)/2
  tstat <- qt(alpha, adj_ddf)
  
  # Construct 95% KR confidence intervals
  theta_CI_KR_low <- REML_theta_est - tstat * adj_SE
  theta_CI_KR_high <- REML_theta_est + tstat * adj_SE
  
  # Get KR-adjusted standard error for theta
  theta_SE_KR <- adj_SE

  # Get unadjusted 95% confidence interval for theta
  
  # Get standard normal test statistic
  zstat <- qnorm(alpha, 0, 1)
  
  # Construct 95% KR confidence intervals
  theta_CI_low <- REML_theta_est - zstat * REML_theta_stderr
  theta_CI_high <- REML_theta_est + zstat * REML_theta_stderr
  
  # Combine post-warmup posterior draws across chains
  if (CAC==1.0) {
    varcomps <- as.data.frame(sumreml$varcor)
    REML_sig_sq_c <- varcomps$vcov[varcomps$grp=='clust']
    REML_sig_sq_e <- varcomps$vcov[varcomps$grp=='Residual']
    REML_WPICC <- REML_sig_sq_c / (REML_sig_sq_c + REML_sig_sq_e)
    
    est <- c(REML_theta_est, REML_WPICC, REML_sig_sq_e, REML_sig_sq_c)
    
    # Flag if cluster variance estimated as 0
    if (REML_sig_sq_c==0) {
      zerovar <- TRUE
    } else {
      zerovar <- FALSE
    }
  } else {
    varcomps <- as.data.frame(sumreml$varcor)
    REML_sig_sq_cp <- varcomps$vcov[varcomps$grp=='clustper']
    REML_sig_sq_c <- varcomps$vcov[varcomps$grp=='clust']
    REML_sig_sq_e <- varcomps$vcov[varcomps$grp=='Residual']
    sum_c_cp <- (REML_sig_sq_c + REML_sig_sq_cp)
    REML_WPICC <- sum_c_cp / (sum_c_cp + REML_sig_sq_e)
    REML_CAC <- REML_sig_sq_c / sum_c_cp
    REML_BPICC <- REML_WPICC * REML_CAC
    
    est=c(REML_theta_est, REML_WPICC, REML_CAC, REML_BPICC, REML_sig_sq_e,
          REML_sig_sq_c, REML_sig_sq_cp)
    
    # Flag if cluster or cluster-period variance estimated as 0
    if (REML_sig_sq_c==0 | REML_sig_sq_cp==0) {
      zerovar <- TRUE
    } else {
      zerovar <- FALSE
    }
  }
  
  pars <- colnames(draws)
  
  REML_res <- data.frame(
    est=est,
    stderr=c(REML_theta_stderr, rep(NA, length(pars)-1)),
    CI_lower=c(theta_CI_low, rep(NA, length(pars)-1)),
    CI_upper=c(theta_CI_high, rep(NA, length(pars)-1)),
    stderrKR=c(theta_SE_KR, rep(NA, length(pars)-1)),
    CI_KR_lower=c(theta_CI_KR_low, rep(NA, length(pars)-1)),
    CI_KR_upper=c(theta_CI_KR_high, rep(NA, length(pars)-1)),
    parameter=pars
  )

  # Posterior means
  MCMC_means <- apply(draws, 2, mean)

  # Posterior medians
  MCMC_medians <- apply(draws, 2, quantile, probs=0.5)
  
  # Standard deviation of marginal posterior distributions
  MCMC_sds <- apply(draws, 2, sd)
  
  # 95% credible intervals from 2.5% and 97.5% percentiles of posterior draws
  MCMC_credints <- apply(draws, 2, quantile, probs=c(0.025, 0.975))
  
  # Posterior probability: P(theta > 0)
  # Calculate proportion of posterior draws for theta that are greater than 0
  prop <- sum(draws$theta > 0)/length(draws$theta)
  # Is this probability greater than cutoff C=0.975?
  MCMC_greaterthanC <- (prop > 0.975)
  
  MCMC_res <- data.frame(
    post_mean=MCMC_means,
    post_median=MCMC_medians,
    post_sd=MCMC_sds,
    post_025=MCMC_credints['2.5%',],
    post_975=MCMC_credints['97.5%',],
    theta_greaterthanC=c(MCMC_greaterthanC, rep(NA, length(pars)-1)),
    parameter=pars
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

# Average interval length
avgintlength <- function(sim_lower, sim_upper, measure_name, method_name) {
  intlens <- sim_upper - sim_lower
  avgintlen <- colMeans(intlens)
  avgintlen_df <- as.data.frame(t(avgintlen))
  avgintlen_df$measure <- measure_name
  avgintlen_df$method <- method_name
  return(avgintlen_df)
}

# MCSE of average interval length (placeholder, not calculated)
MCSE_avgintlength <- function(intlen_ests, measure_name, method_name) {
  intlen_ests <- select(intlen_ests, -c('measure', 'method'))
  intlen_ests[1:length(intlen_ests)] <- NA
  MCSE_intlen_df <- as.data.frame(intlen_ests)
  MCSE_intlen_df$measure <- measure_name
  MCSE_intlen_df$method <- method_name
  return(MCSE_intlen_df)
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

pcterr_modSE <- function(avgmodSE_ests, empSE_ests, measure_name, method_name) {
  avgmodSE_ests <- select(avgmodSE_ests, -c('measure', 'method'))
  empSE_ests <- select(empSE_ests, -c('measure', 'method'))
  
  pcterr <- 100 * ((avgmodSE_ests/empSE_ests) - 1)
  pcterr_df <- as.data.frame(pcterr)
  pcterr_df$measure <- measure_name
  pcterr_df$method <- method_name
  return(pcterr_df)
}

MCSE_pcterr_modSE <- function(sim_estimates, avgmodSE_ests, empSE_ests, measure_name, method_name) {
  nsim <- dim(sim_estimates)[1]
  avgmodSE_ests <- select(avgmodSE_ests, -c('measure', 'method'))
  empSE_ests <- select(empSE_ests, -c('measure', 'method'))
  
  quad_avgmodSE_ests <- as.matrix(avgmodSE_ests)^4
  sq_ests <- as.matrix(sim_estimates)^2
  varvar <- apply(sq_ests, 2, var)
  MCSE_pcterr_modSE_a <- 100 * (avgmodSE_ests/empSE_ests)
  MCSE_pcterr_modSE_b <- sqrt((varvar/(4*nsim*(quad_avgmodSE_ests))) + (1/(2*(nsim-1))))
  MCSE_pcterr_modSE <- MCSE_pcterr_modSE_a * MCSE_pcterr_modSE_b
  MCSE_pcterr_modSE_df <- as.data.frame(MCSE_pcterr_modSE)
  MCSE_pcterr_modSE_df$measure <- measure_name
  MCSE_pcterr_modSE_df$method <- method_name
  return(MCSE_pcterr_modSE_df)
}
