# Main procedure: data generation and model fitting
# with both MCMC and REML estimation
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(rstan)
library(lme4)
library(pbkrtest)
library(parameters)
library(lmerTest)

dir.create('results')

fit_models <- function(params, stanmod, seedval){
  # Inputs:
  #   params - a data frame containing input parameters for one particular
  #            simulation configuration
  #   stanmod - compiled Stan model
  #   seedval - random number seed
  # Outputs:
  #   MCMC_[config].Rda - saves an R data file containing the MCMC posterior
  #                       draws (draws), true parameter values (parvals), and
  #                       diagnostics (results of diagnostic checks and number
  #                       of divergent transitions)
  #   REML_[config].Rda - saves an R data file containing the summary of the
  #                       fitted REML model (sumreml), true parameter values
  #                       (parvals), adj_SE and adj_ddf (KR-adjusted for all
  #                       but largest configuration)
  
  set.seed(seedval)
  
  # ---- simulate data ----
  
  start.time <- Sys.time()
  
  # fixed params
  
  clust_per_seq <- params$clust_per_seq     # number of clusters per treatment sequence
  periods <- params$periods                 # number of periods
  subjects <- params$subjects               # subjects per cluster-period
  
  sig_sq_subject <- 1                       # between-subject variance
  WPICC <- params$WPICC                     # within-period intracluster correlation
  CAC <- params$CAC                         # cluster autocorrelation
  
  theta <- params$theta                     # treatment effect
  per <- seq(1, periods)/periods            # period effects
  
  # stepped wedge design matrix
  # Note: clusters must be a multiple of (periods-1)
  X <- matrix(data=0, nrow=periods-1, ncol=periods)
  for(i in 1:(periods-1)){
    X[i,(i+1):periods] <- 1
  }
  X <- X %x% rep(1, clust_per_seq)
  
  # implied between-cluster variance
  sig_sq_cluster <- CAC * sig_sq_subject * WPICC / (1 - WPICC)
  
  # implied between-cluster-period variance
  sig_sq_cp <- sig_sq_subject * WPICC / (1 - WPICC) - sig_sq_cluster
  
  # implied between-period intracluster correlation
  BPICC <- WPICC * CAC
  
  # implied number of clusters
  clusters <- clust_per_seq*(periods-1)
  
  # random variables
  
  C <- rnorm(n=clusters, mean=0, sd=sqrt(sig_sq_cluster))
  CP <- rnorm(n=clusters * periods, mean=0, sd=sqrt(sig_sq_cp))
  Y <- rep(NA, clusters * periods * subjects)
  e <- rnorm(n=clusters * periods * subjects, mean = 0, sd = sqrt(sig_sq_subject))
  for (i in 1:clusters) {
    for (j in 1:periods) {
      q <- (i - 1) * periods + j
      for (k in 1:subjects) {
        # l is the index into the data
        l <- (i - 1) * periods * subjects + (j - 1) * subjects + k
        Y[l] <- per[j] + X[i,j] * theta + C[i] + CP[q] + e[l]
      }
    }
  }
  
  parvals <- data.frame(
    clust_per_seq = clust_per_seq,
    periods = periods,
    subjects = subjects,
    theta = theta,
    WPICC = WPICC,
    CAC = CAC,
    BPICC = BPICC,
    sig_sq_subject = sig_sq_subject,
    sig_sq_cluster = sig_sq_cluster,
    sig_sq_cp = sig_sq_cp
  )
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('Simulated data took:', format(time.taken, digits=4)))
  
  # ---- fit models ----
  
  if (CAC==1.0) {
    pars <- c('theta', 'WPICC', 'sig_sq_subject', 'sig_sq_cluster')
  } else {
    pars <- c('theta', 'WPICC', 'CAC', 'BPICC', 'sig_sq_subject', 'sig_sq_cluster', 'sig_sq_cp')
  }
  
  # Bayesian model
  
  start.time <- Sys.time()
  
  data <- list(
    clusters = clusters,
    subjects = subjects,
    periods = periods,
    N = clusters*subjects*periods,
    Q = clusters*periods,
    X = X,
    Y = Y
  )
  
  fit <- sampling(stanmod, data = data, pars = pars, iter = 6e3, warmup = 1e3,
                  seed=seedval, control = list(adapt_delta = 0.95), cores=1)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('MCMC sampling took:', format(time.taken, digits=4)))
  
  # GLMM
  
  start.time <- Sys.time()
    
  Xvec <- as.vector(kronecker(as.vector(t(X)), matrix(rep(1, subjects), ncol=1)))
  Cvec <- as.factor(matrix(rep(1:clusters, each=periods*subjects), ncol=1))
  CPvec <- as.factor(matrix(rep(1:(clusters*periods), each=subjects), ncol=1))
  Pvec <- as.factor(matrix(rep(1:periods, each=subjects, times=clusters), ncol=1))
  
  dat <- data.frame(
    Y = Y,
    treat = Xvec,
    clust = Cvec,
    per = Pvec,
    clustper = CPvec
  )
  
  if (CAC==1.0) {
    remlfit <- lmer(Y ~ treat + per + (1 | clust), data=dat, REML=TRUE)
  } else {
    remlfit <- lmer(Y ~ treat + per + (1 | clust) + (1 | clustper), data=dat, REML=TRUE)
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('REML fit took:', format(time.taken, digits=4)))
  
  # Diagnostics
  
  start.time <- Sys.time()
  
  draws <- as.array(fit, pars=pars)
  mon <- monitor(draws, warmup=0, print=FALSE)
  vals <- mon[1:(dim(mon)[1]), c('Rhat', 'Bulk_ESS', 'Tail_ESS')]
  print(vals)
  convergence_check <- NULL
  if (any(vals$Rhat > 1.01) | any(vals$Bulk_ESS < 400) | any(vals$Tail_ESS < 400)){
    convergence_check <- 'FAIL'
    print('Suspected convergence issue. Needs further investigation')
  } else{
    convergence_check <- 'PASS'
  }
  
  div <- 0
  chains <- dim(fit)[2] # Extract number of chains from fit
  for (ch in 1:chains) {
    divergent <- get_sampler_params(fit, inc_warmup=FALSE)[[ch]][,'divergent__']
    div <- div + sum(divergent)
  }
  div_trans_check <- NULL
  if (div > 0){
    div_trans_check <- 'FAIL'
  } else{
    div_trans_check <- 'PASS'
  }
  print(paste('Divergent transitions:', div))
  
  diagnostics <- data.frame(convergence_check=convergence_check,
                            div_trans_check=div_trans_check,
                            ndiv_trans=div)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('Diagnostics checks took:', format(time.taken, digits=4)))

  # Inference
  
  start.time <- Sys.time()
  
  print(fit, pars=pars)
  
  # Get posterior draws for parameters of interest
  draws <- as.data.frame(fit, pars=pars)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('MCMC inference took:', format(time.taken, digits=4)))

  start.time <- Sys.time()
  
  sumreml <- summary(remlfit)
  print(sumreml)
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('REML inference took:', format(time.taken, digits=4)))
  
  
  if (clusters*periods*subjects>15000) {
    print("Skipping KR adjustment")
    
    adj_SE <- coef(sumreml)['treat','Std. Error']
    
    # Get Satterthwaite adjusted ddf using parameters (and lmerTest) packages
    # Note: Could not resolve error when using pbkrtest's SATmodcomp() 
    adj_ddf <- dof_satterthwaite(remlfit)['treat'][[1]]
  } else {
    start.time <- Sys.time()
    
    # Get KR adjustments with pbkrtest package
    vcov0 <- vcov(remlfit)
    vcovKR <- vcovAdj(remlfit)
    adj_SE <- sqrt(vcovKR['treat','treat'])
    
    L <- rep(0, length(fixef(remlfit)))
    L[which(names(fixef(remlfit))=='treat')] <- 1
    adj_ddf <- Lb_ddf(L, vcov0, vcovKR)
    
    end.time <- Sys.time()
    time.taken <- end.time - start.time
    print(paste('KR adjustment (pbkrtest) took:', format(time.taken, digits=4)))
  }
  
  # Save only essential data and results
  
  start.time <- Sys.time()
  
  # Convert WPICC and CAC to whole numbers
  WPICC100 <- WPICC*100
  CAC100 <- CAC*100
  config <- paste0(
    '_S', clust_per_seq,
    '_T', periods,
    '_m', subjects,
    '_WPICC', WPICC100,
    '_CAC', CAC100,
    '_theta', theta,
    sprintf('_seed%04d', seedval)
  )
  
  # Save MCMC results
  save(
    list=c('draws', 'parvals', 'diagnostics'),
    file=paste0(
      'results/',
      'MCMC',
      config,
      '.Rda'
    )
  )
  
  # Save full results for MCMC fits that fail diagnostics checks for review
  if (diagnostics$convergence_check=='FAIL' | diagnostics$div_trans_check=='FAIL'){
    dir.create('review')
    save(
      list=c('fit', 'parvals', 'diagnostics'),
      file=paste0(
        'review/',
        'MCMC',
        config,
        '.Rda'
      )
    )
  }

  # Save REML results
  save(
    list=c('sumreml', 'parvals', 'adj_SE', 'adj_ddf'),
    file=paste0(
      'results/',
      'REML',
      config,
      '.Rda'
    )
  )
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('Save results took:', format(time.taken, digits=4)))
}
