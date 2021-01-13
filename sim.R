# Main procedure: data generation and model fitting
# with both MCMC and REML estimation
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

library(rstan)
library(lme4)

fit_models <- function(params, stanmod, seedval){
  # Inputs:
  #   params - a data frame containing input parameters for one particular
  #            simulation configuration
  # Outputs:
  #   MCMC_[config].Rda - saves an R data file containing the fitted MCMC model
  #                       (fit) and true parameter values (parvals)
  #   REML_[config].Rda - saves an R data file containing the fitted REML model
  #                       (remlfit) and true parameter values (parvals)
  
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
  
  parvals <- list(
    clust_per_seq = clust_per_seq,
    periods = periods,
    subjects = subjects,
    sig_sq_subject = sig_sq_subject,
    WPICC = WPICC,
    CAC = CAC,
    theta = theta,
    per = per,
    sig_sq_cluster = sig_sq_cluster,
    sig_sq_cp = sig_sq_cp
  )
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('Simulated data took:', format(time.taken, digits=4)))
  
  # ---- fit models ----
  
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
  
  fit <- sampling(stanmod, data = data, iter = 6e3, warmup = 1e3,
                  seed=seedval, control = list(adapt_delta = 0.9), cores=1)
  
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
    remlfit <- lmer(Y ~ treat + per + (1 | clust), data=dat, REML = TRUE)
    pars <- c('theta', 'WPICC', 'sig_sq_subject', 'sig_sq_cluster')
  } else {
    remlfit <- lmer(Y ~ treat + per + (1 | clust) + (1 | clustper), data=dat, REML = TRUE)
    pars <- c('theta', 'WPICC', 'CAC', 'sig_sq_subject', 'sig_sq_cluster', 'sig_sq_cp')
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('REML fit took:', format(time.taken, digits=4)))
  
  # Inference
  
  start.time <- Sys.time()
  
  print(fit, pars=pars)

  div <- 0
  chains <- dim(fit)[2] # Extract number of chains from fit
  for (ch in 1:chains) {
    divergent <- get_sampler_params(fit, inc_warmup=FALSE)[[ch]][,'divergent__']
    div <- div + sum(divergent)
  }
  print(paste('Divergent transitions:', div))
  
  print(summary(remlfit))
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  print(paste('Inference took:', format(time.taken, digits=4)))
  
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
    list=c('fit', 'parvals'),
    file=paste0(
      'results/',
      'MCMC',
      config,
      '.Rda'
    )
  )  

  # Save REML results
  save(
    list=c('remlfit', 'parvals'),
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
