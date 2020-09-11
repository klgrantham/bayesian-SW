# install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

# ---- simulate data ----

# fixed params

clusters <- 4          # number of clusters
subjects <- 100          # subjects per cluster-period
periods <- 5            # number of periods

sig_sq_subject <- 0.9   # between-subject variance
WPICC <- 0.5            # within-period intracluster correlation
CAC <- 0.5              # cluster autocorrelation

theta <- 2           # treatment effect

# stepped wedge design matrix
# Note: clusters must be a multiple of (periods-1)
X <- matrix(data=0, nrow=periods-1, ncol=periods)
for(i in 1:(periods-1)){
  X[i,(i+1):periods] <- 1
}
X <- X %x% rep(1, clusters/(periods-1))

# implied between-cluster variance
sig_sq_cluster <- CAC * sig_sq_subject * WPICC / (1 - WPICC)

# implied between-cluster-period variance
sig_sq_cp <- sig_sq_subject * WPICC / (1 - WPICC) - sig_sq_cluster

# random variables

C <- rnorm(n=clusters, mean=0, sd=sqrt(sig_sq_cluster))
CP <- rnorm(n=clusters * periods, mean=0, sd=sqrt(sig_sq_cp))
Y <- rep(NA, clusters * periods * subjects)
e <- rnorm(n=clusters * periods * subjects, mean = 0, sd = sqrt(sig_sq_subject))
for (i in 1:clusters) {
  for (j in 1:periods) {
    for (k in 1:subjects) {
      # l is the index into the data
      l <- (i - 1) * periods * subjects + (j - 1) * subjects + k
      Y[l] <- X[i,j] * theta + C[i] + CP[(i - 1) * periods + j] + e[l]
    }
  }
}

# ---- fit model ----

library(rstan)

data <- list(
  clusters = clusters,
  subjects = subjects,
  periods = periods,
  cluster_periods = clusters * periods,
  A = clusters * periods * subjects,
  X = as.vector(kronecker(as.vector(t(X)), matrix(rep(1, subjects), ncol=1))),
  Cselect = kronecker(diag(clusters), matrix(rep(1, periods * subjects), ncol=1)),
  CPselect = kronecker(diag(clusters * periods), matrix(rep(1, subjects), ncol=1)),
  Y = Y
)
fit <- stan(file='SW-HG-corr-vec.stan', data = data, iter = 6e3, warmup = 1e3, cores=1)

# ---- get diagnostics ----

pars <- c('theta', 'WPICC', 'CAC', 'sig_sq_subject', 'sig_sq_cluster', 'sig_sq_cp', 'C')

# Inference summary table
print(fit, pars=pars, digits_summary=4)

# Check whether chains appear stationary
stan_trace(fit, pars=pars) + ggtitle('Trace plots of MCMC draws')

# Check that there is low autocorrelation for each parameter
stan_ac(fit, pars=pars) + ggtitle('Autocorrelograms of posterior draws')

## Posterior marginals
# Compare posterior marginals to true values
plot(fit, pars=pars) + ggtitle('Posterior marginal credible intervals')

stan_dens(fit, pars=pars) + ggtitle('Marginal posterior density estimates')

suppressMessages(library(bayesplot))
draws <- as.data.frame(fit, pars=pars)
true.vals <- c(theta, WPICC, CAC, sig_sq_subject, sig_sq_cluster, sig_sq_cp, C)
mcmc_recover_hist(draws, true.vals) + ggtitle('Marginal posterior histograms vs true values')

save.image('mcmc_results_vec.Rda')
