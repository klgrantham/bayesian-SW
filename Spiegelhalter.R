# install.packages("rstan", repos = "https://cloud.r-project.org/", dependencies = TRUE)

# ---- simulate data ----

# fixed params

clusters <- 10          # number of clusters
subjects <- 50          # subjects per cluster-period
periods <- 4            # number of periods

sig_sq_cluster <- 0.5   # between-cluster variance
sig_sq_subject <- 1.0   # between-subject variance

theta <- 1.32           # treatment effect

X <- matrix(data=0, nrow=clusters, ncol=periods)
X[1:clusters/2,] <- 1  # simple parallel design

rho <- sig_sq_cluster / (sig_sq_cluster + sig_sq_subject)

# random variables

set.seed(123)
C <- rnorm(n=clusters, mean=0, sd=sqrt(sig_sq_cluster))

Y <- rep(NA, clusters * periods * subjects)
e <- rnorm(n=clusters * periods * subjects, mean = 0,
           sd = sqrt(sig_sq_subject))
for (i in 1:clusters) {
  for (j in 1:periods) {
    for (k in 1:subjects) {
      # l is the index into the data
      l <- (i - 1) * periods * subjects + (j - 1) * subjects + k
      Y[l] <- X[i,j] * theta + C[i] + e[l]
    }
  }
}

# ---- fit model ----

library(rstan)
options(mc.cores = parallel::detectCores())

data <- list(
  clusters = clusters,
  subjects = subjects,
  periods = periods,
  N = clusters*subjects*periods,
  X = X,
  Y = Y
)
fit <- stan(file='Spiegelhalter.stan', data = data,
            chain = 4, iter = 11e3, warmup = 1e3, thin=2)

# ---- get diagnostics ----

pars <- c('theta', 'sig_sq_cluster', 'sig_sq_subject', 'rho', 'C')

# Check whether chains appear stationary
stan_trace(fit, pars=pars) + ggtitle('Trace plots of MCMC draws', 'Spiegelhalter model')

# Check that there is low autocorrelation for each parameter
stan_ac(fit, pars=pars) + ggtitle('Autocorrelograms of posterior draws', 'Spiegelhalter model')

# Inference summary table
print(fit, pars=pars)

## Posterior marginals
# Compare posterior marginals to true values
plot(fit, pars=pars) + ggtitle('Posterior marginal credible intervals', 'Spiegelhalter model')

stan_dens(fit, pars=pars) + ggtitle('Marginal posterior density estimates')

suppressMessages(library(bayesplot))
draws <- as.data.frame(fit, pars=pars)
true.vals <- c(theta, sig_sq_cluster, sig_sq_subject, rho, C)
mcmc_recover_hist(draws, true.vals) + ggtitle('Marginal posterior histograms vs true values')
