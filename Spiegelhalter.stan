// Stan program for Spiegelhalter model
data {
  int<lower=0> clusters;
  int<lower=0> subjects;
  int<lower=0> periods;
  int<lower=0> N;
  matrix[clusters, periods] X;
  vector[N] Y;
}
parameters {
  // random variables in our model
  real theta;                    // treatment effect
  vector[clusters] C;            // vector of cluster effects
  real<lower=0> sig_sq_cluster;  // between-cluster variance
  real<lower=0> sig_sq_subject;  // between-subject variance
}
transformed parameters {
  // we also want draws for our implied correlation parameter
  real<lower=0> rho = sig_sq_cluster / (sig_sq_cluster + sig_sq_subject);
}
model {
  // this section computes the log joint likelihood conditional on the
  // parameters, in two parts
  // part 1: log prior density contribution
  theta ~ uniform(-1e4, 1e4);
  sig_sq_cluster ~ inv_gamma(1, 1);
  sig_sq_subject ~ inv_gamma(1, 1);
  // part 2: log likelihood contribution
  C ~ normal(0, sqrt(sig_sq_cluster));
  for (i in 1:clusters) {
    for (j in 1:periods) {
      for (k in 1:subjects) {
        int l;  // index into the data
        l = (i - 1) * periods * subjects + (j - 1) * subjects + k;
        Y[l] ~ normal(X[i,j] * theta + C[i], sqrt(sig_sq_subject));
      }
    }
  }
}
