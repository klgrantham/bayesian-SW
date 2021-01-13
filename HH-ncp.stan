// Stan program for non-centered SW-HH-corr-time model
data {
  int<lower=0> clusters;
  int<lower=0> subjects;
  int<lower=0> periods;
  int<lower=0> N;
  int<lower=0> Q;
  matrix[clusters, periods] X;
  vector[N] Y;
}
parameters {
  // random variables in our model
  real theta;                       // treatment effect
  vector[clusters] C_tilde;         // vector of cluster effects
  real<lower=0> sig_sq_subject;     // between-subject variance
  real<lower=0,upper=1> WPICC;      // within-period intracluster correlation
  vector[periods] per;              // vector of period effects
}
transformed parameters {
  // we also want draws for our implied variance components
  real<lower=0> sig_sq_cluster = sig_sq_subject * WPICC / (1 - WPICC);
  vector[Q] mu;
  for (i in 1:clusters){
    for (j in 1:periods){
      int q = (i - 1) * periods + j;
      mu[q] = (per[j] + X[i,j] * theta) + sqrt(sig_sq_cluster) * C_tilde[i];
    }
  }
}
model {
  // this section computes the log joint likelihood conditional on the
  // parameters, in two parts
  // part 1: log prior density contribution
  theta ~ normal(0, 1e2);
  per ~ normal(0, 1e2); // normal(0, 1)
  sig_sq_subject ~ cauchy(0, 1) T[0,];
  WPICC ~ beta(1.5, 10.5); // peaked around 0.05, some mass around 0, divergent transitions occasionally
  C_tilde ~ normal(0, 1);
  // part 2: log likelihood contribution
  for (i in 1:clusters) {
    for (j in 1:periods) {
      int q = (i - 1) * periods + j; 
      for (k in 1:subjects) {
        int l;  // index into the data
        l = (i - 1) * periods * subjects + (j - 1) * subjects + k;
        Y[l] ~ normal(mu[q], sqrt(sig_sq_subject));
      }
    }
  }
}
