// Stan program for Spiegelhalter model
data {
  int cluster_periods;
  int clusters;
  int periods;
  vector[cluster_periods] subjects;
  vector[cluster_periods] invSubjects;
  matrix[cluster_periods, clusters] Cselect;
  vector[cluster_periods] means;
  vector[cluster_periods] X;
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
  means ~ normal(X * theta + Cselect * C, sqrt(sig_sq_subject * invSubjects));
}
