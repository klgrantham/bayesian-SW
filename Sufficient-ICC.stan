// Stan program for Sufficient-ICC model
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
  real theta;                       // treatment effect
  vector[clusters] C;               // vector of cluster effects
  real<lower=0,upper=1> rho;        // intra-cluster correlation
  real<lower=0.01> sig_sq_subject;  // between-subject variance
}
transformed parameters {
  // we also want draws for our implied between-cluster variance
  real<lower=0> sig_sq_cluster = (rho/(1-rho))*sig_sq_subject;
}
model {
  // this section computes the log joint likelihood conditional on the
  // parameters, in two parts
  // part 1: log prior density contribution
  theta ~ uniform(-1e4, 1e4);
  log(sig_sq_subject) ~ uniform(-5, 5);
  target += -log(sig_sq_subject);
  rho ~ uniform(0, 1);
  //rho ~ beta(2, 38);
  // part 2: log likelihood contribution
  C ~ normal(0, sqrt(sig_sq_cluster));
  means ~ normal(X * theta + Cselect * C, sqrt(sig_sq_subject * invSubjects));
}
