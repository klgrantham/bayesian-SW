// Stan program for Sufficient-SW-HG-corr model
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
  vector[cluster_periods] CP;       // vector of cluster-period effects
  real<lower=0> sig_sq_subject;  // between-subject variance
  real<lower=0,upper=1> WPICC;      // within-period intracluster correlation
  real<lower=0,upper=1> CAC;        // cluster autocorrelation
}
transformed parameters {
  // we also want draws for our implied variance components
  real<lower=0> sig_sq_cluster = CAC * sig_sq_subject * WPICC / (1 - WPICC);
  real<lower=0> sig_sq_cp = sig_sq_subject * WPICC / (1 - WPICC) - sig_sq_cluster;
}
model {
  // this section computes the log joint likelihood conditional on the
  // parameters, in two parts
  // part 1: log prior density contribution
  theta ~ uniform(-1e4, 1e4);
  sig_sq_subject ~ inv_gamma(2, 2);
  WPICC ~ beta(1.2, 5);
  CAC ~ beta(5, 2);
  // part 2: log likelihood contribution
  C ~ normal(0, sqrt(sig_sq_cluster));
  CP ~ normal(0, sqrt(sig_sq_cp));
  means ~ normal(X * theta + Cselect * C + CP, sqrt(sig_sq_subject * invSubjects));
}
