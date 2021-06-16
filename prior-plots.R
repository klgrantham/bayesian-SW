load("~/src/bayesian-SW/results/results_S1_T5_m10_WPICC5_CAC80_theta0_seed0001.Rda")

library(ggplot2)
library(tidyverse)
library(extraDistr)

params <- c(
  "theta", "WPICC", "CAC", "sig_sq_subject", "sig_sq_cluster", "sig_sq_cp",
  paste0('C', 1:clusters)
)

plot_dens_true <- function(df, dfvals, dftrue, dftruevals, fillcolor){
  # Plots a density (prior or posterior) against the true parameter value
  # for all parameters
  
  orderedparams <- params
  df %>%
    mutate(
      parameters = ordered(parameters, levels=orderedparams)
    ) %>%
    ggplot(aes(x=dfvals)) +
    geom_density(aes(x=dfvals), fill=fillcolor, alpha=0.5) +
    geom_vline(aes(xintercept=dftruevals), color="steelblue", size=1, data=dftrue) +
    facet_wrap(~parameters, scales="free") +
    theme_bw()
}

sig2c_implied <- function(sig2e, WPICC, CAC){
  # Calculates implied cluster variance using draws from
  # subject variance and within-cluster correlation components
  
  return(CAC * sig2e * WPICC / (1 - WPICC))
}

sig2cp_implied <- function(sig2e, WPICC, sig2c){
  # Calculates implied between-cluster-period variance using
  # draws from subject and cluster variances and within-period ICC
  
  return(sig2e * WPICC / (1 - WPICC) - sig2c)
}

plot_prior <- function(df) {
  p <- ggplot(df) +
        geom_line(aes(x=x, y=y), color="darkorange", size=1) +
        ylab("") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size=20),
              axis.title=element_text(size=18),
              axis.text=element_text(size=18))
}

plot_dens <- function(densvals) {
  p <- ggplot() +
        geom_density(aes(x=densvals), color="slateblue", size=1) +
        ylab("") +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5, size=20),
              axis.title=element_text(size=18),
              axis.text=element_text(size=18))
}

# Specify prior distributions for plotting
theta_prior <- data.frame(x=seq(-400, 400, length.out=100))
theta_prior$y <- dnorm(theta_prior$x, 0, 100)
p_theta <- plot_prior(theta_prior) +
  geom_vline(aes(xintercept=0), color="steelblue", linetype="dashed", size=1) +
  xlab(expression(theta)) +
  ggtitle(expression(paste(theta, " ~ N(0, ", 10^4, ")")))
#ggsave("plots/prior_theta.jpg", width=6, height=6, units="in", dpi=600)

betaj_prior <- data.frame(x=seq(-400, 400, length.out=100))
betaj_prior$y <- dnorm(theta_prior$x, 0, 100)
p_betaj <- plot_prior(betaj_prior) +
  xlab(expression(beta[j])) +
  ggtitle(expression(paste(beta[j], " ~ N(0, ", 10^4, ")")))

WPICC_prior <- data.frame(x=seq(0, 1, length.out=100))
WPICC_prior$y <- dbeta(WPICC_prior$x, 1.5, 10.5)
p_WPICC <- plot_prior(WPICC_prior) +
  geom_vline(aes(xintercept=0.05), color="steelblue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=0.1), color="darkseagreen", linetype="dashed", size=1) +
  xlab(expression(rho[1])) +
  ggtitle(expression(paste(rho[1], " ~ Beta(1.5, 10.5)")))
#ggsave("plots/prior_WPICC.jpg", width=6, height=6, units="in", dpi=600)

CAC_prior <- data.frame(x=seq(0, 1, length.out=100))
CAC_prior$y <- dbeta(CAC_prior$x, 5, 2)
p_CAC <- plot_prior(CAC_prior) +
  geom_vline(aes(xintercept=0.8), color="steelblue", linetype="dashed", size=1) +
  xlab(expression(r)) +
  ggtitle(expression(paste(r, " ~ Beta(5, 2)")))
#ggsave("plots/prior_CAC.jpg", width=6, height=6, units="in", dpi=600)

sig2e_prior <- data.frame(x=seq(0, 50, length.out=100))
sig2e_prior$y <- dhcauchy(sig2e_prior$x, 1)
p_sig2e <- plot_prior(sig2e_prior) +
  geom_vline(aes(xintercept=1.0), color="steelblue", linetype="dashed", size=1) +
  xlab(expression(sigma[e]^2)) +
  ggtitle(expression(paste(sigma[e]^2, " ~ Half-Cauchy(0, 1)")))
#ggsave("plots/prior_sig2e.jpg", width=6, height=6, units="in", dpi=600)

p_all <- grid.arrange(arrangeGrob(p_theta, p_betaj, p_sig2e, p_WPICC, p_CAC, nrow=2))
ggsave("plots/allpriors.jpg", p_all, width=12, height=7, units="in", dpi=400)

# Implied distributions
sig2e_vals <- rhcauchy(1e5, 1)
WPICC_vals <- rbeta(1e5, 1.5, 10.5)
CAC_vals <- rbeta(1e5, 5, 2)
sig2c_prior <- sig2c_implied(sig2e_vals, WPICC_vals, CAC_vals)
sig2cp_prior <- sig2cp_implied(sig2e_vals, WPICC_vals, sig2c_prior)
BPICC_prior <- WPICC_vals * CAC_vals

# implied between-cluster variance
sig_sq_c1 <- sig2c_implied(sig2e=1, WPICC=0.05, CAC=0.8)
sig_sq_c2 <- sig2c_implied(sig2e=1, WPICC=0.1, CAC=0.8)
sig_sq_c3 <- sig2c_implied(sig2e=1, WPICC=0.05, CAC=1.0)
sig_sq_c4 <- sig2c_implied(sig2e=1, WPICC=0.1, CAC=1.0)

# implied between-cluster-period variance
sig_sq_cp1 <- sig2cp_implied(sig2e=1, WPICC=0.05, sig2c=sig_sq_c1)
sig_sq_cp2 <- sig2cp_implied(sig2e=1, WPICC=0.1, sig2c=sig_sq_c2)

# implied between-period ICC
BPICC1 <- 0.05 * 0.8
BPICC2 <- 0.1 * 0.8

plot_dens(BPICC_prior) +
  geom_vline(aes(xintercept=BPICC1), color="steelblue", size=1) +
  geom_vline(aes(xintercept=BPICC2), color="darkseagreen", size=1) +
  xlab(expression(rho[2])) +
  ggtitle(expression(paste("Implied distribution of ", rho[2]))) +
  xlim(c(0, 1))
ggsave("plots/priors_BPICC_implied.jpg", width=6, height=6, units="in", dpi=600)

p_sig2c <- plot_dens(sig2c_prior) +
  xlab(expression(sigma[C]^2)) +
  ggtitle(expression(paste("Implied distribution of ", sigma[C]^2))) +
  xlim(c(0, 5))
ggsave("plots/priors_sig2c_implied.jpg", width=6, height=6, units="in", dpi=600)

p_sig2cp <- plot_dens(sig2cp_prior) +
  xlab(expression(sigma[CP]^2)) +
  ggtitle(expression(paste("Implied distribution of ", sigma[CP]^2))) +
  xlim(c(0, 5))
ggsave("plots/priors_sig2cp_implied.jpg", width=6, height=6, units="in", dpi=600)

p_vars <- grid.arrange(arrangeGrob(p_sig2c, p_sig2cp, nrow=1))
ggsave("plots/impliedpriors.jpg", p_vars, width=8, height=4, units="in", dpi=400)
