# Plot performance measures across all combinations of
# trial configurations and true parameter values
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

source('results_tables.R')

plotminigrid <- function(df, ylims, title) {
  p <- ggplot(df, aes(x=S, y=value, group=method)) +
    geom_line(aes(color=method)) +
    geom_point(aes(color=method), show.legend=FALSE) +
    facet_grid(
      m ~ Tp,
      labeller = labeller(Tp = Tp.labs, m = m.labs)
    ) +
    expand_limits(y=ylims) +
    xlab("Clusters per sequence") +
    ylab("") +
    labs(title=title) +
    theme_bw()  +
    theme(
      strip.background = element_rect(
        color="white", fill="white", linetype="solid"
      ),
      plot.title=element_text(hjust=0.5, size=12),
      axis.title=element_text(size=10), axis.text=element_text(size=10),
      legend.key.width = unit(1.5, "cm"),
      legend.title=element_text(size=12), legend.text=element_text(size=12),
      legend.position="bottom"
    )
  return(p)
}

# Extract legend
g_legend <- function(a.gplot) {
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Multiplot framework
make_2x2_multiplot <- function(p1, p2, p3, p4, legend, title){
  p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2 + theme(legend.position="none"),
                                p3 + theme(legend.position="none"),
                                p4 + theme(legend.position="none"),
                                ncol=2),
                    legend, nrow=2, heights=c(10,1),
                    top=textGrob(title,
                                 gp=gpar(fontsize=18)))
  return(p)
}

make_1x2_multiplot <- function(p1, p2, legend, title){
  p <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                                p2 + theme(legend.position="none"),
                                ncol=2),
                    legend, nrow=2, heights=c(10,1),
                    top=textGrob(title,
                                 gp=gpar(fontsize=18)))
  return(p)
}

# Generate grid plot labels
Tp.labs <- c("T = 5", "T = 9")
names(Tp.labs) <- c("5", "9")

m.labs <- c("m = 10", "m = 100")
names(m.labs) <- c("10", "100")


## Plot bias for theta

# Determine common y-axis limits
minbiasHG <- signif(
  biasdfHG %>%
    filter(parameters=='theta') %>%
    summarise(min(value)),
  2)
minbiasHH <- signif(
  biasdfHG %>%
    filter(parameters=='theta') %>%
    summarise(min(value)),
  2)
maxbiasHG <- signif(
  biasdfHG %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
maxbiasHH <- signif(
  biasdfHH %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
ylims <- c(min(minbiasHG, minbiasHH), max(maxbiasHG, maxbiasHH))

# Get separate results blocks for plotting
biasdf_WPICC5_CAC80_theta <- biasdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta')
biasdf_WPICC10_CAC80_theta <- biasdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta')
biasdf_WPICC5_CAC100_theta <- biasdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta')
biasdf_WPICC10_CAC100_theta <- biasdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta')

# Generate grid plots
p1 <- plotminigrid(
  biasdf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dashed')
p2 <- plotminigrid(
  biasdf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dashed')
p3 <- plotminigrid(
  biasdf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dashed')
p4 <- plotminigrid(
  biasdf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dashed')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Bias of ", hat(theta)))
p_bias_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/bias_theta.jpg", p_bias_theta, width=9, height=7, units="in", dpi=800)


## Plot bias for WPICC (rho1)

# Determine common y-axis limits
minbiasHG <- signif(
  biasdfHG %>%
    filter(parameters=='WPICC') %>%
    summarise(min(value)),
  2)
minbiasHH <- signif(
  biasdfHG %>%
    filter(parameters=='WPICC') %>%
    summarise(min(value)),
  2)
maxbiasHG <- signif(
  biasdfHG %>%
    filter(parameters=='WPICC') %>%
    summarise(max(value)),
  2)
maxbiasHH <- signif(
  biasdfHH %>%
    filter(parameters=='WPICC') %>%
    summarise(max(value)),
  2)
ylims <- c(min(minbiasHG, minbiasHH), max(maxbiasHG, maxbiasHH))

# Get separate results blocks for plotting
biasdf_WPICC5_CAC80_WPICC <- biasdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='WPICC')
biasdf_WPICC10_CAC80_WPICC <- biasdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='WPICC')
biasdf_WPICC5_CAC100_WPICC <- biasdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='WPICC')
biasdf_WPICC10_CAC100_WPICC <- biasdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='WPICC')

# Generate grid plots
p1 <- plotminigrid(
  biasdf_WPICC5_CAC80_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dashed')
p2 <- plotminigrid(
  biasdf_WPICC10_CAC80_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dashed')
p3 <- plotminigrid(
  biasdf_WPICC5_CAC100_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dashed')
p4 <- plotminigrid(
  biasdf_WPICC10_CAC100_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dashed')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Bias of ", hat(rho[1])))
p_bias_WPICC <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/bias_WPICC.jpg", p_bias_WPICC, width=9, height=7, units="in", dpi=800)


## Plot bias for BPICC (rho2)
# Note: Only for HG model, not HH

# Determine common y-axis limits
minbiasHG <- signif(
  biasdfHG %>%
    filter(parameters=='BPICC') %>%
    summarise(min(value)),
  2)
maxbiasHG <- signif(
  biasdfHG %>%
    filter(parameters=='BPICC') %>%
    summarise(max(value)),
  2)
ylims <- c(min(minbiasHG), max(maxbiasHG))

# Get separate results blocks for plotting
biasdf_WPICC5_CAC80_BPICC <- biasdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='BPICC')
biasdf_WPICC10_CAC80_BPICC <- biasdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='BPICC')

# Generate grid plots
p1 <- plotminigrid(
  biasdf_WPICC5_CAC80_BPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dashed')
p2 <- plotminigrid(
  biasdf_WPICC10_CAC80_BPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dashed')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Bias of ", hat(rho[2])))
p_bias_BPICC <- make_1x2_multiplot(p1, p2, mylegend, title)
ggsave("plots/bias_BPICC.jpg", p_bias_BPICC, width=9, height=7, units="in", dpi=800)


## Plot MSE for theta

# Determine common y-axis limits
maxMSEHG <- signif(
  MSEdfHG %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
maxMSEHH <- signif(
  MSEdfHH %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
ylims <- c(0, max(maxMSEHG, maxMSEHH))

# Get separate results blocks for plotting
MSEdf_WPICC5_CAC80_theta <- MSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta')
MSEdf_WPICC10_CAC80_theta <- MSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta')
MSEdf_WPICC5_CAC100_theta <- MSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta')
MSEdf_WPICC10_CAC100_theta <- MSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta')

# Generate grid plots
p1 <- plotminigrid(
  MSEdf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
)
p2 <- plotminigrid(
  MSEdf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
)
p3 <- plotminigrid(
  MSEdf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
)
p4 <- plotminigrid(
  MSEdf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
)

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("MSE of ", hat(theta)))
p_MSE_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/MSE_theta.jpg", p_MSE_theta, width=9, height=7, units="in", dpi=800)


## Plot MSE for WPICC

# Determine common y-axis limits
maxMSEHG <- signif(
  MSEdfHG %>%
    filter(parameters=='WPICC') %>%
    summarise(max(value)),
  2)
maxMSEHH <- signif(
  MSEdfHH %>%
    filter(parameters=='WPICC') %>%
    summarise(max(value)),
  2)
ylims <- c(0, max(maxMSEHG, maxMSEHH))

# Get separate results blocks for plotting
MSEdf_WPICC5_CAC80_WPICC <- MSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='WPICC')
MSEdf_WPICC10_CAC80_WPICC <- MSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='WPICC')
MSEdf_WPICC5_CAC100_WPICC <- MSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='WPICC')
MSEdf_WPICC10_CAC100_WPICC <- MSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='WPICC')

# Generate grid plots
p1 <- plotminigrid(
  MSEdf_WPICC5_CAC80_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
)
p2 <- plotminigrid(
  MSEdf_WPICC10_CAC80_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
)
p3 <- plotminigrid(
  MSEdf_WPICC5_CAC100_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
)
p4 <- plotminigrid(
  MSEdf_WPICC10_CAC100_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
)

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("MSE of ", hat(rho[1])))
p_MSE_WPICC <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/MSE_WPICC.jpg", p_MSE_WPICC, width=9, height=7, units="in", dpi=800)


## Plot MSE for BPICC

# Determine common y-axis limits
maxMSEHG <- signif(
  MSEdfHG %>%
    filter(parameters=='BPICC') %>%
    summarise(max(value)),
  2)
maxMSEHH <- signif(
  MSEdfHH %>%
    filter(parameters=='BPICC') %>%
    summarise(max(value)),
  2)
ylims <- c(0, max(maxMSEHG, maxMSEHH))

# Get separate results blocks for plotting
MSEdf_WPICC5_CAC80_BPICC <- MSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='BPICC')
MSEdf_WPICC10_CAC80_BPICC <- MSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='BPICC')
MSEdf_WPICC5_CAC100_BPICC <- MSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='BPICC')
MSEdf_WPICC10_CAC100_BPICC <- MSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='BPICC')

# Generate grid plots
p1 <- plotminigrid(
  MSEdf_WPICC5_CAC80_BPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
)
p2 <- plotminigrid(
  MSEdf_WPICC10_CAC80_BPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
)

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("MSE of ", hat(rho[2])))
p_MSE_BPICC <- make_1x2_multiplot(p1, p2, mylegend, title)
ggsave("plots/MSE_BPICC.jpg", p_MSE_BPICC, width=9, height=7, units="in", dpi=800)


## Plot 95% confidence/credible interval coverage for theta

# Determine common y-axis limits
mincovHG <- signif(
  covdfHG %>%
    filter(parameters=='theta') %>%
    summarise(min(value)),
  2)
mincovHH <- signif(
  covdfHH %>%
    filter(parameters=='theta') %>%
    summarise(min(value)),
  2)
ylims <- c(min(mincovHG, mincovHH), 1.0)

# Get separate results blocks for plotting
covdf_WPICC5_CAC80_theta <- covdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta')
covdf_WPICC10_CAC80_theta <- covdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta')
covdf_WPICC5_CAC100_theta <- covdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta')
covdf_WPICC10_CAC100_theta <- covdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta')

# Generate grid plots
p1 <- plotminigrid(
  covdf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0.95, linetype='dashed')
p2 <- plotminigrid(
  covdf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0.95, linetype='dashed')
p3 <- plotminigrid(
  covdf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
) + geom_hline(yintercept=0.95, linetype='dashed')
p4 <- plotminigrid(
  covdf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
) + geom_hline(yintercept=0.95, linetype='dashed')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Confidence/credible interval coverage of ", hat(theta)))
p_cov_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/cov_theta.jpg", p_cov_theta, width=9, height=7, units="in", dpi=800)


## Plot empirical SE for theta

# Determine common y-axis limits
maxempSEHG <- signif(
  empSEdfHG %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
maxempSEHH <- signif(
  empSEdfHH %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
ylims <- c(0, max(maxempSEHG, maxempSEHH))

# Get separate results blocks for plotting
empSEdf_WPICC5_CAC80_theta <- empSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta')
empSEdf_WPICC10_CAC80_theta <- empSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta')
empSEdf_WPICC5_CAC100_theta <- empSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta')
empSEdf_WPICC10_CAC100_theta <- empSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta')

# Generate grid plots
p1 <- plotminigrid(
  empSEdf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
)
p2 <- plotminigrid(
  empSEdf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
)
p3 <- plotminigrid(
  empSEdf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
)
p4 <- plotminigrid(
  empSEdf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
)

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Empirical SE of ", hat(theta)))
p_empSE_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/empSE_theta.jpg", p_empSE_theta, width=9, height=7, units="in", dpi=800)


## Plot model-based SE for theta

# Determine common y-axis limits
maxmodSEHG <- signif(
  modSEdfHG %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
maxmodSEHH <- signif(
  modSEdfHH %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
ylims <- c(0, max(maxmodSEHG, maxmodSEHH))

# Get separate results blocks for plotting
empSEdf_WPICC5_CAC80_theta <- modSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta')
modSEdf_WPICC10_CAC80_theta <- modSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta')
modSEdf_WPICC5_CAC100_theta <- modSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta')
modSEdf_WPICC10_CAC100_theta <- modSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta')

# Generate grid plots
p1 <- plotminigrid(
  modSEdf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
)
p2 <- plotminigrid(
  modSEdf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
)
p3 <- plotminigrid(
  modSEdf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
)
p4 <- plotminigrid(
  modSEdf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
)

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Model-based SE of ", hat(theta)))
p_modSE_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/modSE_theta.jpg", p_modSE_theta, width=9, height=7, units="in", dpi=800)


## Plot empirical versus model-based SE for theta

# Determine common y-axis limits
maxmodSEHG <- signif(
  modSEdfHG %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
maxmodSEHH <- signif(
  modSEdfHH %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
maxempSEHG <- signif(
  empSEdfHG %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
maxempSEHH <- signif(
  empSEdfHH %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
ylims <- c(0, max(maxmodSEHG, maxmodSEHH, maxempSEHG, maxempSEHH))

# Combine empirical and model-based results
empSEdfHH$method <- as.factor(empSEdfHH$method)
levels(empSEdfHH$method) <- c('MCMC-empirical', 'REML-empirical')
empSEdfHG$method <- as.factor(empSEdfHG$method)
levels(empSEdfHG$method) <- c('MCMC-empirical', 'REML-empirical')
modSEdfHH$method <- as.factor(modSEdfHH$method)
levels(modSEdfHH$method) <- c('MCMC-model', 'REML-model')
modSEdfHG$method <- as.factor(modSEdfHG$method)
levels(modSEdfHG$method) <- c('MCMC-model', 'REML-model')
bothSEdfHH <- rbind(empSEdfHH, modSEdfHH)
bothSEdfHG <- rbind(empSEdfHG, modSEdfHG)

# Get separate results blocks for plotting
bothSEdf_WPICC5_CAC80_theta <- bothSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta')
bothSEdf_WPICC10_CAC80_theta <- bothSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta')
bothSEdf_WPICC5_CAC100_theta <- bothSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta')
bothSEdf_WPICC10_CAC100_theta <- bothSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta')

# Generate grid plots
p1 <- plotminigrid(
  bothSEdf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
)
p2 <- plotminigrid(
  bothSEdf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
)
p3 <- plotminigrid(
  bothSEdf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
)
p4 <- plotminigrid(
  bothSEdf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
)

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Empirical versus model-based SE of ", hat(theta)))
p_bothSE_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/bothSE_theta.jpg", p_bothSE_theta, width=9, height=7, units="in", dpi=800)
