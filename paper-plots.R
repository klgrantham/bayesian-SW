# Plot performance measures across all combinations of
# trial configurations and true parameter values
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

source('results_tables.R')

plotminigrid <- function(df, ylims, title) {
  p <- ggplot(df, aes(x=Sfac, y=value, group=Method,
                      color=Method, linetype=Method)) +
    geom_line(size=1) +
    geom_point(size=2) +
    facet_grid(
      m ~ Tp,
      labeller = labeller(Tp = Tp.labs, m = m.labs)
    ) +
    expand_limits(y=ylims) +
    xlab("Clusters per sequence") +
    ylab("") +
    labs(title=title) +
    scale_x_discrete(breaks=c("1", "2", "5")) +
    theme_bw() +
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
  filter(rho1==0.05 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
biasdf_WPICC10_CAC80_theta <- biasdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
biasdf_WPICC5_CAC100_theta <- biasdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
biasdf_WPICC10_CAC100_theta <- biasdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))

# Generate grid plots
p1 <- plotminigrid(
  biasdf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')
p2 <- plotminigrid(
  biasdf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')
p3 <- plotminigrid(
  biasdf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dotted')
p4 <- plotminigrid(
  biasdf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dotted')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Bias for the treatment effect, ", hat(theta)))
p_bias_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/bias_theta.jpg", p_bias_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/bias_theta.pdf", p_bias_theta, width=9, height=7, units="in", dpi=600)


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
  filter(rho1==0.05 & r==0.8 & parameters=='WPICC') %>%
  mutate(Sfac=as.factor(S))
biasdf_WPICC10_CAC80_WPICC <- biasdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='WPICC') %>%
  mutate(Sfac=as.factor(S))
biasdf_WPICC5_CAC100_WPICC <- biasdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='WPICC') %>%
  mutate(Sfac=as.factor(S))
biasdf_WPICC10_CAC100_WPICC <- biasdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='WPICC') %>%
  mutate(Sfac=as.factor(S))

# Generate grid plots
p1 <- plotminigrid(
  biasdf_WPICC5_CAC80_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')
p2 <- plotminigrid(
  biasdf_WPICC10_CAC80_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')
p3 <- plotminigrid(
  biasdf_WPICC5_CAC100_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dotted')
p4 <- plotminigrid(
  biasdf_WPICC10_CAC100_WPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dotted')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Bias for the within-period intracluster correlation, ", hat(rho[1])))
p_bias_WPICC <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/bias_WPICC.jpg", p_bias_WPICC, width=9, height=7, units="in", dpi=400)
ggsave("plots/bias_WPICC.pdf", p_bias_WPICC, width=9, height=7, units="in", dpi=600)


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
  filter(rho1==0.05 & r==0.8 & parameters=='BPICC') %>%
  mutate(Sfac=as.factor(S))
biasdf_WPICC10_CAC80_BPICC <- biasdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='BPICC') %>%
  mutate(Sfac=as.factor(S))

# Generate grid plots
p1 <- plotminigrid(
  biasdf_WPICC5_CAC80_BPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')
p2 <- plotminigrid(
  biasdf_WPICC10_CAC80_BPICC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Bias for the between-period intracluter correlation, ", hat(rho[2])))
p_bias_BPICC <- make_1x2_multiplot(p1, p2, mylegend, title)
ggsave("plots/bias_BPICC.jpg", p_bias_BPICC, width=9, height=4, units="in", dpi=400)
ggsave("plots/bias_BPICC.pdf", p_bias_BPICC, width=9, height=4, units="in", dpi=600)


## Plot bias for CAC (r)
# Note: Only for HG model, not HH

# Determine common y-axis limits
minbiasHG <- signif(
  biasdfHG %>%
    filter(parameters=='CAC') %>%
    summarise(min(value)),
  2)
maxbiasHG <- signif(
  biasdfHG %>%
    filter(parameters=='CAC') %>%
    summarise(max(value)),
  2)
ylims <- c(min(minbiasHG), max(maxbiasHG))

# Get separate results blocks for plotting
biasdf_WPICC5_CAC80_CAC <- biasdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='CAC') %>%
  mutate(Sfac=as.factor(S))
biasdf_WPICC10_CAC80_CAC <- biasdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='CAC') %>%
  mutate(Sfac=as.factor(S))

# Generate grid plots
p1 <- plotminigrid(
  biasdf_WPICC5_CAC80_CAC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')
p2 <- plotminigrid(
  biasdf_WPICC10_CAC80_CAC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Bias for the cluster autocorrelation, ", hat(r)))
p_bias_CAC <- make_1x2_multiplot(p1, p2, mylegend, title)
ggsave("plots/bias_CAC.jpg", p_bias_CAC, width=9, height=4, units="in", dpi=400)
ggsave("plots/bias_CAC.pdf", p_bias_CAC, width=9, height=4, units="in", dpi=600)


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
  filter(rho1==0.05 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
MSEdf_WPICC10_CAC80_theta <- MSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
MSEdf_WPICC5_CAC100_theta <- MSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
MSEdf_WPICC10_CAC100_theta <- MSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))

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
title <- expression(paste("MSE for the treatment effect, ", hat(theta)))
p_MSE_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/MSE_theta.jpg", p_MSE_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/MSE_theta.pdf", p_MSE_theta, width=9, height=7, units="in", dpi=600)


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
  filter(rho1==0.05 & r==0.8 & parameters=='WPICC') %>%
  mutate(Sfac=as.factor(S))
MSEdf_WPICC10_CAC80_WPICC <- MSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='WPICC') %>%
  mutate(Sfac=as.factor(S))
MSEdf_WPICC5_CAC100_WPICC <- MSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='WPICC') %>%
  mutate(Sfac=as.factor(S))
MSEdf_WPICC10_CAC100_WPICC <- MSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='WPICC') %>%
  mutate(Sfac=as.factor(S))

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
title <- expression(paste("MSE for the within-period intracluster correlation, ", hat(rho[1])))
p_MSE_WPICC <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/MSE_WPICC.jpg", p_MSE_WPICC, width=9, height=7, units="in", dpi=400)
ggsave("plots/MSE_WPICC.pdf", p_MSE_WPICC, width=9, height=7, units="in", dpi=600)


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
  filter(rho1==0.05 & r==0.8 & parameters=='BPICC') %>%
  mutate(Sfac=as.factor(S))
MSEdf_WPICC10_CAC80_BPICC <- MSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='BPICC') %>%
  mutate(Sfac=as.factor(S))

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
title <- expression(paste("MSE for the between-period intracluster correlation, ", hat(rho[2])))
p_MSE_BPICC <- make_1x2_multiplot(p1, p2, mylegend, title)
ggsave("plots/MSE_BPICC.jpg", p_MSE_BPICC, width=9, height=4, units="in", dpi=400)
ggsave("plots/MSE_BPICC.pdf", p_MSE_BPICC, width=9, height=4, units="in", dpi=600)


## Plot MSE for CAC

# Determine common y-axis limits
maxMSEHG <- signif(
  MSEdfHG %>%
    filter(parameters=='CAC') %>%
    summarise(max(value)),
  2)
maxMSEHH <- signif(
  MSEdfHH %>%
    filter(parameters=='CAC') %>%
    summarise(max(value)),
  2)
ylims <- c(0, max(maxMSEHG, maxMSEHH))

# Get separate results blocks for plotting
MSEdf_WPICC5_CAC80_CAC <- MSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='CAC') %>%
  mutate(Sfac=as.factor(S))
MSEdf_WPICC10_CAC80_CAC <- MSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='CAC') %>%
  mutate(Sfac=as.factor(S))

# Generate grid plots
p1 <- plotminigrid(
  MSEdf_WPICC5_CAC80_CAC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
)
p2 <- plotminigrid(
  MSEdf_WPICC10_CAC80_CAC,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
)

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("MSE for the cluster autocorrelation, ", hat(r)))
p_MSE_CAC <- make_1x2_multiplot(p1, p2, mylegend, title)
ggsave("plots/MSE_CAC.jpg", p_MSE_CAC, width=9, height=4, units="in", dpi=400)
ggsave("plots/MSE_CAC.pdf", p_MSE_CAC, width=9, height=4, units="in", dpi=600)


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
  filter(rho1==0.05 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
covdf_WPICC10_CAC80_theta <- covdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
covdf_WPICC5_CAC100_theta <- covdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
covdf_WPICC10_CAC100_theta <- covdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))

# Generate grid plots
p1 <- plotminigrid(
  covdf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0.95, linetype='dotted')
p2 <- plotminigrid(
  covdf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0.95, linetype='dotted')
p3 <- plotminigrid(
  covdf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
) + geom_hline(yintercept=0.95, linetype='dotted')
p4 <- plotminigrid(
  covdf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
) + geom_hline(yintercept=0.95, linetype='dotted')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Confidence/credible interval coverage for the treatment effect, ", hat(theta)))
p_cov_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/cov_theta.jpg", p_cov_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/cov_theta.pdf", p_cov_theta, width=9, height=7, units="in", dpi=600)


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
  filter(rho1==0.05 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
empSEdf_WPICC10_CAC80_theta <- empSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
empSEdf_WPICC5_CAC100_theta <- empSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
empSEdf_WPICC10_CAC100_theta <- empSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))

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
title <- expression(paste("Empirical SE for the treatment effect, ", hat(theta)))
p_empSE_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/empSE_theta.jpg", p_empSE_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/empSE_theta.pdf", p_empSE_theta, width=9, height=7, units="in", dpi=600)


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
modSEdf_WPICC5_CAC80_theta <- modSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
modSEdf_WPICC10_CAC80_theta <- modSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
modSEdf_WPICC5_CAC100_theta <- modSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
modSEdf_WPICC10_CAC100_theta <- modSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))

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
title <- expression(paste("Model-based SE for the treatment effect, ", hat(theta)))
p_modSE_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/modSE_theta.jpg", p_modSE_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/modSE_theta.pdf", p_modSE_theta, width=9, height=7, units="in", dpi=600)


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
empSEdfHH$Method <- as.factor(empSEdfHH$Method)
levels(empSEdfHH$Method) <- c('MCMC-empirical', 'REML-empirical')
empSEdfHG$Method <- as.factor(empSEdfHG$Method)
levels(empSEdfHG$Method) <- c('MCMC-empirical', 'REML-empirical')
modSEdfHH$Method <- as.factor(modSEdfHH$Method)
levels(modSEdfHH$Method) <- c('MCMC-model', 'REML-model', 'REML (KR)-model')
modSEdfHG$Method <- as.factor(modSEdfHG$Method)
levels(modSEdfHG$Method) <- c('MCMC-model', 'REML-model', 'REML (KR)-model')
bothSEdfHH <- rbind(empSEdfHH, modSEdfHH)
bothSEdfHG <- rbind(empSEdfHG, modSEdfHG)

# Get separate results blocks for plotting
bothSEdf_WPICC5_CAC80_theta <- bothSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
bothSEdf_WPICC10_CAC80_theta <- bothSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
bothSEdf_WPICC5_CAC100_theta <- bothSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
bothSEdf_WPICC10_CAC100_theta <- bothSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))

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
title <- expression(paste("Empirical and model-based SE for the treatment effect, ", hat(theta)))
p_bothSE_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/bothSE_theta.jpg", p_bothSE_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/bothSE_theta.pdf", p_bothSE_theta, width=9, height=7, units="in", dpi=600)


## Plot relative % error in model-based SE for theta

# Determine common y-axis limits
maxpcterrmodSEHG <- signif(
  pcterrmodSEdfHG %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
maxpcterrmodSEHH <- signif(
  pcterrmodSEdfHH %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
ylims <- c(0, max(maxpcterrmodSEHG, maxpcterrmodSEHH))

# Get separate results blocks for plotting
pcterrmodSEdf_WPICC5_CAC80_theta <- pcterrmodSEdfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
pcterrmodSEdf_WPICC10_CAC80_theta <- pcterrmodSEdfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
pcterrmodSEdf_WPICC5_CAC100_theta <- pcterrmodSEdfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
pcterrmodSEdf_WPICC10_CAC100_theta <- pcterrmodSEdfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))

# Generate grid plots
p1 <- plotminigrid(
  pcterrmodSEdf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')
p2 <- plotminigrid(
  pcterrmodSEdf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
) + geom_hline(yintercept=0, linetype='dotted')
p3 <- plotminigrid(
  pcterrmodSEdf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dotted')
p4 <- plotminigrid(
  pcterrmodSEdf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
) + geom_hline(yintercept=0, linetype='dotted')

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Relative % error in model-based SE for the treatment effect, ", hat(theta)))
p_pcterrmodSE_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/pcterrmodSE_theta.jpg", p_pcterrmodSE_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/pcterrmodSE_theta.pdf", p_pcterrmodSE_theta, width=9, height=7, units="in", dpi=600)


## Plot average 95% confidence/credible interval length for theta

# Determine common y-axis limits
minintlenHG <- signif(
  intlendfHG %>%
    filter(parameters=='theta') %>%
    summarise(min(value)),
  2)
minintlenHH <- signif(
  intlendfHH %>%
    filter(parameters=='theta') %>%
    summarise(min(value)),
  2)
ylims <- c(min(minintlenHG, minintlenHH), 1.0)

# Get separate results blocks for plotting
intlendf_WPICC5_CAC80_theta <- intlendfHG %>%
  filter(rho1==0.05 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
intlendf_WPICC10_CAC80_theta <- intlendfHG %>%
  filter(rho1==0.1 & r==0.8 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
intlendf_WPICC5_CAC100_theta <- intlendfHH %>%
  filter(rho1==0.05 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))
intlendf_WPICC10_CAC100_theta <- intlendfHH %>%
  filter(rho1==0.1 & r==1.0 & parameters=='theta') %>%
  mutate(Sfac=as.factor(S))

# Generate grid plots
p1 <- plotminigrid(
  intlendf_WPICC5_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 0.8"))
)
p2 <- plotminigrid(
  intlendf_WPICC10_CAC80_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 0.8"))
)
p3 <- plotminigrid(
  intlendf_WPICC5_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.05", ", r = 1.0"))
)
p4 <- plotminigrid(
  intlendf_WPICC10_CAC100_theta,
  ylims=ylims,
  title=expression(paste(rho[1], " = 0.1", ", r = 1.0"))
)

# Combine grid plots
mylegend <- g_legend(p1)
title <- expression(paste("Average confidence/credible interval length for the treatment effect, ", hat(theta)))
p_intlen_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/intlen_theta.jpg", p_intlen_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/intlen_theta.pdf", p_intlen_theta, width=9, height=7, units="in", dpi=600)
