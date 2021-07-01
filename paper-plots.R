# Plot performance measures across all combinations of
# trial configurations and true parameter values
#
# Kelsey Grantham (kelsey.grantham@monash.edu)
#

library(ggplot2)

source('process-results.R')

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


## Figure 3: Bias for theta

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
ggsave("plots/figure3_bias_theta.jpg", p_bias_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/figure3_bias_theta.pdf", p_bias_theta, width=9, height=7, units="in", dpi=600)


## Figure 4: MSE for theta

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
ggsave("plots/figure4_MSE_theta.jpg", p_MSE_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/figure4_MSE_theta.pdf", p_MSE_theta, width=9, height=7, units="in", dpi=600)


## Figure 5: 95% confidence/credible interval coverage for theta

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
ggsave("plots/figure5_cov_theta.jpg", p_cov_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/figure5_cov_theta.pdf", p_cov_theta, width=9, height=7, units="in", dpi=600)


## Figure 6: Relative % error in model-based SE for theta

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
ggsave("plots/figure6_pcterrmodSE_theta.jpg", p_pcterrmodSE_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/figure6_pcterrmodSE_theta.pdf", p_pcterrmodSE_theta, width=9, height=7, units="in", dpi=600)


## Figure 7: Average 95% confidence/credible interval length for theta

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
maxintlenHG <- signif(
  intlendfHG %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
maxintlenHH <- signif(
  intlendfHH %>%
    filter(parameters=='theta') %>%
    summarise(max(value)),
  2)
ylims <- c(0, max(maxintlenHG, maxintlenHH))

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
title <- expression(paste("Average confidence/credible interval width for the treatment effect, ", hat(theta)))
p_intlen_theta <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/figure7_intlen_theta.jpg", p_intlen_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/figure7_intlen_theta.pdf", p_intlen_theta, width=9, height=7, units="in", dpi=600)


## Figure 8: Bias for WPICC (rho1)

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
ggsave("plots/figure8_bias_WPICC.jpg", p_bias_WPICC, width=9, height=7, units="in", dpi=400)
ggsave("plots/figure8_bias_WPICC.pdf", p_bias_WPICC, width=9, height=7, units="in", dpi=600)


## Figure 9: MSE for WPICC

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
ggsave("plots/figure9_MSE_WPICC.jpg", p_MSE_WPICC, width=9, height=7, units="in", dpi=400)
ggsave("plots/figure9_MSE_WPICC.pdf", p_MSE_WPICC, width=9, height=7, units="in", dpi=600)


## Figure 10: Bias for CAC (r)
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
ggsave("plots/figure10_bias_CAC.jpg", p_bias_CAC, width=9, height=4, units="in", dpi=400)
ggsave("plots/figure10_bias_CAC.pdf", p_bias_CAC, width=9, height=4, units="in", dpi=600)


## Figure 11: MSE for CAC

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
ggsave("plots/figure11_MSE_CAC.jpg", p_MSE_CAC, width=9, height=4, units="in", dpi=400)
ggsave("plots/figure11_MSE_CAC.pdf", p_MSE_CAC, width=9, height=4, units="in", dpi=600)


# Figures D1 and D2: Histogram plots of MCMC posterior medians, comparing replicates with divergent
# transitions to replicates without divergent transitions
load(file='estimates/estimates_S1_T5_m100_WPICC5_CAC100_theta0.Rda')

# Extract theta and WPICC posterior medians, divergent and non-divergent
MCMC_medians_div <- est$MCMC_medians %>%
  filter(div==1) %>%
  select(-div) %>%
  select(c(theta, WPICC)) %>%
  add_column(Replicates_type="Divergent")
MCMC_medians_nodiv <- est$MCMC_medians %>%
  filter(div==0) %>%
  select(-div) %>%
  select(c(theta, WPICC)) %>%
  add_column(Replicates_type="Non-divergent")

MCMC_medians <- rbind(MCMC_medians_nodiv, MCMC_medians_div)

# theta posterior medians
p_theta <- ggplot(MCMC_medians, aes(x=theta, fill=Replicates_type, colour=Replicates_type)) +
  geom_histogram(position="identity", alpha=0.5) +
  geom_vline(aes(xintercept=0), colour="black", linetype="dashed") +
  xlab(expression(hat(theta))) +
  ylab("count") +
  labs(title=expression(paste("Posterior medians for the treatment effect, ", hat(theta))),
       subtitle=expression(paste("S=1, T=5, m=100, ", rho[1], "=0.05, r=1.0")),
       fill="Replicate type", colour="Replicate type") +
  theme_bw()  +
  theme(
    plot.title=element_text(hjust=0.5, size=18),
    plot.subtitle=element_text(hjust=0.5, size=16),
    axis.title=element_text(size=16), axis.text=element_text(size=14),
    legend.key.width = unit(1.5, "cm"),
    legend.title=element_text(size=12), legend.text=element_text(size=12),
    legend.position="bottom"
  )
ggsave("plots/figureD1_theta_medians_div_nodiv.jpg", p_theta, width=9, height=7, units="in", dpi=400)
ggsave("plots/figureD1_theta_medians_div_nodiv.pdf", p_theta, width=9, height=7, units="in", dpi=600)

# WPICC posterior medians
p_WPICC <- ggplot(MCMC_medians, aes(x=WPICC, fill=Replicates_type, colour=Replicates_type)) +
  geom_histogram(position="identity", alpha=0.5) +
  geom_vline(aes(xintercept=0.05), colour="black", linetype="dashed") +
  xlab(expression(hat(rho[1]))) +
  ylab("count") +
  labs(title=expression(paste("Posterior medians for the within-period intracluster correlation, ", hat(rho[1]))),
       subtitle=expression(paste("S=1, T=5, m=100, ", rho[1], "=0.05, r=1.0")),
       fill="Replicate type", colour="Replicate type") +
  theme_bw()  +
  theme(
    plot.title=element_text(hjust=0.5, size=18),
    plot.subtitle=element_text(hjust=0.5, size=16),
    axis.title=element_text(size=16), axis.text=element_text(size=14),
    legend.key.width = unit(1.5, "cm"),
    legend.title=element_text(size=12), legend.text=element_text(size=12),
    legend.position="bottom"
  )
ggsave("plots/figureD2_WPICC_medians_div_nodiv.jpg", p_WPICC, width=9, height=7, units="in", dpi=400)
ggsave("plots/figureD2_WPICC_medians_div_nodiv.pdf", p_WPICC, width=9, height=7, units="in", dpi=600)
