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
make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/bias_theta.jpg", p, width=9, height=7, units="in", dpi=800)


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
p <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/bias_WPICC.jpg", p, width=9, height=7, units="in", dpi=800)


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
p <- make_1x2_multiplot(p1, p2, mylegend, title)
ggsave("plots/bias_BPICC.jpg", p, width=9, height=7, units="in", dpi=800)


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
p <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/MSE_theta.jpg", p, width=9, height=7, units="in", dpi=800)


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
p <- make_2x2_multiplot(p1, p2, p3, p4, mylegend, title)
ggsave("plots/MSE_WPICC.jpg", p, width=9, height=7, units="in", dpi=800)

