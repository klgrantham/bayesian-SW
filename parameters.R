# Generate .csv file of parameter configurations for simulation study
#
# Kelsey Grantham (kelsey.grantham@monash.edu)

# Trial configuration parameters
clust_per_seq <- c(1, 2, 5)
periods <- c(5, 9)
subjects <- c(10, 100)

# True parameter values
WPICC <- c(0.05, 0.1)
CAC <- c(1.0, 0.8)
theta <- c(0)

parvals <- expand.grid(
  subjects=subjects,
  periods=periods,
  clust_per_seq=clust_per_seq,
  WPICC=WPICC,
  CAC=CAC,
  theta=theta
)

parvals <- parvals %>% arrange(subjects) %>% arrange(periods)

write.csv(parvals, 'parameters.csv', row.names=FALSE)
