# Chapter 2 of Rasmussen and Williams's book `Gaussian Processes
# for Machine Learning' provides a detailed explanation of the
# math for Gaussian process regression.

##===============================
# Set the working directory and #
# load any required libraries   #
##===============================
cur.dir <- dirname(parent.frame(2)$ofile)
setwd(cur.dir)
require(MASS)
require(plyr)
require(reshape2)
require(ggplot2)
require(mvtnorm)
source("gp-regression.R")
source("solve-cholesky.R")
source("covSE-iso.R")
source("sq-dist.R")

##=======================
# Initialize parameters #
##=======================
set.seed(12345)   # Set a seed for repeatable plots
N       <- 50     # Number of samples
l       <- 1.1    # Length-scale parameter
sf2     <- 1.1    # Singal variance
sn2     <- .03    # Noise variance

Xs      <- seq(-8, 8, len=100)  # Test data points
covFunc <- "covSE.iso"  # Covariance function to be used
covFunc <- get(covFunc) # Set the string as a variable
method  <- "cholesky"
theta   <- list(lambda=l, sf2=sf2, sn2=sn2)

# Assume that we have some known data points
n       <- 20
x       <- as.vector(15 * (runif(n) - 0.5))
y       <- as.vector(chol(covFunc(theta, x, x)) %*% rnorm(n))
f       <- data.frame(x=x,
                      y=y)
f       <- data.frame(x=c(-4,-3,-1,0,2,4, 5),
                      y=c(-1,1,1,1,-1,-2, 0))

##=================================
# Call the GP Regression function #
##=================================
GP      <- GP.fit(theta, covFunc, f, Xs, method=method)
mu      <- GP$E.f
S2      <- GP$C.f

# If we computed only diagonal covariance instead of all test point covariances
if (dim(S2)[2] == 1){
  S2 <- diag(as.vector(S2))
}

# Create a lot of samples.  We could of course
# simply plot a +/- 2 standard deviation confidence interval.
values <- matrix(rep(0,length(Xs)*N), ncol=N)
for (i in 1:N) { # Sample from a multivariate normal distribution
  #values[,i] <- S2 %*% rnorm(length(mu)) + mu
  values[,i] <- rmvnorm(1, mean=mu, sigma=S2, method="svd")
  #values[,i] <- t(rmvnorm(1, mean=rep(0,length(mu)),sigma=S2,method="svd"))+mu
}
values <- cbind(x=Xs,as.data.frame(values))
values <- melt(values,id="x")

# Plot the result, including error bars on the observed points
gg2 <- ggplot(values, aes(x=x,y=value)) + 
  geom_line(aes(group=variable), colour="grey80") +
  geom_line(data=NULL,aes(x=Xs,y=mu),colour="red", size=1) + 
  geom_errorbar(data=f,aes(x=x,y=NULL,ymin=y-2*sn2, ymax=y+2*sn2), width=0.2) +
  geom_point(data=f,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(lim=c(-5,5), name="output, f(x)") +
  xlab("input, x")


GPM      <- GP.fit(theta, covFunc, f, method=method)