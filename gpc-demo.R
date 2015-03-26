# Chapter 3 of Rasmussen and Williams's book `Gaussian Processes
# for Machine Learning' provides a detailed explanation of the
# math for Gaussian process classification

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
source("gpc-Laplace.R")
source("solve-cholesky.R")
source("covSE-iso.R")
source("sq-dist.R")
source("meshgrid.R")
source("cumGauss.R")
source("newton-optimization.R")

##=======================
# Generate data set     #
##=======================
set.seed(12345)   # Set a seed for repeatable plots
n1  <- 80         # Number of data points from each class
n2  <- 40
mu1 <- c(1,0)     # Two mean vectors
mu2 <- c(-1,0)
S1  <- diag(2)    # Two covatriance matrices
S2  <- matrix(c(1, 0.95, 0.95, 1), nrow=2, ncol=2)
x1  <- mvrnorm(n1, mu=mu1, Sigma=S1)
x2  <- mvrnorm(n2, mu=mu2, Sigma=S2)

x   <- rbind(x1, x2)                # Inputs of training data
y   <- c(rep(-1, n1), rep(1, n2))   # Outputs of training data
t1  <- seq(from=-4, to=4, by=0.1)   # Test data
t   <- meshgrid(t1, t1)
Xs  <- cbind(as.vector(t$x),  as.vector(t$y))


##=======================
# Initialize parameters #
##=======================
l       <- 1      # Length-scale parameter
sf2     <- 1      # Singal variance
sn2     <- 0.01   # Noise variance
theta   <- list(lambda=l, sf2=sf2, sn2=sn2)

tol     <- 1e-6;        # Tolerance for when to stop Newton iterations
covFunc <- "covSE.iso"  # Covariance function to be used
covFunc <- get(covFunc) # Set the string as a variable

lik     <- "cumGauss"   # Likelihood function
lik     <- get(lik)     # Set the string as a variable

##========================================#
# Call the Laplace approximation function #
# for GP Classification                   #
##=========================================
GP <- gpc.Laplace(theta, covfunc, lik, x, y, Xs, tol)