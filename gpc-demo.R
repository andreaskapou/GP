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
library(R.utils)
source("gpc-Laplace.R")
source("newton-optimization.R")
sourceDirectory("lib", modifiedOnly=FALSE) # Source the 'lib' directory


##=======================
# Generate data set     #
##=======================
set.seed(1235)    # Set a seed for repeatable plots
n1  <- 80         # Number of data points from each class
n2  <- 40
mu1 <- c(3,2)     # Two mean vectors
mu2 <- c(-3,-2)
S1  <- diag(2)    # Two covariance matrices
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
theta   <- list(lambda=l, sf2=sf2)

tol     <- 1e-6;        # Tolerance for when to stop Newton iterations
covFunc <- "covSE.iso"  # Covariance function to be used
covFunc <- get(covFunc) # Set the string as a variable

lik     <- "cumGauss"   # Likelihood function
lik     <- get(lik)     # Set the string as a variable

##=========================================
# Call the Laplace approximation function #
# for GP Classification                   #
##=========================================
GP <- gpc.Laplace(theta=theta, covfunc=covfunc, lik=lik, tol=tol, x=x, y=y, Xs=Xs)




##=========================================
# Use the binCumGauss Likelihood function #
##=========================================
m1      <- rbinom(80, 13, prob=0.95)
k1      <- rbinom(80, 6,  prob=0.95)
m2      <- rbinom(40, 13, prob=0.95)
k2      <- rbinom(40, 10, prob=0.95)

m       <- c(m1, m2)
k       <- c(k1, k2)
y       <- list(m=m, k=k)
lik     <- "binCumGauss"   # Likelihood function
lik     <- get(lik)     # Set the string as a variable

##=========================================
# Call the Laplace approximation function #
# for GP Classification                   #
##=========================================
GPB <- gpc.Laplace(theta=theta, covfunc=covfunc, lik=lik, tol=tol, x=x, y=y, Xs=Xs)