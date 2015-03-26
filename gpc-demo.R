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
source("gp-classification.R")
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
l       <- 1     # Length-scale parameter
sf2     <- 1     # Singal variance
theta   <- list(lambda=l, sf2=sf2)
covFunc <- "covSE.iso"  # Covariance function to be used
covFunc <- get(covFunc) # Set the string as a variable

n       <- NROW(x)      # Length of the training data
I       <- diag(1, n)   # Identity matrix

tol     <- 1e-6;        # Tolerance for when to stop the Newton iterations
K       <- covFunc(theta, x, x)   # Covariance matrix of the training data

##===================================
# Call the Newton's method function #
##===================================
A.mode  <- newton.optimization(K, y, tol)
Phi     <- A.mode$Phi
a       <- A.mode$a
f       <- A.mode$f


##=======================
# Laplace Approximation #
##=======================
W       <- (-Phi$d2lp)  # W = -DDlog p(y|f)
sW      <- sqrt(W)      # Compute W^1/2
L       <- t(chol(I + sW %*% t(sW) * K))   # B = I + W^1/2*K*W^1/2

NLML    <- t(a)%*%f/2 - Phi$lp + sum(log(diag(L)))  # Approx neg log marginal lik 

# Compute predictive probabilities
k.star  <- covFunc(theta,x,Xs)
E.f     <- t(k.star) %*% Phi$d1lp
sW.k    <- matrix(sW, nrow=length(sW), ncol=NROW(Xs)) * k.star
v       <- solve(t(L), sW.k)

#C.f     <- covFunc(theta,Xs,Xs) - t(v)%*%v #impractical for large datasets
Kss     <- rep(theta$sf2^2 + 1, NROW(Xs))
C.f     <- as.matrix(Kss) - as.matrix(colSums(v * v))

