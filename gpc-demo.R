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


##=======================
# Generate data set #
##=======================
set.seed(12345)   # Set a seed for repeatable plots
n1 <- 80          # Number of data points from each class
n2 <- 40

mu1 <- c(1,0)     # Two mean vectors
mu2 <- c(-1,0)

S1 <- diag(2)     # Two covatriance matrices
S2 <- matrix(c(1, 0.95, 0.95, 1), nrow=2, ncol=2)

x1 <- mvrnorm(n1, mu=mu1, Sigma=S1)
x2 <- mvrnorm(n2, mu=mu2, Sigma=S2)

x <- rbind(x1, x2)  # Inputs of training data
y <- c(rep(-1, n1), rep(1, n2)) # Outputs of training data
  
t1 <- seq(from=-4, to=4, by=0.5)   # Test data
t <- meshgrid(t1, t1)
t <- cbind(as.vector(t$x),  as.vector(t$y))

##=======================
# Initialize parameters #
##=======================
l       <- 1     # Length-scale parameter
sf2     <- 1     # Singal variance

covFunc <- "covSE.iso"  # Covariance function to be used
covFunc <- get(covFunc) # Set the string as a variable

theta   <- list(lambda=l, sf2=sf2)



cumGauss <- function(f, y){
  Phi <- pnorm(f*y) + 1e-10 # Cumulative Density Function of N(0,1)
  N   <- dnorm(f)           # Density for function f under N(0,1)
  lp  <- sum(log(Phi))      # log p(y|f)  log likelihood
  d1lp <- y*N / Phi         # 1st derivative
  d2lp <- - N^2 / Phi^2 - (y*f*N)/Phi  # 2nd derivative
  d3lp <- 3*f*N^2/Phi^2 + 2*y*N^3/Phi^3 + y*N*(f^2-1)/Phi # 3rd derivative
  
  return(list(lp=lp, d1lp=d1lp, d2lp=d2lp, d3lp=d3lp))
}




















