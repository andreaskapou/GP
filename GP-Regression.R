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
source("GP-fit.R")

##=======================
# Initialize parameters #
##=======================
set.seed(12345)   # Set a seed for repeatable plots
N       <- 50     # Number of samples
l       <- 1      # Length-scale parameter
s.f     <- 1      # Singal variance
s.n     <- .1     # Noise variance
# Define the points at which we want to define the functions
X       <- seq(-5, 5, len=100)
# Assume that we have some known data points
f       <- data.frame(x=c(-4,-3,-1,0,2,4, 5),
                      y=c(-1,1,1,1,-1,-2, 0))
method <- "normal" 

##=================================
# Call the GP Regression function #
##=================================
GP      <- GP.fit(f, X, N, l, s.f, s.n, method)
m.star  <- GP$E.f
values  <- GP$values

# Plot the result, including error bars on the observed points
gg2 <- ggplot(values, aes(x=x,y=value)) + 
  geom_line(aes(group=variable), colour="grey80") +
  geom_line(data=NULL,aes(x=X,y=m.star),colour="red", size=1) + 
  geom_errorbar(data=f,aes(x=x,y=NULL,ymin=y-2*s.n, ymax=y+2*s.n), width=0.2) +
  geom_point(data=f,aes(x=x,y=y)) +
  theme_bw() +
  scale_y_continuous(lim=c(-3,3), name="output, f(x)") +
  xlab("input, x")