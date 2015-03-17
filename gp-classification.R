# Load in the required libraries for data manipulation
# and multivariate normal distribution
require(MASS)
require(plyr)
require(reshape2)
require(ggplot2)

# Set a seed for repeatable plots
#set.seed(123457)
n.samples <- 1
l       <- 1      # Length-scale parameter
s.f     <- 1      # Singal variance
s.n     <- .1     # Noise variance

##===============================================
# Functions to compute the covariance function  #
##===============================================
# Calculate Squared Exponential Function
SE  <- function(X1, X2, l, s.f){
  (s.f^2) * exp(-.5 * ((X1 - X2)^2) / (l^2))
}
# Calculate covariance function
Cov <- function(X, Y) outer(X, Y, FUN=SE, l, s.f)

sigmoid <- function(z) {
  return(1/(1 + exp(-z)))
}

# Define the points at which we want to define the functions
x.star <- seq(-5,5,len=100)
# Calculate the covariance matrix using Gaussian kernel
K <- Cov(x.star,x.star)

values <- matrix(rep(0,length(x.star)*n.samples), ncol=n.samples)
for (i in 1:n.samples) {
  # Each column represents a sample from a multivariate normal distribution
  # with zero mean and covariance sigma
  values[,i] <- rmvnorm(1, mean=rep(0, length(x.star)), sigma=K, method="svd")
  #values[,i] <- mvrnorm(1, rep(0, length(x.star)), sigma)
}
cl.values <- sigmoid(values[,1])

values  <- cbind(x=x.star,as.data.frame(values))
values  <- melt(values,id="x")

cl.values <- cbind(x=x.star,as.data.frame(cl.values))
cl.values <- melt(cl.values,id="x")

# Plot the result
fig2a <- ggplot(values,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-2, ymax=2, fill="grey80") +
  geom_line(aes(group=variable)) +
  theme_bw() +
  scale_y_continuous(lim=c(-2.5,2.5), name="output, f(x)") +
  xlab("input, x")


# Plot the result
fig2b <- ggplot(cl.values,aes(x=x,y=value)) +
  geom_rect(xmin=-Inf, xmax=Inf, ymin=-.5, ymax=1.5, fill="grey80") +
  geom_line(aes(group=variable)) +
  theme_bw() +
  scale_y_continuous(lim=c(-.2,1.2), name="output, f(x)") +
  xlab("input, x")