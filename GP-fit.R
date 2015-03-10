GP.fit <- function(f, X, N=50, l=1, s.f=1, s.n=0.1, method="cholesky"){
  ##=======================================================
  # Fits GP Regression to the observed data               #
  # Notation mainly adopted from Rasmussen  and Williams  #
  # book `Gaussian Processes for Machine Learning'        #
  ##=======================================================
  
  # Calculate covariance function
  Cov <- function(X, Y) outer(X, Y, FUN=SE, l, s.f)
  #Cov <- function(X, Y) outer(X, Y, FUN=Wiener)
  ##=======================
  # Initialize parameters #
  ##=======================
  n       <- length(f$x)  # Length of the training data
  I       <- diag(1, n)   # Identity matrix
  logLik  <- 0            # Log margninal likelihood
  
  ##=================
  # Fit GP to data  #
  ##=================
  K <- Cov(f$x,  f$x)  # Covariance matrix of the training inputs
  
  if (identical(method,"normal")){ # Compute using direct equations
    invCovXX  <- solve(K + s.n^2 * I)
    E.f       <- Cov(X, f$x) %*% invCovXX %*% f$y
    C.f       <- Cov(X, X) - Cov(X, f$x) %*% invCovXX %*% Cov(f$x, X)
    logLik    <- -.5 * t(f$y) %*% invCovXX %*% f$y - .5 * log(det(invCovXX)) - 
                                                                n * log(2*pi)/2
  }else if (identical(method,"cholesky")){ # Compute using Cholesky decomposition
    L       <- t(chol(K + s.n^2 * I))
    a       <- solve(t(L), solve(L, f$y))
    k.star  <- Cov(f$x, X)
    E.f     <- t(k.star) %*% a
    v       <- solve(L, k.star)
    C.f     <- Cov(X, X) - t(v) %*% v
    logLik  <- -.5 * t(f$y) %*% a - sum(log(diag(L))) - n * log(2*pi)/2
  }
  
  # Create a lot of samples.  We could of course
  # simply plot a +/- 2 standard deviation confidence interval.
  values <- matrix(rep(0,length(X)*N), ncol=N)
  for (i in 1:N) {
    
    # Sample from a multivariate normal distribution
    #values[,i] <- C.f %*% rnorm(length(E.f)) + E.f
    C.f <- C.f
    values[,i] <- rmvnorm(1, mean=E.f, sigma=C.f, method="svd")
    #values[,i] <- t(rmvnorm(1, mean=rep(0, length(E.f)), sigma=C.f, method="svd")) + E.f
  }
  values <- cbind(x=X,as.data.frame(values))
  values <- melt(values,id="x")
  return(list(E.f=E.f, C.f=C.f, logLik=logLik, values=values, K=K))
}

##===============================================
# Functions to compute the covariance function  #
##===============================================
# Calculate Squared Exponential Function
SE  <- function(X1, X2, l, s.f){
  (s.f^2) * exp(-.5 * ((X1 - X2)^2) / (l^2))
}
Wiener <- function(X1, X2){
  min(X1,X2) - X1*X2
}