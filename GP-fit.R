GP.fit <- function(theta=list(lambda=1,sf2=1,sn2=0.05), covFunc, f, Xs, method="cholesky"){
  ##=======================================================
  # Fits GP Regression to the observed data               #
  # Notation mainly adopted from Rasmussen  and Williams  #
  # book 'Gaussian Processes for Machine Learning'        #
  ##=======================================================
  
  x       <- f$x
  y       <- f$y
  n       <- NROW(x)              # Length of the training data
  I       <- diag(1, n)           # Identity matrix
  
  K       <- covFunc(theta, x, x) # Covariance matrix of the training inputs
  noise   <- theta$sn2^2 * I      # Calculate the white noise variance
  
  if (identical(method,"normal")){  # Compute using direct equations
    invXX   <- solve(K + noise)
    E.f     <- covFunc(theta,Xs,x) %*% invXX %*% y
    C.f     <- covFunc(theta,Xs,Xs) - covFunc(theta,Xs,x) %*% invXX %*% covFunc(theta,x,Xs)
    NLML    <- 0.5*t(y) %*% invXX %*% y + 0.5*log(det(invXX)) + 0.5*n*log(2*pi)
  }else if (identical(method,"cholesky")){  # Compute using Cholesky decomposition
    L         <- t(chol(K + noise))
    a         <- solve.cholesky(L, y)       # solve(t(L), solve(L, f$y))
    k.star    <- covFunc(theta,x,Xs)
    E.f       <- t(k.star) %*% a
    v         <- solve(L, k.star)
    C.f       <- covFunc(theta,Xs,Xs) - t(v) %*% v
    NLML      <- 0.5*t(f$y) %*% a + sum(log(diag(L))) + 0.5*n*log(2*pi)
  }
  return(list(E.f=E.f, C.f=C.f, NLML=NLML))
}