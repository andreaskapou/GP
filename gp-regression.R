GP.fit <- function(theta=list(lambda=1,sf2=1,sn2=0.05),covFunc,f,Xs,method="cholesky"){
  ##=================================================================
  # Gaussian process regression implementation. A valid covariance  #
  # function can be given as a parameter. Two modes are possible:   # 
  # 1) training: if no test data (Xs) are given, function returns   #
  #       the negative log likelihood and its partial derivatives   #
  #       with respect to the hyperparameters(NOT DONE YET); this   #
  #       mode is used to fit the hyperparameters.                  #
  # 2) prediction: If test data are given, then (marginal) Gaussian #
  #       predictions are computed, whose mean and variance are     #
  #       returned. Note that in cases where covariance function    # 
  #       has noise contributions, the variance returned in S2 is   #
  #       for noisy test targets; if you want the variance of the   #
  #       noise-free latent function, you must substract the noise  #
  #       variance.                                                 #
  # Also, there are two modes of implementing the matrix inversion  #
  # one is direct computation and the other uses the Cholesky       #
  # factorization. Better to use the Cholesky factorization         # 
  #                                                                 #
  # usage: GP$NLML <- gpr(theta, covFunc, f, method)                #
  #    or: (GP$E.f, GP$C.f) <- gpr(theta, covFunc, f, Xs, method)   #
  #                                                                 #
  # where:                                                          #
  #   theta   is a list of hyperparameters                          #
  #   covFunc is the covariance function                            #
  #   f       is a list of D-dim training inputs and 1-dim targets  #
  #   y       is a (column) vector (of size n) of targets           #
  #   Xs      is a nn by D matrix of test inputs                    #
  #   method  is the method used to compute the matrix inversion    #
  #                                                                 #
  #   NLML    is the value of the negative log marginal likelihood  #
  #   DNLML   is a (column) vector of partial derivatives of the    #
  #               negative log marginal likelihood wrt each log     #
  #               hyperparameter                                    #
  #   E.f     is a (column) vector (of size nn) of prediced means   #
  #   C.f     is a (column) vector (of size nn) of predicted var    #
  #                                                                 #
  # Adapted from (C) copyright 2006 by Carl Edward Rasmussen        #
  # version in matlab.                                              #
  # Notation mainly follows from Rasmussen and Williams's book      #
  # 'Gaussian Processes for Machine Learning'                       #
  #                                                                 #
  ##=================================================================
  x       <- f$x
  y       <- f$y
  n       <- NROW(x)                # Length of the training data
  I       <- diag(1, n)             # Identity matrix
  
  K       <- covFunc(theta, x, x)   # Covariance matrix of the training inputs
  noise   <- theta$sn2^2 * I        # Calculate the white noise variance
  
  if (identical(method,"normal")){  # Compute using direct equations
    invXX     <- solve(K + noise)
    if (missing(Xs)){ # If no test points, just compute the marginal
      NLML    <- -0.5*t(y) %*% invXX %*% y - 0.5*log(det(invXX)) - 0.5*n*log(2*pi)
    }else{
      E.f     <- covFunc(theta,Xs,x) %*% invXX %*% y
      C.f     <- covFunc(theta,Xs,Xs)-covFunc(theta,Xs,x)%*%invXX%*%covFunc(theta,x,Xs)
    }
  }else if (identical(method,"cholesky")){  # Compute using Cholesky decomposition
    L         <- t(chol(K + noise))
    a         <- solve.cholesky(L, y)       # solve(t(L), solve(L, f$y))
    NLML      <- -0.5*t(y) %*% a - sum(log(diag(L))) - 0.5*n*log(2*pi)
    if (missing(Xs)){ # If no test points, just compute the marginal
      ##### NOT IMPLEMENTED YET
      DE      <- vector(length=length(theta))
      W       <- 
      for (j in 1: (length(theta)-1)){ # Compute only the first two parameter derivatives
        C     <- covFunc(theta, x, x, j)
      }
    }else{
      k.star  <- covFunc(theta,x,Xs)
      E.f     <- t(k.star) %*% a            # Latent means
      v       <- solve(L, k.star)
      #C.f    <- covFunc(theta,Xs,Xs) - t(v)%*%v #impractical for large datasets
      Kss     <- rep(theta$sn2^2 + 1, NROW(Xs))
      C.f     <- as.matrix(Kss) - as.matrix(colSums(v * v)) # Latent variances
    }
  }
  if (missing(Xs))
    return(list(NLML=NLML))
  else
    return(list(E.f=E.f, C.f=C.f))
}