gpc.Laplace <- function(theta=list(lambda=1,sf2=1,sn2=0.001),covfunc,lik,x,y,Xs,tol=1e-6){
  ##=======================================================================
  # gpc.Laplace - Laplace's approximation for binary Gaussian process     #
  # classification. Two modes are possible: training or testing:          #
  #                                                                       #
  # 1) If no test cases are supplied, then the approximate negative log   #
  #  marginal likelihood and its partial derivatives wrt the hyperparam   #
  #  is computed; this mode is used to fit the hyperparameters.           #
  # 2) If test cases are given, then test set predictive probabilities    #
  #  are returned.                                                        #
  #                                                                       #
  # The program is flexible in allowing several different likelihood      #
  # functions and a multitude of covariance functions.                    #
  #                                                                       #
  # usage: [NLML dnlml] <- gpc.Laplace(theta,covfunc,lik,x,y, tol)        #
  #    or: [p E.f C.f NLML] <- gpc.Laplace(theta,covfunc,lik,x,y,Xs, tol) #
  #                                                                       #
  # where:                                                                #
  #   theta     is a (column) vector of hyperparameters                   #
  #   covfunc   is the name of the covariance function (see below)        #
  #   lik       is the name of the likelihood function (see below)        #
  #   x         is a n by D matrix of training inputs                     #
  #   y         is a (column) vector (of size n) of binary +1/-1 targets  #
  #   xstar     is a nn by D matrix of test inputs                        #
  #   nlml      is the returned value of the neg log marg likelihood      #
  #   dnlml     is a (column) vector of partial derivatives of the neg    #
  #               log marginal likelihood wrt each log hyperparameter     #
  #   p         is a (column) vector (of length nn) of predictive probabs #
  #   E.f       is a vector (of length nn) of predictive latent means     #
  #   C.f       is a vector (of length nn) of predictive latent variances #
  #                                                                       #
  # The shape of the likelihood function is given by the "lik" input to   #
  # the function, specifying the name of the likelihood function. At      #
  # present the only implemented likelihood function is:                  #
  #                                                                       #
  #   cumGauss   the cumulative Gaussian (error function)                 #
  #                                                                       #
  # The function can conveniently be used with the "minimize" function to #
  # train a Gaussian process (NOT IMPLEMENTED YET).                       #
  #                                                                       #
  # Implementation follows Algorithm 3.2 p.47 from Rasmussen & Williams   #
  # book 'Gaussian Processes for Machine Learning'                        #
  # Adapted from copyright 2006 Carl Edward Rasmussen version in matlab.  #
  ##=======================================================================
  n       <- NROW(x)                # Length of the training data
  I       <- diag(1, n)             # Identity matrix
  K       <- covFunc(theta, x, x)   # Covariance matrix of the training data
  
  ##=========================================
  # Call the Newton's optimization function #
  ##=========================================
  NO      <- newton.optimization(K, y, lik, tol)
  
  ##=========================================
  # Laplace Approximation to the posterior  #
  ##=========================================
  W       <- (-NO$Phi$d2lp)                 # W = -DDlog p(y|f)
  sW      <- sqrt(W)                        # Compute W^1/2
  L       <- t(chol(I + sW %*% t(sW) * K))  # B = I + W^1/2*K*W^1/2
  # Approximate negative log marginal likelihood 
  NLML    <- t(NO$a)%*%NO$f/2 - NO$Phi$lp + sum(log(diag(L)))
  if (missing(Xs)){
    DE    <- vector(length=length(theta))
            # Z = W^{1/2}*(I + W^{1/2}*K*W^{1/2})^{-1}*W^{1/2} using Cholesky
    Z     <- matrix(sW, nrow=length(sW), ncol=n)*solve.cholesky(L, diag(as.vector(sW)))
            # C = L\(W^{1/2}*K)
    C     <- solve(L, matrix(sW, nrow=length(sW), ncol=n)*K)
            # s2 = -1/2 * [(K^{-1} + W)^{-1}]_{ii} * der^{3}(log p(y|f))
    s2    <- -0.5*(diag(K) - as.matrix(colSums(C^2)))*(-NO$Phi$d3lp) # Should the '-' be there?
    for (j in 1: length(theta)){
      C     <- covFunc(theta, x, x, j)   # Derivative matrix wrt hyperparameters
            # s1 = 1/2*f'*K^{-1}*C*K^{-1}*f - 1/2*tr[(W^{-1}+K)^{-1} * C]
      s1    <- 0.5*t(NO$a) %*% C %*% NO$a - 0.5*sum(colSums(Z*C))
            # b = C * der(log p(y|f))
      b     <- C %*% NO$Phi$d1lp
            # (I + KW)^{-1} * b
      s3    <- b - K %*% (Z %*% b)
      DE[j] <- -s1 - t(s2) %*% s3     # Take the negative Derivative of LML
    }
    return(list(NLML=NLML, DE=DE))
  }else{
    # Compute predictive probabilities
    k.star  <- covFunc(theta,x,Xs)
    E.f     <- t(k.star) %*% NO$Phi$d1lp      # Latent means
              # v = L\(W^{1/2}*k(x*))
    v       <- solve(L, matrix(sW,nrow=length(sW),ncol=NROW(Xs))*k.star)
    #C.f    <- covFunc(theta,Xs,Xs)-t(v)%*%v  # Impractical for large datasets
    Kss     <- rep(1, NROW(Xs))
    C.f     <- as.matrix(Kss) - as.matrix(colSums(v * v)) # Latent variances
    p       <- lik(E.f, C.f, 1)               # Average predictive probability
    
    return(list(NLML=NLML, E.f=E.f, C.f=C.f, p=p))
  }
}