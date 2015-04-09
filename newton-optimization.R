newton.optimization <- function(K, y, lik, tol){
  ##=====================================================================
  # Newton's Method for finding the mode of the posterior distribution  #
  # which is log-concave and thus contains a unique maximum.            #
  #                                                                     #
  # Thus, we need to find f^{hat} = argmax_{f} p(f|X,y), that is in     #
  # the Newton iterations f^{new} := f - (DD Psi(f))^{-1} D Psi(f).     #
  #                                                                     #
  # The found mode will be used as the mean of the Gaussian when doing  #
  # Laplace approximation and covariance matrix will be approximated by #
  # the curvature at the mode, which is equal to the inverse Hessian    #
  #                                                                     #
  # Implementation follows Algorithm 3.1 p.46 from Rasmussen & Williams #
  # book 'Gaussian Processes for Machine Learning'                      #
  #                                                                     #
  # Input:                                                              #
  #     K   is covariance matrix on training points                     #
  #     y   is an object of size n of target values and the structure   #
  #           depends on the likelihood function that will be used.     #
  #           E.g if 'cumGauss' is used y will be a (column) vector     #
  #                 (of size n) of binary +1/-1 targets                 #
  #               if 'binCumGauss is used y will be an object of size   #
  #                 n with each entry being a tuple (m, k)              #
  #     lik is the likelihood function                                  #    
  #     tol is the tolerance for when to stop the Newton iterations     #
  # Output:                                                             #
  #     Phi is an object containing the likelihood and its derivatives  #
  #     a   needed for later computations since a = K^{-1} * f          #
  #     f   is the posterior mode                                       #
  #                                                                     #
  ##=====================================================================
  n       <- NROW(K)        # Length of the training data
  I       <- diag(1, n)     # Identity matrix
  a = f   <- rep(0, n)      # Initial points for f and a
  
  Phi     <- lik(f, y)        # Compute likelihood function and its derivatives
  Psi_new <- -((n*log(2))^2)  # Objective initial value
  Psi_old <- (-Inf)           # Make sure Newton iteration starts
  while (Psi_new-Psi_old > tol){  # Begin Newton's iterations
    Psi_old <- Psi_new
    a_old   <- a                  # In case objective does not increase
    
    W       <- (-Phi$d2lp)                            # W = -DDlog p(y|f)
    sW      <- sqrt(W)                                # Compute W^{1/2}
    L       <- t(chol(I + sW %*% t(sW) * K))          # B = I + W^{1/2}*K*W^{1/2}
    b       <- W*f + Phi$d1lp                         # b = Wf + Dlog p(y|f)
    a       <- b - sW * solve.cholesky(L, sW*(K%*%b)) # b - W^{1/2}*B^-1*W^{1/2}*K*b
    f       <- K %*% a                                # Update f (a = K^{-1}*f)
    Phi     <- lik(f, y)                              # Update Phi using f_new
    
    Psi_new <- (-t(a)%*%f/2) + Phi$lp # objective: -1/2 a'*f + log p(y|f)
                                      #       i.e. -1/2 f'*K^{-1}*f + log p(y|f)
    i         <- 0
    while ((i < 10) & (Psi_new < Psi_old)){ # If objective didn't increase ...
      a       <- (a_old + a)/2                # Reduce step size by half
      f       <- K %*% a
      Phi     <- lik(f, y)
      Psi_new <- (-t(a) %*% f/2) + Phi$lp
      i       <- i+1                          # Repeat at most 10 times
    }
  }
  return(list(Phi=Phi, a=a, f=f))
}