cumGauss <- function(f, y, pred){
  ##=====================================================================
  # cumGauss - Cumulative Gaussian likelihood function. The expression  #
  # for likelihood is cumGauss(t) = normcdf(t) = (1+erf(t/sqrt(2)))/2.  #
  #                                                                     #
  # Two modes are provided, for computing likelihoods and derivatives   #
  # and for computing the approximate predictive probability.           #
  #                                                                     #
  # Usage: Phi <- cumGaus(f, y)                                         #
  #    or: Phi <- cumGaus(f, y, pred)                                   #
  # where in Mode 1:                                                    #
  #     f     is a vector of inputs where the CDF needs to be evaluated #
  #     y     is a vector of outputs (-1 or +1)                         #
  #                                                                     #
  #     lp    is log p(y|f) log likelihood                              #
  #     d1lp  is the 1st derivative                                     #
  #     d2lp  is the 2nd derivative                                     #
  #     d3lp  is the 3rd derivative                                     #
  # where in Mode 2:                                                    #
  #     f     is the mean of a Gaussian                                 #
  #     y     is the variance of the Gaussian                           #
  #     pred  is a dummy argument                                       #
  #                                                                     #
  #     pi.pred is the approximate predictive probability               #
  #                                                                     #
  ##=====================================================================
  if (missing(pred)){
    Phi   <- pnorm(f*y) + 1e-10 # Cumulative Density Function of N(0,1)
    N     <- dnorm(f)           # Density for function f under N(0,1)
    lp    <- sum(log(Phi))      # log p(y|f) log likelihood
    d1lp  <- y*N / Phi          # 1st derivative
    d2lp  <- (- N^2 / Phi^2) - (y*f*N)/Phi  # 2nd derivative
    d3lp  <- 3*f*N^2/Phi^2 + 2*y*N^3/Phi^3 + y*N*(f^2-1)/Phi # 3rd derivative
    return(list(lp=lp, d1lp=d1lp, d2lp=d2lp, d3lp=d3lp))
  }else{
    pi.pred <- pnorm(f)   # Need to compute the average likelihood
    return(pi.pred)
  }
}