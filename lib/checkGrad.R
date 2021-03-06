checkGrad <- function(f, df, e, X){
  ##===============================================================
  # checkGrad checks the derivatives in a function, by comparing  #
  # them to finite difference approximations. We use:             #
  # (f(x+e/2) - f(x-e/2))/e = df + O(e^2)                         #
  #                                                               #
  # Usage:                                                        #
  #   d <- checkgrad(f, df, e, X)                                 #
  # where:                                                        #
  #     f     is the function                                     #
  #     df    is the derivative of the function                   #
  #     e     is small perturbation used for finite differences   #
  #     X     is a vector of data points to compute f and df      #
  #                                                               #
  #     d     differences of the values between f(X) and df(X)    #
  #                                                               #
  ##===============================================================
  
  if (is.function(f) & is.function(df)){
    x <- (f(X+e/2)-f(X-e/2))/e
    y <- df(X)
  }else{
    stop("f and df need to be Objects of Type Function")
  }
  return(abs(x-y))
}