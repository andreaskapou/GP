binCumGauss <- function(f, m, k, pred){
  ##=====================================================================
  # binCumGauss - Binomial Distributed Cumulative Gaussian likelihood   #
  # function. The expression for log likelihood is                      #
  # log p(y|f) = log(m over k) + k*log Phi(f) + (m-k)*log(1-Phi(f))     #
  #                                                                     #
  # Two modes are provided, for computing likelihoods and derivatives   #
  # and for computing the approximate predictive probability.           #
  #                                                                     #
  # Usage: Phi <- cumGaus(f, m, k)                                      #
  #    or: Phi <- cumGaus(f, m, k, pred)                                #
  # where in Mode 1:                                                    #
  #     f     is a vector of inputs where the CDF needs to be evaluated #
  #     m     is a vector of the total number of trials                 #
  #     k     is a vector of the corresponding successes on each trial  #
  #                                                                     #
  #     lp    is log p(y|f) log likelihood                              #
  #     d1lp  is the 1st derivative                                     #
  #     d2lp  is the 2nd derivative                                     #
  #     d3lp  is the 3rd derivative                                     #
  # where in Mode 2:   ???? not implemented yet                         #
  #     f     is the mean of a Gaussian                                 #
  #     y     is the variance of the Gaussian                           #
  #     pred  is a dummy argument                                       #
  #                                                                     #
  #     pi.pred is the approximate predictive probability               #
  #                                                                     #
  ##=====================================================================
  if (missing(pred)){
    Phi   <- pnorm(f) + 1e-10   # Cumulative Density Function of N(0,1)
    N     <- dnorm(f)           # Density for function f under N(0,1)
    # Compute -> log(choose(m, k)) + k*log(Phi) + (m-k)*log(1-Phi)
    llik  <- dbinom(k, m, Phi, log=TRUE)
    lp    <- sum(llik)          # log p(y|f) log likelihood
    
    d1lp  <- k*N/Phi - (m-k)*N/(1-Phi)    # 1st derivative
    d2lp  <- (-k*(N*(f*Phi+N)/Phi^2)) - (m-k)*(N*(N-f*(1-Phi))/(1-Phi)^2)
    d3lp  <- (-k*(N*(-2*(N^2)-3*N*x*Phi-(Phi^2)*((x^2)-1)))/(Phi^3)) - 
                (m-k)*(N*(N*(2*N-3*x+3*x*Phi)+(x^2)*((Phi^2)-2*Phi+1)-
                                            (Phi^2)+2*Phi-1)/((1-Phi)^3))
    return(list(lp=lp, d1lp=d1lp, d2lp=d2lp, d3lp=d3lp))
  }else{
    pi.pred <- pnorm(f/sqrt(1 + y)) # Approximate predictive probability for 
    # the probit likelihood i.e.
    # p(y*=1|D,t,x*) = N.CDF((m*/sqrt(1+s2*)))
    return(pi.pred)
  }
}