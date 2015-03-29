covSE.iso  <- function(theta, X1, X2, der){
  ##=========================================================
  # Squared Exponential covariance function with isotropic  #
  # distance measure. The covariance is parameterized as:   #
  # k(x,x') = sf2 * exp(-t((x-x'))*inv(P)*(x-x')/2)         #
  # where the P matrix is lambda^2 times the unit matrix    #
  # and sf2 is the signal variance.                         #
  # The hyperparameters are:                                #
  #     theta = [ lambda, sf2 ]                             #
  #                                                         #
  # If 'der' parameter is passed, the derivatives of the SE #
  # function are computed w.r.t. the hyperparameters.       #
  # So when:                                                #
  #     der == 1: derivatives wrt lambda                    #
  #     der == 2: derivatives wrt sf2                       #
  ##=========================================================
  
  # Unpack the parameters
  if (length(theta) < 2)
    stop("Theta needs to be provided at least with 2 hyperparameters.")
  if (nargs() < 3)
    stop("Function requires 3 parameters.")
  lambda    <- theta$lambda^2 # Length-scale parameter
  sf2       <- theta$sf2^2    # Signal variance parameter
  
  # Divide by the square of length-scale lambda
  X1        <- as.matrix(X1) / lambda
  X2        <- as.matrix(X2) / lambda
  
  # Split by rows (i.e. to take D-dimensional data as a unit)
  X1.rows   <- split(X1, row(X1))
  X2.rows   <- split(X2, row(X2))
  
  # Vectorize the 'sq.dist' function so as to be able to call
  # the 'outer' function on sets of vectors.
  Vec.Dist  <- Vectorize(sq.dist)
  
  # Compute distance of the D-dimensional data points
  K <- outer(X1.rows, X2.rows, FUN=Vec.Dist)
  if (missing(der)){
    return(sf2 * exp(-0.5*K))
  }else{
    if (der==1){
      return(sf2 * exp(-0.5*K) * K)
    }else if (der==2){
      return(2 * sf2 * exp(-0.5*K))
    }else{
      stop("Unknown hyperparameter derivative")
    }
  }
}