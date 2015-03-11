sq.dist <- function(x1, x2){
  ##=======================================================
  # Function to compute the squared distances between two #
  # sets of points, vectors or matrices.                  #
  # Example: X <- sq.dist(2, 4) would return X=4          #
  ##=======================================================
  return( sum((x1 - x2)^2) )
}