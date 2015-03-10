solve.cholesky <- function(L, B){
  ##=========================================================
  # Solve linear equations from the Cholesky factorization. # 
  # Solve A*X = B for X, where A is square, symmetric,      #
  # positive definite. The input to the function is L the   #
  # Cholesky decomposition of A and the matrix B.           #
  # Example: X <- solve.cholesky(t(chol(A)), B)             #
  ##=========================================================
  if (is.matrix(L)){
    if ( (nrow(L) != ncol(L)) | (nrow(L)) != length(B) )
      stop("Wrong sizes of matrix arguments.")
  }else{
    stop("First argument should be a square matrix.")
  }
  x <- solve(t(L), solve(L, B))
  return(x)
}