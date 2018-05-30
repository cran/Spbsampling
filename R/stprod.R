#' Standardize a symmetric matrix (distances) to fixed row (column) products.
#'
#' \code{stprod} standardizes the distance matrix to fixed rows and columns products to use \code{\link{pwd}}.
#' The function iteratively constrains a logarithmic transformed matrix to know products, and in order to keep the simmetry of the matrix, at each iteration performs an average with its transpose.
#' When the known products are all equal to a constant (e.g. 1), this method provides a simple and accurate way to scale a distance matrix to a doubly stochastic matrix. The new matrix will not be affected by problems arising from units with different inclusion probabilities, due to not required features of the spatial distribution of the population, such as edge effects and isolated points.
#'
#' @param mat A distance matrix size NxN.
#' @param vec A vector of row (column) constraints.
#' @param differ A scalar with the maximum accepted difference with the constraint.
#' @param niter An integer with the maximum number of iterations.
#' @return Return a distance matrix constrained size NxN.
#' @references
#' \insertRef{BIMJ:BIMJ1785}{Spbsampling}
#'
#' @examples
#' dis<-as.matrix(dist(cbind(simul1$x,simul1$y))) #distance matrix
#' vec<-rep(1,nrow(dis)) #vector of constraints
#' stand_dist<-stsum(dis,vec,1e-15,1000) #standardized matrix


#' @useDynLib Spbsampling
#' @importFrom Rcpp sourceCpp


#'@export
stprod <- function (mat, vec,differ,niter)
{
  mat <- log(mat)
  diag(mat) <- 0
  nsize  <- nrow(mat)
  v      <- 0
  dif    <- 1
  while (dif > differ && (v <- v+1)<niter)
  {
    mat <- (mat-rowSums(mat)/(nsize-1))+vec/(nsize-1)
    mat <- (mat+t(mat))/2
    diag(mat) <- 0
    dif <- max(abs(rowSums(mat)-vec))
    print(c(v,dif))
  }
  mat <- exp(mat)
  diag(mat) <- 0
  mat
}
