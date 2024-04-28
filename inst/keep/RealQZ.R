#' Real QZ decomposition
#' @description Real QZ decomposition of a pair of square matrices.
#'
#' @param A,B real square matrices with the same size
#'
#' @return A list with the \code{Q}, \code{Z}, \code{S} and \code{T} matrices.
#' @export
#' @details See \href{https://eigen.tuxfamily.org/dox/classEigen_1_1RealQZ.html}{Eigen::RealQZ}.
#'
#' @examples
#' library(EigenR)
#' A <- toeplitz(c(1, 2, 3))
#' B <- cbind(c(3, 2, 3), c(1, 1, 1), c(5, 0, -2))
#' qz <- Eigen_realQZ(A, B)
#' Q <- qz$Q
#' Z <- qz$Z
#' S <- qz$S
#' T <- qz$T
#' # check decomposition:
#' A - Q %*% S %*% Z # should be zero matrix
#' B - Q %*% T %*% Z # should be zero matrix
#' # check orthogonality of Q and Z:
#' tcrossprod(Q) # should be identity matrix
#' tcrossprod(Z) # should be identity matrix
Eigen_realQZ <- function(A, B) {
  stopifnot(isSquareMatrix(A), isSquareMatrix(B))
  stopifnot(isReal(A), isReal(B))
  if(nrow(A) != nrow(B)) {
    stop("The matrices `A` and `B` must have the same size.")
  }
  EigenR_realQZ(A, B)
}

