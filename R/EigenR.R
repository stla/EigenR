#' @useDynLib EigenR
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
NULL

#' Determinant
#'
#' @param M a square matrix, real or complex
#'
#' @return The determinant of \code{M}.
#' @export
Eigen_det <- function(M){
  stopifnot(is.matrix(M) && (nrow(M) == ncol(M)))
  stopifnot(is.numeric(M) || is.complex(M))
  if(is.complex(M)){
    EigenR_det_cplx(Re(M), Im(M))
  }else{
    EigenR_det_real(M)
  }
}
