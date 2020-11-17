#' @useDynLib EigenR
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
NULL

#' Determinant of a matrix
#' @description Determinant of a real or complex matrix.
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

#' Rank of a matrix
#' @description Rank of a real or complex matrix.
#'
#' @param M a matrix, real or complex
#'
#' @return The rank of \code{M}.
#' @export
Eigen_rank <- function(M){
  stopifnot(is.matrix(M))
  stopifnot(is.numeric(M) || is.complex(M))
  if(is.complex(M)){
    EigenR_rank_cplx(Re(M), Im(M))
  }else{
    EigenR_rank_real(M)
  }
}

Eigen_kernel <- function(M, method = "COD"){
  stopifnot(is.matrix(M))
  stopifnot(is.numeric(M) || is.complex(M))
  method <- match.arg(method, c("COD", "LU"))
  if(is.complex(M)){
    if(method == "COD"){
      parts <- EigenR_kernel_COD_cplx(M)
    }else{
      parts <- EigenR_kernel_LU_cplx(M)
    }
    parts[["real"]] + 1i * parts[["image"]]
  }else{
    if(method == "COD"){
      EigenR_kernel_COD_real(M)
    }else{
      EigenR_kernel_LU_real(M)
    }
  }
}