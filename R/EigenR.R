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

#' Kernel of a matrix
#' Kernel (null-space) of a real or complex matrix.
#'
#' @param M a matrix, real or complex
#' @param method one of \code{"COD"} or \code{"LU"}; the faster method depends 
#'   on the size of the matrix
#'
#' @return A basis of the kernel of \code{M}. With \code{method = "COD"}, the 
#'   basis is orthonormal, while it is not with \code{method = "LU"}.
#' @export
#'
#' @examples xx
Eigen_kernel <- function(M, method = "COD"){
  stopifnot(is.matrix(M))
  stopifnot(is.numeric(M) || is.complex(M))
  method <- match.arg(method, c("COD", "LU"))
  if(is.complex(M)){
    if(method == "COD"){
      parts <- EigenR_kernel_COD_cplx(Re(M), Im(M))
    }else{
      parts <- EigenR_kernel_LU_cplx(Re(M), Im(M))
    }
    parts[["real"]] + 1i * parts[["imag"]]
  }else{
    if(method == "COD"){
      EigenR_kernel_COD_real(M)
    }else{
      EigenR_kernel_LU_real(M)
    }
  }
}

#' Range of a matrix
#' Range (column-space, image, span) of a real or complex matrix.
#'
#' @param M a matrix, real or complex
#' @param method one of \code{"LU"} or \code{"QR"}
#'
#' @return A basis of the range of \code{M}. With \code{method = "LU"}, the 
#'   basis is not orthonormal, while it is with \code{method = "QR"}.
#' @export
#'
#' @examples xx
Eigen_range <- function(M, method = "QR"){
  stopifnot(is.matrix(M))
  stopifnot(is.numeric(M) || is.complex(M))
  method <- match.arg(method, c("LU", "QR"))
  if(is.complex(M)){
    if(method == "QR"){
      parts <- EigenR_image_QR_cplx(Re(M), Im(M))
    }else{
      parts <- EigenR_image_LU_cplx(Re(M), Im(M))
    }
    parts[["real"]] + 1i * parts[["imag"]]
  }else{
    if(method == "QR"){
      EigenR_image_QR_real(M)
    }else{
      EigenR_image_LU_real(M)
    }
  }
}
