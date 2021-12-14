#' @useDynLib EigenR
#' @importFrom Rcpp evalCpp
NULL

#' Determinant of a matrix
#' @description Determinant of a real or complex matrix.
#' 
#' @param M a square matrix or \code{\link{SparseMatrix}}, real or complex
#'
#' @return The determinant of \code{M}.
#' @export
#' 
#' @examples set.seed(666)
#' M <- matrix(rpois(25, 1), 5L, 5L)
#' Eigen_det(M)
#' # determinants of complex matrices are supported:
#' Eigen_det(M + 1i * M)
#' # as well as determinants of sparse matrices:
#' Eigen_det(asSparseMatrix(M))
#' Eigen_det(asSparseMatrix(M + 1i * M))
Eigen_det <- function(M){
  if(inherits(M, "SparseMatrix")){
    stopifnot(M[["nrows"]] == M[["ncols"]])
    if(is.complex(M[["Mij"]])){
      EigenR_det_sparse_cplx(
        M[["i"]], M[["j"]], M[["Mij"]], M[["nrows"]], M[["ncols"]]
      )
    }else{
      EigenR_det_sparse_real(
        M[["i"]], M[["j"]], M[["Mij"]], M[["nrows"]], M[["ncols"]]
      )
    }
  }else{
    stopifnot(isSquareMatrix(M))
    stopifnot(isRealOrComplex(M))
    if(is.complex(M)){
      EigenR_det_cplx(Re(M), Im(M))
    }else{
      EigenR_det_real(M)
    }
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
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    EigenR_rank_cplx(Re(M), Im(M))
  }else{
    EigenR_rank_real(M)
  }
}

#' Inverse of a matrix
#' 
#' @description Inverse of a real or complex matrix.
#'
#' @param M an invertible square matrix, real or complex
#'
#' @return The inverse matrix of \code{M}.
#' @export
Eigen_inverse <- function(M){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    parts <- EigenR_inverse_cplx(Re(M), Im(M))
    Minv <- parts[["real"]] + 1i * parts[["imag"]]
  }else{
    Minv <- EigenR_inverse_real(M)
  }
  Minv
}

#' Kernel of a matrix
#' 
#' @description Kernel (null-space) of a real or complex matrix.
#'
#' @param M a matrix, real or complex
#' @param method one of \code{"COD"} or \code{"LU"}; the faster method depends 
#'   on the size of the matrix
#'
#' @return A basis of the kernel of \code{M}. With \code{method = "COD"}, the 
#'   basis is orthonormal, while it is not with \code{method = "LU"}.
#' @export
#'
#' @examples set.seed(666)
#' M <- matrix(rgamma(30L, 12, 1), 10L, 3L)
#' M <- cbind(M, M[,1]+M[,2], M[,2]+2*M[,3])
#' # basis of the kernel of `M`:
#' Eigen_kernel(M, method = "LU")
#' # orthonormal basis of the kernel of `M`:
#' Eigen_kernel(M, method = "COD")
Eigen_kernel <- function(M, method = "COD"){
  stopifnot(is.matrix(M))
  stopifnot(isRealOrComplex(M))
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
#' 
#' @description Range (column-space, image, span) of a real or complex matrix.
#'
#' @param M a matrix, real or complex
#' @param method one of \code{"LU"}, \code{"QR"}, or \code{"COD"}; the 
#'   \code{"LU"} method is faster
#'
#' @return A basis of the range of \code{M}. With \code{method = "LU"}, the 
#'   basis is not orthonormal, while it is with \code{method = "QR"} and 
#'   \code{method = "COD"}.
#' @export
Eigen_range <- function(M, method = "QR"){
  stopifnot(is.matrix(M))
  stopifnot(isRealOrComplex(M))
  method <- match.arg(method, c("LU", "QR", "COD"))
  if(is.complex(M)){
    if(method == "QR"){
      parts <- EigenR_image_QR_cplx(Re(M), Im(M))
    }else if(method == "LU"){
      parts <- EigenR_image_LU_cplx(Re(M), Im(M))
    }else{
      parts <- EigenR_image_COD_cplx(Re(M), Im(M))
    }
    parts[["real"]] + 1i * parts[["imag"]]
  }else{
    if(method == "QR"){
      EigenR_image_QR_real(M)
    }else if(method == "LU"){
      EigenR_image_LU_real(M)
    }else{
      EigenR_image_COD_real(M)
    }
  }
}

#' QR decomposition of a matrix
#' 
#' @description QR decomposition of a real or complex matrix.
#'
#' @param M a matrix, real or complex
#'
#' @return A list with the \code{Q} matrix and the \code{R} matrix.
#' @export
#'
#' @examples M <- cbind(c(1,2,3), c(4,5,6))
#' x <- Eigen_QR(M)
#' x$Q %*% x$R
Eigen_QR <- function(M){
  stopifnot(is.matrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    QRparts <- EigenR_QR_cplx(Re(M), Im(M))
    lapply(QRparts, function(parts) parts[["real"]] + 1i * parts[["imag"]])
  }else{
    EigenR_QR_real(M)
  }
}

#' Cholesky decomposition of a matrix
#' 
#' @description Cholesky decomposition of a symmetric or Hermitian matrix.
#'
#' @param M a square symmetric/Hermitian positive-definite matrix or 
#'   \code{\link{SparseMatrix}}, real/complex
#'
#' @return The upper triangular factor of the Cholesky decomposition of 
#'   \code{M}.
#' @export
#' 
#' @details Symmetry is not checked; only the lower triangular part of 
#'   \code{M} is used.
#'
#' @examples M <- rbind(c(5,1), c(1,3))
#' U <- Eigen_chol(M)
#' t(U) %*% U # this is `M`
#' # a Hermitian example:
#' A <- rbind(c(1,1i), c(1i,2))
#' ( M <- A %*% t(Conj(A)) )
#' try(chol(M)) # fails
#' U <- Eigen_chol(M) 
#' t(Conj(U)) %*% U # this is `M`
#' # a sparse example
#' M <- asSparseMatrix(diag(1:5))
#' Eigen_chol(M)
Eigen_chol <- function(M){
  if(inherits(M, "SparseMatrix")){
    stopifnot(M[["nrows"]] == M[["ncols"]])
    if(is.complex(M[["Mij"]])){
      EigenR_chol_sparse_cplx(
        M[["i"]], M[["j"]], M[["Mij"]], M[["nrows"]], M[["ncols"]]
      )
    }else{
      EigenR_chol_sparse_real(
        M[["i"]], M[["j"]], M[["Mij"]], M[["nrows"]], M[["ncols"]]
      )
    }
  }else{
    stopifnot(isSquareMatrix(M))
    stopifnot(isRealOrComplex(M))
    if(is.complex(M)){
      EigenR_chol_cplx(Re(M), Im(M))
    }else{
      EigenR_chol_real(M)
    }
  }
}

#' 'UtDU' decomposition of a matrix
#' 
#' @description Cholesky-'UtDU' decomposition of a symmetric or Hermitian matrix.
#'
#' @param M a square symmetric/Hermitian positive or negative semidefinite 
#'   matrix, real/complex
#'
#' @return The Cholesky-'UtDU' decomposition of \code{M} in a list 
#'   (see example).
#' @export
#' 
#' @details Symmetry is not checked; only the lower triangular part of 
#'   \code{M} is used.
#'
#' @examples x <- matrix(c(1:5, (1:5)^2), 5, 2)
#' x <- cbind(x, x[, 1] + 3*x[, 2])
#' M <- crossprod(x)
#' UtDU <- Eigen_UtDU(M)
#' U <- UtDU$U
#' D <- UtDU$D
#' perm <- UtDU$perm
#' UP <- U[, perm]
#' t(UP) %*% diag(D) %*% UP # this is `M`
Eigen_UtDU <- function(M){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    utdu <- EigenR_UtDU_cplx(Re(M), Im(M))
    utdu[["U"]] <- utdu[["U"]][["real"]] + 1i * utdu[["U"]][["imag"]]
    utdu[["D"]] <- utdu[["D"]][["real"]] + 1i * utdu[["D"]][["imag"]]
    utdu
  }else{
    EigenR_UtDU_real(M)
  }
}

#' Linear least-squares problems
#' 
#' @description Solves a linear least-squares problem.
#'
#' @param A a \code{n*p} matrix, real or complex
#' @param b a vector of length \code{n} or a matrix with \code{n} rows, 
#'   real or complex
#'
#' @return The solution \code{X} of the least-squares problem \code{AX ~= b} 
#'   (similar to \code{lm.fit(A, b)$coefficients}). This is a matrix if 
#'   \code{b} is a matrix, or a vector if \code{b} is a vector.
#' @export
#'
#' @examples set.seed(129)
#' n <- 7; p <- 2
#' A <- matrix(rnorm(n * p), n, p)
#' b <- rnorm(n)
#' lsfit <- Eigen_lsSolve(A, b)
#' b - A %*% lsfit # residuals
Eigen_lsSolve <- function(A, b){
  stopifnot(is.matrix(A)) 
  stopifnot(is.atomic(b))
  b <- cbind(b)
  stopifnot(nrow(A) == nrow(b))
  stopifnot(is.numeric(A) || is.complex(A))
  stopifnot(is.numeric(b) || is.complex(b))
  if(is.complex(A) || is.complex(b)){
    parts <- EigenR_lsSolve_cplx(Re(A), Im(A), Re(b), Im(b))
    parts[["real"]][,] + 1i * parts[["imag"]][,]
  }else{
    EigenR_lsSolve_real(A, b)[,]
  }
}

#' Exponential of a matrix
#' 
#' @description Exponential of a real or complex square matrix.
#'
#' @param M a square matrix, real or complex
#'
#' @return The exponential of \code{M}.
#' @export
Eigen_exp <- function(M){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    parts <- EigenR_exp_cplx(Re(M), Im(M))
    Mexp <- parts[["real"]] + 1i * parts[["imag"]]
  }else{
    Mexp <- EigenR_exp_real(M)
  }
  Mexp
}

#' Logarithm of a matrix
#' 
#' @description Logarithm of a real or complex square matrix, when possible.
#'
#' @param M a square matrix, real or complex
#'
#' @return The logarithm of \code{M}.
#' @details The logarithm of a matrix does not always exist. 
#'   See \href{https://eigen.tuxfamily.org/dox/unsupported/group__MatrixFunctions__Module.html#title7}{matrix logarithm}.
#' @export
Eigen_log <- function(M){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    parts <- EigenR_log_cplx(Re(M), Im(M))
    Mlog <- parts[["real"]] + 1i * parts[["imag"]]
  }else{
    Mlog <- EigenR_log_real(M)
  }
  Mlog
}

#' Matrix cosine
#' 
#' @description Matrix cosine of a real or complex square matrix.
#'
#' @param M a square matrix, real or complex
#'
#' @return The matrix cosine of \code{M}.
#' @export
#' @examples library(EigenR)
#' M <- toeplitz(c(1,2,3))
#' cosM <- Eigen_cos(M) 
#' sinM <- Eigen_sin(M)
#' cosM %*% cosM + sinM %*% sinM # identity matrix
Eigen_cos <- function(M){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    parts <- EigenR_cos_cplx(Re(M), Im(M))
    Mcos <- parts[["real"]] + 1i * parts[["imag"]]
  }else{
    Mcos <- EigenR_cos_real(M)
  }
  Mcos
}

#' Matrix sine
#' 
#' @description Matrix sine of a real or complex square matrix.
#'
#' @param M a square matrix, real or complex
#'
#' @return The matrix sine of \code{M}.
#' @export
Eigen_sin <- function(M){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    parts <- EigenR_sin_cplx(Re(M), Im(M))
    Msin <- parts[["real"]] + 1i * parts[["imag"]]
  }else{
    Msin <- EigenR_sin_real(M)
  }
  Msin
}

#' Matrix hyperbolic cosine
#' 
#' @description Matrix hyperbolic cosine of a real or complex square matrix.
#'
#' @param M a square matrix, real or complex
#'
#' @return The matrix hyperbolic cosine of \code{M}.
#' @export
#' @examples library(EigenR)
#' M <- toeplitz(c(1,2,3))
#' Eigen_cosh(M)
#' (Eigen_exp(M) + Eigen_exp(-M)) / 2 # identical
Eigen_cosh <- function(M){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    parts <- EigenR_cosh_cplx(Re(M), Im(M))
    Mcosh <- parts[["real"]] + 1i * parts[["imag"]]
  }else{
    Mcosh <- EigenR_cosh_real(M)
  }
  Mcosh
}

#' Matrix hyperbolic sine
#' 
#' @description Matrix hyperbolic sine of a real or complex square matrix.
#'
#' @param M a square matrix, real or complex
#'
#' @return The matrix hyperbolic sine of \code{M}.
#' @export
#' @examples library(EigenR)
#' M <- toeplitz(c(1,2,3))
#' Eigen_sinh(M)
#' (Eigen_exp(M) - Eigen_exp(-M)) / 2  # identical
Eigen_sinh <- function(M){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    parts <- EigenR_sinh_cplx(Re(M), Im(M))
    Msinh <- parts[["real"]] + 1i * parts[["imag"]]
  }else{
    Msinh <- EigenR_sinh_real(M)
  }
  Msinh
}

#' Matricial power
#' 
#' @description Matricial power of a real or complex square matrix, when possible.
#'
#' @param M a square matrix, real or complex
#' @param p a number, real or complex, the power exponent
#'
#' @return The matrix \code{M} raised at the power \code{p}.
#' @details The power is defined with the help of the exponential and the 
#'   logarithm. See \href{https://eigen.tuxfamily.org/dox/unsupported/group__MatrixFunctions__Module.html#title8}{matrix power}.
#' @export
Eigen_pow <- function(M, p){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  stopifnot(isRealOrComplexScalar(p))
  if(is.complex(p) && !is.complex(M)){
    Zeros <- matrix(0, nrow = nrow(M), ncol = ncol(M))
    parts <- EigenR_pow_cplx(Re(M), Zeros, p)
    Mpow <- parts[["real"]] + 1i * parts[["imag"]]
  }else if(is.complex(M) && !is.complex(p)){
    p <- as.complex(p)
  }
  if(is.complex(M) && is.complex(p)){
    parts <- EigenR_pow_cplx(Re(M), Im(M), p)
    Mpow <- parts[["real"]] + 1i * parts[["imag"]]
  }else if(!is.complex(M) && !is.complex(p)){
    Mpow <- EigenR_pow_real(M, p)
  }
  Mpow
}


#' Square root of a matrix
#' 
#' @description Square root of a real or complex square matrix, when possible.
#'
#' @param M a square matrix, real or complex
#'
#' @return A square root of \code{M}.
#' @details See \href{https://eigen.tuxfamily.org/dox/unsupported/group__MatrixFunctions__Module.html#title12}{matrix square root}.
#' @export
#' @examples # Rotation matrix over 60 degrees:
#' M <- cbind(c(cos(pi/3), sin(pi/3)), c(-sin(pi/3), cos(pi/3)))
#' # Its square root, the rotation matrix over 30 degrees:
#' Eigen_sqrt(M)
Eigen_sqrt <- function(M){
  stopifnot(isSquareMatrix(M))
  stopifnot(isRealOrComplex(M))
  if(is.complex(M)){
    parts <- EigenR_sqrt_cplx(Re(M), Im(M))
    Msqrt <- parts[["real"]] + 1i * parts[["imag"]]
  }else{
    Msqrt <- EigenR_sqrt_real(M)
  }
  Msqrt
}
