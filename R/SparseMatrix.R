#' Sparse matrix
#' 
#' @description Constructs a sparse matrix, real or complex.
#'
#' @param i,j indices of the non-zero coefficients
#' @param Mij values of the non-zero coefficients; must be a vector of the same 
#'   length as \code{i} and \code{j} or a single number which will be recycled
#' @param nrows,ncols dimensions of the matrix
#' @param x a \code{SparseMatrix} object
#' @param ... ignored
#'
#' @return A list with the class \code{SparseMatrix}.
#' @export
#' 
#' @name SparseMatrix
#'
#' @examples xx
SparseMatrix <- function(i, j, Mij, nrows, ncols){
  stopifnot(is.atomic(i), is.atomic(j), is.atomic(Mij))
  stopifnot(isStrictPositiveInteger(i), isStrictPositiveInteger(j))
  stopifnot(is.numeric(Mij) || is.complex(Mij))
  stopifnot(length(nrows) == 1L, length(ncols) == 1L)
  stopifnot(isStrictPositiveInteger(nrows), isStrictPositiveInteger(ncols))
  stopifnot(max(i) <= nrows, max(j) <= ncols)
  stopifnot(length(i) == length(j))
  stopifnot(length(Mij) == 1L || length(Mij) == length(i))
  if(length(Mij) == 1L){
    Mij <- rep(Mij, times = length(i))
  }
  out <- list(i = i, j = j, Mij = Mij, nrows = nrows, ncols = ncols)
  class(out) <- "SparseMatrix"
  out
}

#' @rdname SparseMatrix
#' @export
print.SparseMatrix <- function(x, ...){
  M <- matrix(".", nrow = x$nrows, ncol = x$ncols)
  Mij <- formatC(x$Mij)
  M[cbind(x$i,x$j)] <- Mij
  print(M, quote = FALSE)
}
