# perm <- function(transpositions){
#   x <- seq_along(transpositions)
#   for(i in x){
#     x[c(i,transpositions[i])] <- x[c(transpositions[i],i)]
#   }
#   x
# }

isStrictPositiveInteger <- function(x){
  all(x > 0 & (floor(x) == x))
}

isSquareMatrix <- function(M){
  is.matrix(M) && (nrow(M) == ncol(M))
}

isRealOrComplex <- function(M){
  is.numeric(M) || is.complex(M)
}

isRealOrComplexScalar <- function(x){
  (is.numeric(x) || is.complex(x)) && length(x) == 1L && !is.na(x)
}
