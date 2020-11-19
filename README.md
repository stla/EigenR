EigenR
================

``` r
library(EigenR)
library(microbenchmark)
```

## Determinant

``` r
set.seed(666L)
M <- matrix(rnorm(300L*300L, mean = 1), 300L, 300L)
M[sample.int(300L*300L, 300L*270L)] <- 0 # 90% of zeros
Ms <- asSparseMatrix(M)
microbenchmark(
  base          = det(M),
  EigenR        = Eigen_det(M),
  EigenR_sparse = Eigen_det(Ms), # :-(
  times = 200L
)
## Unit: milliseconds
##           expr       min        lq      mean    median        uq      max neval
##           base  2.303390  2.556973  3.691258  2.944695  3.749370 14.04845   200
##         EigenR  3.287136  6.601224  7.150831  7.018299  7.438121 12.15553   200
##  EigenR_sparse 11.120235 12.330039 13.044584 12.715541 13.244792 23.46231   200
##  cld
##  a  
##   b 
##    c
```

Determinants of complex matrices are supported.

## Cholesky decomposition

``` r
set.seed(666L)
M <- matrix(rpois(10000L, 25), 100L, 100L)
A <- crossprod(M)
microbenchmark(
  base   = chol(A),
  EigenR = Eigen_chol(A), # :-(
  times = 500L
)
## Unit: microseconds
##    expr     min       lq     mean   median       uq       max neval cld
##    base 144.696 155.6465 362.2259 166.5185 201.3765 12887.677   500   a
##  EigenR 161.502 324.0675 354.6938 328.8125 338.6975  9421.688   500   a
```

Cholesky decomposition of complex matrices is supported.

## Kernel

``` r
set.seed(666L)
M <- matrix(rpois(10000L, 25), 100L, 100L)
A <- cbind(M, M)
At <- t(A)
library(MASS)
microbenchmark(
  MASS       = Null(At),
  EigenR_LU  = Eigen_kernel(A, method = "LU"),
  EigenR_COD = Eigen_kernel(A, method = "COD"), # :-(
  times = 100L
)
## Unit: milliseconds
##        expr      min       lq     mean   median       uq      max neval cld
##        MASS 4.800403 5.068113 5.407869 5.191852 5.528975 9.353463   100  b 
##   EigenR_LU 2.233529 2.357797 2.568849 2.443474 2.599349 5.727903   100 a  
##  EigenR_COD 5.302908 5.559947 5.924981 5.787556 6.021780 9.626854   100   c
```

## Linear least-squares problems

``` r
set.seed(666L)
n <- 700L; p <- 200L
A <- matrix(rnorm(n * p), n, p)
b <- rnorm(n)
microbenchmark(
  stats   = lm.fit(A, b),
  Eigen_R = Eigen_lsSolve(A, b), # :-(
  times = 10L
)
## Unit: milliseconds
##     expr      min       lq     mean   median      uq      max neval cld
##    stats 15.05648 15.17506 15.63534 15.50915 15.7935 16.62512    10  a 
##  Eigen_R 77.36481 78.50575 79.28756 78.75564 80.9617 81.92589    10   b
```

Complex matrices `A` and `b` are supported.
