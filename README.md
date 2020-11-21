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
##           base  2.309594  2.368946  3.643059  2.612380  3.216231 92.44736   200
##         EigenR  2.996998  6.556045  6.931398  6.830802  7.170283 10.52983   200
##  EigenR_sparse 11.876631 12.351218 12.985456 12.556399 13.067759 43.75062   200
##  cld
##  a  
##   b 
##    c
```

Determinants of complex matrices are supported:

``` r
set.seed(666L)
Mr <- matrix(rnorm(300L*300L, mean = 1), 300L, 300L)
Mi <- matrix(rnorm(300L*300L, mean = 1), 300L, 300L)
M <- Mr + 1i * Mi
library(complexplus)
microbenchmark(
  EigenR      = Eigen_det(M),
  complexplus = Det(M), # :-(
  times = 10L
)
## Unit: milliseconds
##         expr      min        lq      mean    median       uq       max neval
##       EigenR  23.7719  50.45722  50.04908  52.33652  55.5562  56.22574    10
##  complexplus 488.7458 508.31172 517.47328 516.76704 528.3779 534.35087    10
##  cld
##   a 
##    b
```

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
##    expr     min       lq     mean  median       uq       max neval cld
##    base 145.787 156.3220 306.2726 161.925 170.1975 13567.840   500   a
##  EigenR 150.031 324.5875 351.3277 328.020 335.7955  4440.927   500   a
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
##        expr      min       lq     mean   median       uq        max neval cld
##        MASS 4.568793 4.694112 6.125222 4.842747 5.071969 114.094655   100   b
##   EigenR_LU 2.152736 2.193173 2.317315 2.244509 2.364692   4.524372   100  a 
##  EigenR_COD 5.043219 5.144412 5.377307 5.254402 5.454072   7.727485   100   b
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
##     expr      min       lq     mean   median       uq      max neval cld
##    stats 14.28079 14.36905 14.78100 14.57292 15.16738 15.75864    10  a 
##  Eigen_R 73.01995 73.34898 75.11627 73.71472 74.60609 85.31473    10   b
```

Complex matrices `A` and `b` are supported.
