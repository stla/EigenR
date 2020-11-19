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
##           expr      min        lq      mean    median        uq      max neval
##           base 2.325018  2.606910  3.851104  3.177992  3.879467 13.89208   200
##         EigenR 3.560423  6.757584  7.601636  7.120388  7.747838 17.32185   200
##  EigenR_sparse 6.318647 12.694837 13.865138 13.080720 13.771776 24.00013   200
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
##    expr     min       lq     mean  median       uq      max neval cld
##    base 144.756 155.2175 350.0837 166.133 215.9040 6743.635   500   a
##  EigenR 162.359 324.2405 357.6519 329.014 340.4125 6888.029   500   a
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
##        MASS 4.852772 5.129205 5.508211 5.290449 5.598714 8.693911   100  b 
##   EigenR_LU 2.210166 2.378106 2.604436 2.482603 2.613480 5.466532   100 a  
##  EigenR_COD 5.398932 5.667117 6.074435 5.896511 6.261531 9.272335   100   c
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
##    stats 15.29396 15.90042 16.18114 16.20677 16.48451 16.84350    10  a 
##  Eigen_R 78.01373 78.39762 79.67060 79.55024 80.70331 81.91639    10   b
```
