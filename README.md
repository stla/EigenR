EigenR
================

Originally, I entitled this package *Fast Matrix Algebra with ‘Eigen’*,
because I expected it to be faster than R base. But this is not the
case. So I entitled it *Complex Matrix Algebra with ‘Eigen’*, because it
supports some operations on complex matrices which are not supported by
R base: determinant, Cholesky decomposition, and linear least-squares
problems.

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
##           expr      min        lq      mean    median        uq       max neval
##           base 2.348464  4.842405 31.934640  8.706389 29.089927 249.68330   200
##         EigenR 3.587675  6.706318  8.785941  7.363495  9.025337  26.93790   200
##  EigenR_sparse 9.570377 12.717644 15.615663 13.922106 15.723067  42.66454   200
##  cld
##    b
##   a 
##   a
```

Determinants of complex matrices are supported:

``` r
set.seed(666L)
Mr <- matrix(rnorm(300L*300L, mean = 1), 300L, 300L)
Mi <- matrix(rnorm(300L*300L, mean = 1), 300L, 300L)
M <- Mr + 1i * Mi
library(complexplus)
microbenchmark(
  EigenR      = Eigen_det(M), # :-)
  complexplus = Det(M), 
  times = 10L
)
## Unit: milliseconds
##         expr       min        lq      mean    median        uq        max neval
##       EigenR  29.13314  45.90238  57.99029  54.97324  72.55315   87.81947    10
##  complexplus 616.27011 682.67422 850.77522 848.73617 997.65058 1130.64936    10
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
  EigenR = Eigen_chol(A), # :-|
  times = 500L
)
## Unit: microseconds
##    expr     min      lq     mean   median       uq       max neval cld
##    base 199.183 261.304 889.3098 322.1920 654.8255 33523.410   500   b
##  EigenR 162.540 367.595 442.7996 388.7965 475.6985  3355.604   500  a
```

Cholesky decomposition of complex matrices is supported.

## Pivoted Cholesky decomposition

``` r
set.seed(666L)
M <- matrix(rgamma(202L*199L, 10), 202L, 199L)
M <- cbind(M, M[, 1L] + 3*M[, 2L])
A <- crossprod(M)
microbenchmark(
  base   = chol(A, pivot = TRUE),
  EigenR = Eigen_UtDU(A), # :-|
  times = 500L
)
## Unit: milliseconds
##    expr      min       lq     mean   median       uq      max neval cld
##    base 1.107769 1.326954 3.515281 1.782150 3.592810 38.48691   500   b
##  EigenR 1.426556 1.952104 2.608765 2.119216 2.454293 25.39490   500  a
```

Pivoted Cholesky decomposition of complex matrices is supported.

## Kernel

``` r
set.seed(666L)
M <- matrix(rpois(10000L, 25), 100L, 100L)
A <- cbind(M, M)
At <- t(A)
library(MASS)
microbenchmark(
  MASS       = Null(At),
  EigenR_LU  = Eigen_kernel(A, method = "LU"),  # :-)
  EigenR_COD = Eigen_kernel(A, method = "COD"), 
  times = 100L
)
## Unit: milliseconds
##        expr      min       lq     mean   median       uq        max neval cld
##        MASS 5.098672 5.392457 7.114541 5.669072 6.084231 126.900039   100   b
##   EigenR_LU 2.385806 2.511427 2.677017 2.588909 2.692103   5.659926   100  a 
##  EigenR_COD 5.642747 5.977678 6.459983 6.233274 6.695495   9.042926   100   b
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
##    stats 15.82585 15.96568 16.84272 16.49552 17.65445 18.53197    10  a 
##  Eigen_R 82.10392 82.64413 84.44043 84.15784 86.39429 87.31093    10   b
```

Complex matrices `A` and `b` are supported.
