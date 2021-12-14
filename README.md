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
##           expr      min       lq      mean    median       uq     max neval
##           base 8.008401 9.439550 12.360928 11.772001 14.25095 33.3621   200
##         EigenR 2.948801 4.360351  5.575071  5.634151  6.35175 10.9766   200
##  EigenR_sparse 6.001900 7.378351  9.051891  8.702951 10.33490 20.7182   200
```

Determinants of complex matrices are supported:

``` r
set.seed(666L)
Mr <- matrix(rnorm(100L*100L, mean = 1), 100L, 100L)
Mi <- matrix(rnorm(100L*100L, mean = 1), 100L, 100L)
M <- Mr + 1i * Mi
library(complexplus)
microbenchmark(
  EigenR      = Eigen_det(M), # :-)
  complexplus = Det(M), 
  times = 30L
)
## Unit: milliseconds
##         expr       min        lq      mean    median        uq       max neval
##       EigenR  1.441501  1.547502  1.928104  1.733801  2.076301  4.014701    30
##  complexplus 16.981601 18.186601 20.988558 18.945551 21.245501 33.804401    30
```

## Inverse matrix

``` r
set.seed(666L)
M <- matrix(rnorm(100L*100L), 100L, 100L)
microbenchmark(
  base   = solve(M),
  EigenR = Eigen_inverse(M), # :-(
  times = 500L
)
## Unit: microseconds
##    expr      min       lq     mean   median       uq       max neval
##    base 1262.400 1393.401 1696.938 1543.051 1857.151 11106.500   500
##  EigenR  976.101 1099.951 1408.154 1212.501 1680.851  3836.901   500
```

## Cholesky decomposition

``` r
set.seed(666L)
M <- matrix(rpois(10000L, 25), 100L, 100L)
A <- crossprod(M)
microbenchmark(
  base   = chol(A),
  EigenR = Eigen_chol(A), # :-(
  times = 1000L
)
## Unit: microseconds
##    expr     min       lq     mean   median      uq      max neval
##    base 173.501 185.1010 246.3900 215.9515 291.451  782.501  1000
##  EigenR 118.702 131.7515 212.9907 167.4505 245.001 5333.601  1000
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
  EigenR = Eigen_UtDU(A), # :-(
  times = 1000L
)
## Unit: microseconds
##    expr      min       lq     mean   median       uq     max neval
##    base 1441.801 1658.901 2139.929 1940.301 2394.202 11157.9  1000
##  EigenR  653.102 1014.701 1506.732 1388.051 1820.851 14772.8  1000
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
##        expr      min       lq     mean   median        uq       max neval
##        MASS 7.888201 8.646801 9.705527 9.257001 10.013302 19.889201   100
##   EigenR_LU 1.810100 2.235001 2.857512 2.587801  3.381900  5.337701   100
##  EigenR_COD 4.645302 5.082600 6.320570 5.709651  7.366052 16.115100   100
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
  times = 20L
)
## Unit: milliseconds
##     expr     min       lq     mean   median       uq     max neval
##    stats 31.1994 31.68015 34.07953 33.11745 35.75225 41.8597    20
##  Eigen_R 64.6860 69.21850 71.62660 70.56360 74.16900 80.3629    20
```

Complex matrices `A` and `b` are supported.

## Exponential

``` r
set.seed(666L)
M <- matrix(rnorm(30L*30L, mean = 1), 30L, 30L)
microbenchmark(
  expm   = expm::expm(M),
  EigenR = Eigen_exp(M), # :-)
  times = 500L
)
## Unit: microseconds
##    expr     min       lq      mean    median       uq      max neval
##    expm 555.900 998.7005 1185.3213 1123.3515 1250.501 8260.901   500
##  EigenR 122.001 215.5015  276.7422  263.8005  295.801 4423.201   500
```

Exponential of complex matrices is supported:

``` r
set.seed(666L)
Mr <- matrix(rnorm(30L*30L, mean = 1), 30L, 30L)
Mi <- matrix(rnorm(30L*30L, mean = 1), 30L, 30L)
M <- Mr + 1i * Mi
library(complexplus)
microbenchmark(
  EigenR      = Eigen_exp(M), # :-)
  complexplus = matexp(M), 
  times = 500L
)
## Unit: microseconds
##         expr      min       lq      mean    median       uq       max neval
##       EigenR  575.201  678.850  933.5834  823.7505 1192.001  2036.401   500
##  complexplus 2650.301 3016.001 4003.7130 3660.4015 4802.751 10809.201   500
```
