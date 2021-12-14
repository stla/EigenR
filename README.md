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
##           expr      min       lq      mean    median      uq     max neval
##           base 7.678000 9.133401 11.930930 11.023601 13.8092 57.2875   200
##         EigenR 3.385101 4.051051  5.470443  5.264401  6.2005 27.8239   200
##  EigenR_sparse 5.841101 6.773851  8.705014  8.288101 10.1530 21.0749   200
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
##         expr       min        lq      mean   median        uq       max neval
##       EigenR  1.500201  1.698701  2.091281  2.16685  2.422002  2.614001    30
##  complexplus 16.929602 17.893801 20.250934 19.15265 21.995201 26.774302    30
```

## Inverse matrix

``` r
set.seed(666L)
M <- matrix(rnorm(100L*100L), 100L, 100L)
microbenchmark(
  base   = solve(M),
  EigenR = Eigen_inverse(M), 
  times = 500L
)
## Unit: microseconds
##    expr      min       lq     mean   median       uq      max neval
##    base 1256.502 1367.601 1655.012 1508.351 1828.051 6815.901   500
##  EigenR  977.302 1089.050 1378.164 1195.901 1563.450 3187.401   500
```

## Cholesky decomposition

``` r
set.seed(666L)
M <- matrix(rpois(10000L, 25), 100L, 100L)
A <- crossprod(M)
microbenchmark(
  base   = chol(A),
  EigenR = Eigen_chol(A), 
  times = 1000L
)
## Unit: microseconds
##    expr     min       lq     mean   median      uq      max neval
##    base 173.101 191.5515 257.3574 224.1005 297.651 4097.900  1000
##  EigenR 119.802 135.5010 209.8055 171.2510 251.251 3122.201  1000
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
##    expr      min        lq     mean   median       uq      max neval
##    base 1436.701 1606.3010 1901.581 1763.101 2021.552 8589.702  1000
##  EigenR  644.301  926.9015 1320.769 1266.552 1513.551 6806.600  1000
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
##        expr      min       lq     mean   median        uq     max neval
##        MASS 7.910201 8.544351 9.530592 9.133401 10.167551 17.1191   100
##   EigenR_LU 1.779500 2.155702 2.663815 2.447001  2.946551  5.0123   100
##  EigenR_COD 4.861101 5.432752 6.466603 5.972050  7.240001 20.5207   100
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
##    stats 27.6975 28.28305 30.30363 29.18565 33.06465 36.5224    20
##  Eigen_R 52.6666 53.52400 58.19708 55.06110 60.92345 84.8454    20
```

Complex matrices `A` and `b` are supported.

## Exponential

``` r
set.seed(666L)
M <- matrix(rnorm(40L*40L, mean = 1), 40L, 40L)
microbenchmark(
  expm   = expm::expm(M),
  EigenR = Eigen_exp(M), # :-)
  times = 500L
)
## Unit: microseconds
##    expr      min       lq      mean    median       uq      max neval
##    expm 1002.701 1071.901 1330.1628 1143.9000 1427.502 5827.901   500
##  EigenR  213.101  236.601  302.2438  265.3015  329.351 1399.101   500
```

Exponential of complex matrices is supported:

``` r
set.seed(666L)
Mr <- matrix(rnorm(40L*40L, mean = 1), 40L, 40L)
Mi <- matrix(rnorm(40L*40L, mean = 1), 40L, 40L)
M <- Mr + 1i * Mi
library(complexplus)
microbenchmark(
  EigenR      = Eigen_exp(M), # :-)
  complexplus = matexp(M), 
  times = 500L
)
## Unit: milliseconds
##         expr      min       lq     mean   median       uq       max neval
##       EigenR 1.281301 1.419850 1.775880 1.606201 2.043300  3.442001   500
##  complexplus 6.515601 7.076051 8.320601 7.755651 9.039451 16.800501   500
```
