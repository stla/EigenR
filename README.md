EigenR
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/stla/EigenR/workflows/R-CMD-check/badge.svg)](https://github.com/stla/EigenR/actions)
<!-- badges: end -->

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
##           expr    min      lq     mean  median      uq     max neval
##           base 7.6210 8.06630 9.627830 8.63330 9.77255 30.7122   200
##         EigenR 2.7194 3.74275 4.246512 3.96390 4.44920  6.7440   200
##  EigenR_sparse 5.8260 6.28425 7.294666 6.64745 7.84145 14.1825   200
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
##         expr     min     lq      mean  median      uq     max neval
##       EigenR  1.4003  1.494  1.735253  1.5589  1.8048  2.9544    30
##  complexplus 16.4352 17.626 19.880123 18.0317 20.2785 36.2521    30
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
## Unit: milliseconds
##    expr    min      lq     mean median     uq    max neval
##    base 1.2242 1.38025 1.839062 1.6311 2.2023 7.1696   500
##  EigenR 1.0135 1.12585 1.584468 1.3026 2.0330 3.5632   500
```

## Pseudo-inverse matrix

``` r
set.seed(666L)
M <- matrix(rnorm(100L*70L), 100L, 70L)
library(MASS)
library(pracma)
microbenchmark(
  MASS   = ginv(M),
  pracma = pinv(M),
  EigenR = Eigen_pinverse(M), # :-)
  times = 500L
)
## Unit: microseconds
##    expr    min      lq     mean  median      uq     max neval
##    MASS 3175.2 3414.45 4297.759 3987.90 5218.75  8736.1   500
##  pracma 3124.0 3403.20 4245.509 3788.70 5002.15 11380.4   500
##  EigenR  693.1  808.55 1043.352  897.35 1328.60  4042.1   500
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
##    expr   min     lq     mean median     uq    max neval
##    base 173.0 178.15 217.1493  185.0 219.20 5069.7  1000
##  EigenR 131.1 136.70 170.6095  144.3 174.95 3631.7  1000
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
##    expr    min     lq     mean median      uq     max neval
##    base 1423.6 1600.3 1986.700 1739.0 2186.95 10603.7  1000
##  EigenR  692.7 1232.2 1480.108 1397.2 1639.30  8341.2  1000
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
##        expr    min      lq     mean median      uq     max neval
##        MASS 7.8387 8.51355 9.382178 9.0692 9.72050 17.0432   100
##   EigenR_LU 1.8250 2.19795 2.524910 2.3874 2.69615  4.6348   100
##  EigenR_COD 4.8351 5.21975 6.027508 5.6400 6.25525 17.4528   100
```

## Linear least-squares problems

``` r
set.seed(666L)
n <- 700L; p <- 200L
A <- matrix(rnorm(n * p), n, p)
b <- rnorm(n)
microbenchmark(
  stats       = lm.fit(A, b),
  EigenR_svd  = Eigen_lsSolve(A, b, method = "svd"), #  :-(
  EigenR_cod  = Eigen_lsSolve(A, b, method = "cod"), #  :-)
  times = 100L
)
## Unit: milliseconds
##        expr     min       lq     mean   median       uq     max neval
##       stats 26.7695 27.86835 29.11444 28.59330 29.70460 35.8353   100
##  EigenR_svd 50.5333 52.18155 55.78505 54.01155 58.15850 72.2574   100
##  EigenR_cod 12.6617 13.44090 14.70645 13.97585 15.03515 23.6157   100
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
##    expr   min      lq      mean  median      uq    max neval
##    expm 998.5 1070.65 1227.2344 1110.35 1215.40 9754.9   500
##  EigenR 216.7  232.15  271.8416  254.45  285.95  700.6   500
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
##         expr    min      lq     mean  median     uq      max neval
##       EigenR 1.2714 1.39295 1.618107 1.45955 1.6503   4.6186   500
##  complexplus 6.4594 6.88835 7.941399 7.15010 8.1093 113.3130   500
```
