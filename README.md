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
##           expr    min      lq      mean  median       uq     max neval
##           base 7.7281 8.67860 11.096060 10.5468 13.28665 25.6797   200
##         EigenR 3.3535 3.90675  4.933978  4.8636  5.73430  8.0125   200
##  EigenR_sparse 6.0182 6.92185  8.341857  8.2838  9.69950 12.9405   200
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
##         expr     min      lq     mean   median      uq     max neval
##       EigenR  1.3335  1.4878  1.87118  1.72875  2.2763  2.8158    30
##  complexplus 17.0684 19.4424 21.05672 21.03340 22.3398 30.0007    30
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
##    expr    min     lq     mean  median     uq    max neval
##    base 1221.1 1327.8 1532.809 1410.25 1624.8 4418.2   500
##  EigenR  987.4 1065.1 1224.727 1119.85 1263.0 2621.8   500
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
##    expr    min      lq      mean  median      uq     max neval
##    MASS 3175.2 3414.70 3977.1486 3666.55 4217.90 14669.1   500
##  pracma 3173.6 3380.85 4002.3044 3710.35 4282.25  9449.8   500
##  EigenR  679.3  793.65  951.0248  838.10 1005.55  9904.9   500
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
##    expr   min    lq     mean median    uq     max neval
##    base 173.6 216.0 275.4291  230.7 272.4 11841.5  1000
##  EigenR 121.9 164.8 211.8479  178.0 217.1 11746.6  1000
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
##    expr    min     lq     mean  median     uq     max neval
##    base 1518.0 1671.9 2008.952 1764.35 2071.6 15535.5  1000
##  EigenR  765.8 1203.1 1498.988 1403.10 1599.5 12758.5  1000
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
##        expr    min      lq     mean  median     uq     max neval
##        MASS 7.8542 8.38020 9.410370 9.15620 9.8947 16.4025   100
##   EigenR_LU 1.7730 2.01270 2.475684 2.28625 2.6359  7.6111   100
##  EigenR_COD 4.7002 5.18835 5.977337 5.49445 6.1066 18.9573   100
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
##     expr     min      lq     mean   median       uq     max neval
##    stats 27.6842 28.2868 29.78608 28.89755 31.57220 34.6231    20
##  Eigen_R 52.2499 53.8344 58.25450 54.93255 61.64995 76.3218    20
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
##    expr    min      lq      mean  median      uq    max neval
##    expm 1009.1 1090.95 1452.7794 1234.05 1752.85 7207.2   500
##  EigenR  221.0  247.10  343.4398  283.80  455.40 1022.9   500
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
##         expr    min     lq     mean  median     uq      max neval
##       EigenR 1.2758 1.4008 1.639416 1.47615 1.7298   5.5131   500
##  complexplus 6.4326 6.9897 8.076339 7.36210 8.2665 124.0744   500
```
