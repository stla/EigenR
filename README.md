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
##           base 2.280564  2.384481 13.997633  2.773429  3.810345 316.68765   200
##         EigenR 3.106442  6.631640  7.502237  7.025754  7.595042  19.05179   200
##  EigenR_sparse 6.250881 12.635205 14.450814 13.095196 14.612821  28.75310   200
##  cld
##   ab
##   a 
##    b
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
##         expr       min        lq     mean    median        uq        max neval
##       EigenR  2.657419  2.850238  3.37646  3.218996  3.481693   5.176787    30
##  complexplus 36.334199 38.599482 48.57074 39.790650 44.854627 262.863927    30
##  cld
##   a 
##    b
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
##    expr      min       lq     mean   median        uq      max neval cld
##    base  624.262  656.725 1001.600  697.907  813.7565 7161.314   500  a 
##  EigenR 1641.917 2448.156 2468.013 2471.816 2506.2770 2758.072   500   b
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
##    expr     min       lq     mean   median       uq      max neval cld
##    base 142.948 155.6255 288.3603 162.2095 204.4600 7900.534  1000  a 
##  EigenR 147.082 320.7905 398.7887 341.4415 395.6865 9949.041  1000   b
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
## Unit: milliseconds
##    expr      min       lq     mean   median       uq       max neval cld
##    base 1.092026 1.187069 1.582198 1.268445 1.425253 12.261755  1000  a 
##  EigenR 1.287782 1.885200 2.080151 1.942995 2.086418  9.553564  1000   b
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
##        MASS 4.567780 4.796990 7.320759 4.993532 5.319909 115.105350   100   b
##   EigenR_LU 2.244941 2.327626 2.481592 2.397051 2.530926   5.318979   100  a 
##  EigenR_COD 5.011778 5.170457 5.421014 5.327620 5.593926   6.863361   100  ab
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
##     expr      min       lq     mean   median       uq      max neval cld
##    stats 14.26740 14.47446 14.88505 14.74442 15.02579 16.73055    20  a 
##  Eigen_R 73.03921 73.84556 74.57334 74.40150 75.18944 76.46069    20   b
```

Complex matrices `A` and `b` are supported.
