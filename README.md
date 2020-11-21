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
##           expr      min        lq      mean    median        uq      max neval
##           base 2.295802  3.216064  9.310243  5.653305 11.682343 78.60175   200
##         EigenR 3.254337  6.265068  7.024137  6.708800  7.247874 18.33710   200
##  EigenR_sparse 8.546355 11.727958 13.677009 12.496962 13.192315 36.37084   200
##  cld
##   b 
##  a  
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
  EigenR      = Eigen_det(M), # :-)
  complexplus = Det(M), 
  times = 10L
)
## Unit: milliseconds
##         expr       min        lq      mean    median        uq       max neval
##       EigenR  24.46205  52.61349  53.81241  55.55755  59.23399  65.54072    10
##  complexplus 639.46520 801.67035 835.61555 853.51478 882.21947 956.18605    10
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
##    expr     min       lq     mean   median       uq       max neval cld
##    base 149.803 168.2005 666.4576 191.9505 225.3335 24211.285   500   b
##  EigenR 150.218 315.1445 330.0574 331.0205 347.1055   621.644   500  a
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
## Unit: microseconds
##    expr      min       lq     mean   median       uq      max neval cld
##    base  800.678 1269.482 3267.925 1540.066 4094.699 19141.01   500   b
##  EigenR 1005.404 1910.888 2486.643 2021.510 2530.508 13963.72   500  a
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
  EigenR_LU  = Eigen_kernel(A, method = "LU"),
  EigenR_COD = Eigen_kernel(A, method = "COD"), # :-(
  times = 100L
)
## Unit: milliseconds
##        expr      min       lq     mean   median       uq        max neval cld
##        MASS 4.625058 5.109697 6.704468 5.280767 5.767802 124.878257   100   b
##   EigenR_LU 2.192517 2.373585 2.570679 2.479659 2.700899   5.248459   100  a 
##  EigenR_COD 5.286845 5.591827 6.037580 5.786849 6.332180   8.450543   100   b
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
##    stats 15.26072 15.65720 16.13600 15.98812 16.66596 17.50740    10  a 
##  Eigen_R 78.71971 80.71122 81.22111 81.37661 82.11088 82.68067    10   b
```

Complex matrices `A` and `b` are supported.
