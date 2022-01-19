library(EigenR)
library(microbenchmark)

set.seed(666L)
n <- 700L; p <- 200L
A <- matrix(rnorm(n * p), n, p)
b <- rnorm(n)
microbenchmark(
  stats       = lm.fit(A, b),
  EigenR_svd  = Eigen_lsSolve(A, b), 
  EigenR_cod  = Eigen_lsSolve(A, b, method = "cod"), 
  times = 200L
)


