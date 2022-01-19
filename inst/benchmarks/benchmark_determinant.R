library(EigenR)
library(microbenchmark)
set.seed(666L)
M <- matrix(rnorm(300L*300L, mean = 1), 300L, 300L)
M[sample.int(300L*300L, 300L*270L)] <- 0 # 90% of zeros
microbenchmark(
  det          = Eigen_det(M),
  absdet       = Eigen_absdet(M),
  logabsdet    = Eigen_logabsdet(M), 
  times = 200L
)



