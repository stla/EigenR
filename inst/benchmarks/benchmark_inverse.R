library(EigenR)
library(microbenchmark)

set.seed(666)
size <- 100
M <- matrix(rgamma(size*size, 10, 1), nrow = size, ncol = size)

microbenchmark(
  EigenR = Eigen_inverse(M+1i*M),
  base = solve(M+1i*M),
  times = 50
)

