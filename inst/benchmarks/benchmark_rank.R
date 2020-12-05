library(EigenR)
library(microbenchmark)

set.seed(666)
nrows <- 10
ncols <- 8
M <- round(matrix(rgamma(nrows*ncols, 10, 1), nrow = nrows, ncol = ncols), digits = 1)
M <- cbind(M, M[,1]+M[,2], M[,2]+2*M[,3])

microbenchmark(
  EigenR = Eigen_rank(M),
  base = qr(M)$rank,
  times = 50
)

