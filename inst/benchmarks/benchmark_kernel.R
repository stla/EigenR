library(EigenR)
library(microbenchmark)

set.seed(666)
nrows <- 1000
ncols <- 900
M <- round(matrix(rgamma(nrows*ncols, 10, 1), nrow = nrows, ncol = ncols), digits = 1)
M <- cbind(M, M[,1]+M[,2], M[,2]+2*M[,3])

microbenchmark(
  COD = Eigen_kernel(M, "COD"),
  LU = Eigen_kernel(M, "LU"),
  times = 50
)

