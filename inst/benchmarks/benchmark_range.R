library(EigenR)
library(microbenchmark)

set.seed(666)
nrows <- 8
ncols <- 4
M <- round(matrix(rgamma(nrows*ncols, 10, 1), nrow = nrows, ncol = ncols), digits = 1)
M <- cbind(M, M[,1]+M[,2], M[,2]+2*M[,3])

microbenchmark(
  QR = Eigen_range(M, "QR"),
  LU = Eigen_range(M, "LU"),
  LU = Eigen_range(M, "COD"),
  times = 20
)

