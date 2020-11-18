library(EigenR)
library(XRJulia)
library(microbenchmark)


set.seed(666)
M <- round(matrix(rgamma(12, 10, 1), nrow = 4, ncol = 3), digits = 1)
M <- cbind(M, M[,1]+M[,2], M[,2]+2*M[,3])
Mjulia <- sprintf("[%s]", paste0(
  apply(M, 1, function(row) paste0(format(row, nsmall = 1), collapse = " ")),
  collapse = "; "
))

code1 <- paste0(c(
  'function Kernel(M)',
  '  @rput M',
  '  R"K = EigenR::Eigen_kernel(M)"',
  '  return @rget K',
  'end'
), collapse = "\n")
code2 <- sprintf("Kernel(%s)", Mjulia)

juliaEval("using RCall")
juliaEval(code1)
x <- juliaEval(code2)
juliaGet(x) 

microbenchmark(
  Eigen = Eigen_kernel(M),
  Julia = juliaEval(code2),
  times = 2
)
