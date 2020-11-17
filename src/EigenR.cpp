// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil;
// -*-

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

/* -------------------------------------------------------------------------- */
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
    MatrixXc;

MatrixXc matricesToMatrixXc(const Eigen::MatrixXd& Re,
                            const Eigen::MatrixXd& Im) {
  return Re.cast<std::complex<double>>() + 1i * Im.cast<std::complex<double>>();
}

Rcpp::List cplxMatrixToList(const MatrixXc& M) {
  Eigen::MatrixXd realPart(M.rows(), M.cols());
  Eigen::MatrixXd imagPart(M.rows(), M.cols());
  for(auto i = 0; i < M.rows(); i++) {
    for(auto j = 0; j < M.cols(); j++) {
      const std::complex<double> z = M.coeff(i, j);
      realPart(i, j) = real(z);
      imagPart(i, j) = imag(z);
    }
  }
  return Rcpp::List::create(Rcpp::Named("real") = realPart,
                            Rcpp::Named("imag") = imagPart);
}

/* determinant -------------------------------------------------------------- */
template <typename Number>
Number determinant(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  return M.determinant();
}

// [[Rcpp::export]]
double EigenR_det_real(const Eigen::MatrixXd& M) {
  return determinant<double>(M);
}

// [[Rcpp::export]]
std::complex<double> EigenR_det_cplx(const Eigen::MatrixXd& Re,
                                     const Eigen::MatrixXd& Im) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  return determinant<std::complex<double>>(M);
}

/* rank --------------------------------------------------------------------- */
template <typename Number>
unsigned rank(const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  return M.colPivHouseholderQr().rank();
}

// [[Rcpp::export]]
unsigned EigenR_rank_real(const Eigen::MatrixXd& M) {
  return rank<double>(M);
}

// [[Rcpp::export]]
unsigned EigenR_rank_cplx(const Eigen::MatrixXd& Re,
                          const Eigen::MatrixXd& Im) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  return rank<std::complex<double>>(M);
}

/* kernel COD --------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> kernel_COD(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  // https://stackoverflow.com/a/53598471/1100107
  Eigen::CompleteOrthogonalDecomposition<
      Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      cod;
  cod.compute(M);
  // Find URV^T
  unsigned rk = cod.rank();
  Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> P =
      cod.colsPermutation();
  Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> V =
      cod.matrixZ().transpose();
  Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Kernel =
      P * V.block(0, rk, V.rows(), V.cols() - rk);
  return Kernel;
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_kernel_COD_real(const Eigen::MatrixXd& M) {
  return kernel_COD<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_kernel_COD_cplx(const Eigen::MatrixXd& Re,
                                  const Eigen::MatrixXd& Im) {
  MatrixXc M = matricesToMatrixXc(Re, Im);
  MatrixXc Kernel = kernel_COD<std::complex<double>>(M);
  return cplxMatrixToList(Kernel);
}

/* kernel LU ---------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> kernel_LU(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::FullPivLU<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      lu(M);
  return lu.kernel();
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_kernel_LU_real(const Eigen::MatrixXd& M) {
  return kernel_LU<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_kernel_LU_cplx(const Eigen::MatrixXd& Re,
                                 const Eigen::MatrixXd& Im) {
  MatrixXc M = matricesToMatrixXc(Re, Im);
  MatrixXc Kernel = kernel_LU<std::complex<double>>(M);
  return cplxMatrixToList(Kernel);
}

/* image LU ----------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> image_LU(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::FullPivLU<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
    lu(M);
  return lu.image(M);
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_image_LU_real(const Eigen::MatrixXd& M) {
  return image_LU<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_image_LU_cplx(const Eigen::MatrixXd& Re,
                                 const Eigen::MatrixXd& Im) {
  MatrixXc M = matricesToMatrixXc(Re, Im);
  MatrixXc Image = image_LU<std::complex<double>>(M);
  return cplxMatrixToList(Image);
}
