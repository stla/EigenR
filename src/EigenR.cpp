// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

/* -------------------------------------------------------------------------- */
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic> MatrixXc;

MatrixXc matricesToMatrixXc(
  const Eigen::MatrixXd & Re, const Eigen::MatrixXd & Im 
){
  return Re.cast<std::complex<double>>() + 1i * Im.cast<std::complex<double>>();
}

Rcpp::List cplxMatrixToList(
  const MatrixXc & M
) {
  Eigen::MatrixXd realPart(M.rows(), M.cols());
  Eigen::MatrixXd imagPart(M.rows(), M.cols());
  for(auto i = 0; i < M.rows(); i++) {
    for(auto j = 0; j < M.cols(); j++) {
      const std::complex<double> z = M.coeff(i, j);
      realPart(i,j) = real(z);
      imagPart(i,j) = imag(z);
    }
  }
  return Rcpp::List::create(Rcpp::Named("real") = realPart,
                            Rcpp::Named("imag") = imagPart);
}

/* determinant -------------------------------------------------------------- */
template <typename Number>
Number determinant(
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> & M
){
  return M.determinant();
}

// [[Rcpp::export]]
double EigenR_det_real(const Eigen::MatrixXd & M) {
  return determinant<double>(M);
}

// [[Rcpp::export]]
std::complex<double> EigenR_det_cplx(
  const Eigen::MatrixXd & Re, const Eigen::MatrixXd & Im 
) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  return determinant<std::complex<double>>(M);
}
