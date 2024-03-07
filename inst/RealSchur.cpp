#include "EigenR.h"

// [[Rcpp::export]]
Rcpp::List EigenR_realSchur(const Eigen::MatrixXd& M) {
  Eigen::RealSchur<Eigen::MatrixXd> schur(M); 
  const Eigen::MatrixXd U = schur.matrixU();
  const Eigen::MatrixXd T = schur.matrixT();
  return Rcpp::List::create(Rcpp::Named("U") = U,
                            Rcpp::Named("T") = T);
}


// [[Rcpp::export]]
Rcpp::List EigenR_complexSchur(const Eigen::MatrixXd& Re, 
	                           const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  Eigen::ComplexSchur<MatrixXcd> schur(M.rows());
  schur.compute(M);
  const Eigen::MatrixXcd U = schur.matrixU();
  const Eigen::MatrixXcd T = schur.matrixT();
  return Rcpp::List::create(Rcpp::Named("U") = cplxMatrixToList(U),
                            Rcpp::Named("T") = cplxMatrixToList(T));
}


std::complex<double> EigenR_det_cplx(const Eigen::MatrixXd& Re,
                                     const Eigen::MatrixXd& Im) {
  const Eigen::MatrixXcd M = matricesToMatrixXcd(Re, Im);
  return determinant<std::complex<double>>(M);
}


MatrixXd A = MatrixXd::Random(6,6);
cout << "Here is a random 6x6 matrix, A:" << endl << A << endl << endl;
 
RealSchur<MatrixXd> schur(A);
cout << "The orthogonal matrix U is:" << endl << schur.matrixU() << endl;
cout << "The quasi-triangular matrix T is:" << endl << schur.matrixT() << endl << endl;
 
MatrixXd U = schur.matrixU();
MatrixXd T = schur.matrixT();