#include "EigenR.h"

// [[Rcpp::export]]
Rcpp::List EigenR_realQZ(const Eigen::MatrixXd& A, const Eigen::MatrixXd& B) {
  Eigen::RealQZ<Eigen::MatrixXd> qz(A.rows()); // preallocate space 
  qz.compute(A, B);  // A = Q S Z,  B = Q T Z
  const Eigen::MatrixXd Q = qz.matrixQ();
  const Eigen::MatrixXd Z = qz.matrixZ();
  const Eigen::MatrixXd S = qz.matrixS();
  const Eigen::MatrixXd T = qz.matrixT();
  return Rcpp::List::create(Rcpp::Named("Q") = Q,
                            Rcpp::Named("Z") = Z,
                            Rcpp::Named("S") = S,
                            Rcpp::Named("T") = T);
}
