// -*- mode: C++; c-indent-level: 2; c-basic-offset: 2; indent-tabs-mode: nil;
// -*-

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

/* -------------------------------------------------------------------------- */
typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>
    MatrixXc;

typedef Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> VectorXc;

MatrixXc matricesToMatrixXc(const Eigen::MatrixXd& Re,
                            const Eigen::MatrixXd& Im) {
  return Re.cast<std::complex<double>>() + 1i * Im.cast<std::complex<double>>();
}

VectorXc vectorsToVectorXc(const Eigen::VectorXd& Re,
                           const Eigen::VectorXd& Im) {
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

Rcpp::List cplxVectorToList(const VectorXc& V) {
  Eigen::VectorXd realPart(V.size());
  Eigen::VectorXd imagPart(V.size());
  for(auto i = 0; i < V.size(); i++) {
    const std::complex<double> z = V.coeff(i);
    realPart(i) = real(z);
    imagPart(i) = imag(z);
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

/* inverse ------------------------------------------------------------------ */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> inverse(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  return M.inverse();
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_inverse_real(const Eigen::MatrixXd& M) {
  return inverse<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_inverse_cplx(const Eigen::MatrixXd& Re,
                               const Eigen::MatrixXd& Im) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  const MatrixXc Minv = inverse<std::complex<double>>(M);
  return cplxMatrixToList(Minv);
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
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> P =
      cod.colsPermutation();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> V =
      cod.matrixZ().transpose();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Kernel =
      P * V.rightCols(V.cols() - cod.rank());
  return Kernel;
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_kernel_COD_real(const Eigen::MatrixXd& M) {
  return kernel_COD<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_kernel_COD_cplx(const Eigen::MatrixXd& Re,
                                  const Eigen::MatrixXd& Im) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  const MatrixXc Kernel = kernel_COD<std::complex<double>>(M);
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
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  const MatrixXc Kernel = kernel_LU<std::complex<double>>(M);
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
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  const MatrixXc Image = image_LU<std::complex<double>>(M);
  return cplxMatrixToList(Image);
}

/* image QR ----------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> image_QR(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::ColPivHouseholderQR<
      Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      qr = M.colPivHouseholderQr();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Q =
      qr.householderQ().setLength(qr.nonzeroPivots());
  return Q.leftCols(qr.rank());
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_image_QR_real(const Eigen::MatrixXd& M) {
  return image_QR<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_image_QR_cplx(const Eigen::MatrixXd& Re,
                                const Eigen::MatrixXd& Im) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  const MatrixXc Image = image_QR<std::complex<double>>(M);
  return cplxMatrixToList(Image);
}

/* image COD ---------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> image_COD(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  Eigen::CompleteOrthogonalDecomposition<
    Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
    cod(M);
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Q =
    cod.householderQ();
  return Q.leftCols(cod.rank());
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_image_COD_real(const Eigen::MatrixXd& M) {
  return image_COD<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_image_COD_cplx(const Eigen::MatrixXd& Re,
                                  const Eigen::MatrixXd& Im) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  const MatrixXc Image = image_COD<std::complex<double>>(M);
  return cplxMatrixToList(Image);
}

/* QR ----------------------------------------------------------------------- */
template <typename Number>
std::vector<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>> QR(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::HouseholderQR<
      Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      qr = M.householderQr();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> R =
      qr.matrixQR().template triangularView<Eigen::Upper>();
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> Q =
      qr.householderQ();
  return {Q, R};
}

// [[Rcpp::export]]
Rcpp::List EigenR_QR_real(const Eigen::MatrixXd& M) {
  const std::vector<Eigen::MatrixXd> QRdecomp = QR<double>(M);
  return Rcpp::List::create(Rcpp::Named("Q") = QRdecomp[0],
                            Rcpp::Named("R") = QRdecomp[1]);
}

// [[Rcpp::export]]
Rcpp::List EigenR_QR_cplx(const Eigen::MatrixXd& Re,
                          const Eigen::MatrixXd& Im) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  const std::vector<MatrixXc> QRdecomp = QR<std::complex<double>>(M);
  return Rcpp::List::create(Rcpp::Named("Q") = cplxMatrixToList(QRdecomp[0]),
                            Rcpp::Named("R") = cplxMatrixToList(QRdecomp[1]));
}

/* Cholesky ----------------------------------------------------------------- */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> chol(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::LLT<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      lltOfM(M);
  if(lltOfM.info() != Eigen::Success) {
    throw Rcpp::exception("The matrix is not positive definite.");
  }
  return lltOfM.matrixU();
}

// [[Rcpp::export]]
Eigen::MatrixXd EigenR_chol_real(const Eigen::MatrixXd& M) {
  return chol<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_chol_cplx(const Eigen::MatrixXd& Re,
                            const Eigen::MatrixXd& Im) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  const MatrixXc U = chol<std::complex<double>>(M);
  return cplxMatrixToList(U);
}

/* UtDU --------------------------------------------------------------------- */
template <typename Number>
Rcpp::List UtDU(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& M) {
  const Eigen::LDLT<Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>>
      ldltOfM(M);
  if(ldltOfM.info() != Eigen::Success) {
    throw Rcpp::exception("Factorization failed.");
  }
  const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic> U =
      ldltOfM.matrixU();
  const Eigen::Matrix<Number, Eigen::Dynamic, 1> D = ldltOfM.vectorD();
  const Eigen::Transpositions<Eigen::Dynamic> T = ldltOfM.transpositionsP();
  Eigen::VectorXi perm(T.size());
  for(auto i = 0; i < T.size(); i++) {
    perm(i) = i;
  }
  const Rcpp::List out =
      Rcpp::List::create(Rcpp::Named("U") = U, Rcpp::Named("D") = D,
                         Rcpp::Named("perm") = T * perm);
  return out;
}

// [[Rcpp::export]]
Rcpp::List EigenR_UtDU_real(const Eigen::MatrixXd& M) {
  return UtDU<double>(M);
}

// [[Rcpp::export]]
Rcpp::List EigenR_UtDU_cplx(const Eigen::MatrixXd& Re,
                            const Eigen::MatrixXd& Im) {
  const MatrixXc M = matricesToMatrixXc(Re, Im);
  const Rcpp::List utdu = UtDU<std::complex<double>>(M);
  Rcpp::List out =
      Rcpp::List::create(Rcpp::Named("U") = cplxMatrixToList(utdu["U"]),
                         Rcpp::Named("D") = cplxVectorToList(utdu["D"]),
                         Rcpp::Named("perm") = utdu["perm"]);
  return out;
}

/* Least-squares ------------------------------------------------------------ */
template <typename Number>
Eigen::Matrix<Number, Eigen::Dynamic, 1> lsSolve(
    const Eigen::Matrix<Number, Eigen::Dynamic, Eigen::Dynamic>& A,
    const Eigen::Matrix<Number, Eigen::Dynamic, 1>& b) {
  return A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);
}

// [[Rcpp::export]]
Eigen::VectorXd EigenR_lsSolve_real(const Eigen::MatrixXd& A,
                                    const Eigen::VectorXd& b) {
  return lsSolve<double>(A, b);
}

// [[Rcpp::export]]
Rcpp::List EigenR_lsSolve_cplx(const Eigen::MatrixXd& ReA,
                               const Eigen::MatrixXd& ImA,
                               const Eigen::VectorXd& Reb,
                               const Eigen::VectorXd& Imb) {
  const MatrixXc A = matricesToMatrixXc(ReA, ImA);
  const VectorXc b = vectorsToVectorXc(Reb, Imb);
  const VectorXc v = lsSolve<std::complex<double>>(A, b);
  return cplxVectorToList(v);
}