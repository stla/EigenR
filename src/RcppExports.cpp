// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// EigenR_det_real
double EigenR_det_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_det_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_det_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_det_cplx
std::complex<double> EigenR_det_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_det_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_det_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_det_sparse_real
double EigenR_det_sparse_real(const std::vector<size_t>& i, const std::vector<size_t>& j, const std::vector<double>& Mij, const size_t nrows, const size_t ncols);
RcppExport SEXP _EigenR_EigenR_det_sparse_real(SEXP iSEXP, SEXP jSEXP, SEXP MijSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type Mij(MijSEXP);
    Rcpp::traits::input_parameter< const size_t >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_det_sparse_real(i, j, Mij, nrows, ncols));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_det_sparse_cplx
std::complex<double> EigenR_det_sparse_cplx(const std::vector<size_t>& i, const std::vector<size_t>& j, const std::vector<std::complex<double>>& Mij, const size_t nrows, const size_t ncols);
RcppExport SEXP _EigenR_EigenR_det_sparse_cplx(SEXP iSEXP, SEXP jSEXP, SEXP MijSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::complex<double>>& >::type Mij(MijSEXP);
    Rcpp::traits::input_parameter< const size_t >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_det_sparse_cplx(i, j, Mij, nrows, ncols));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_rank_real
unsigned EigenR_rank_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_rank_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_rank_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_rank_cplx
unsigned EigenR_rank_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_rank_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_rank_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_inverse_real
Eigen::MatrixXd EigenR_inverse_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_inverse_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_inverse_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_inverse_cplx
Rcpp::List EigenR_inverse_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_inverse_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_inverse_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_kernel_COD_real
Eigen::MatrixXd EigenR_kernel_COD_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_kernel_COD_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_kernel_COD_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_kernel_COD_cplx
Rcpp::List EigenR_kernel_COD_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_kernel_COD_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_kernel_COD_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_kernel_LU_real
Eigen::MatrixXd EigenR_kernel_LU_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_kernel_LU_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_kernel_LU_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_kernel_LU_cplx
Rcpp::List EigenR_kernel_LU_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_kernel_LU_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_kernel_LU_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_image_LU_real
Eigen::MatrixXd EigenR_image_LU_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_image_LU_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_image_LU_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_image_LU_cplx
Rcpp::List EigenR_image_LU_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_image_LU_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_image_LU_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_image_QR_real
Eigen::MatrixXd EigenR_image_QR_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_image_QR_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_image_QR_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_image_QR_cplx
Rcpp::List EigenR_image_QR_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_image_QR_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_image_QR_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_image_COD_real
Eigen::MatrixXd EigenR_image_COD_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_image_COD_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_image_COD_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_image_COD_cplx
Rcpp::List EigenR_image_COD_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_image_COD_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_image_COD_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_QR_real
Rcpp::List EigenR_QR_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_QR_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_QR_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_QR_cplx
Rcpp::List EigenR_QR_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_QR_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_QR_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_chol_real
Eigen::MatrixXd EigenR_chol_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_chol_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_chol_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_chol_cplx
Rcpp::List EigenR_chol_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_chol_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_chol_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_chol_sparse_real
Eigen::MatrixXd EigenR_chol_sparse_real(const std::vector<size_t>& i, const std::vector<size_t>& j, const std::vector<double>& Mij, const size_t nrows, const size_t ncols);
RcppExport SEXP _EigenR_EigenR_chol_sparse_real(SEXP iSEXP, SEXP jSEXP, SEXP MijSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const std::vector<double>& >::type Mij(MijSEXP);
    Rcpp::traits::input_parameter< const size_t >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_chol_sparse_real(i, j, Mij, nrows, ncols));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_chol_sparse_cplx
Rcpp::List EigenR_chol_sparse_cplx(const std::vector<size_t>& i, const std::vector<size_t>& j, const std::vector<std::complex<double>>& Mij, const size_t nrows, const size_t ncols);
RcppExport SEXP _EigenR_EigenR_chol_sparse_cplx(SEXP iSEXP, SEXP jSEXP, SEXP MijSEXP, SEXP nrowsSEXP, SEXP ncolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type i(iSEXP);
    Rcpp::traits::input_parameter< const std::vector<size_t>& >::type j(jSEXP);
    Rcpp::traits::input_parameter< const std::vector<std::complex<double>>& >::type Mij(MijSEXP);
    Rcpp::traits::input_parameter< const size_t >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< const size_t >::type ncols(ncolsSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_chol_sparse_cplx(i, j, Mij, nrows, ncols));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_UtDU_real
Rcpp::List EigenR_UtDU_real(const Eigen::MatrixXd& M);
RcppExport SEXP _EigenR_EigenR_UtDU_real(SEXP MSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type M(MSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_UtDU_real(M));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_UtDU_cplx
Rcpp::List EigenR_UtDU_cplx(const Eigen::MatrixXd& Re, const Eigen::MatrixXd& Im);
RcppExport SEXP _EigenR_EigenR_UtDU_cplx(SEXP ReSEXP, SEXP ImSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Re(ReSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type Im(ImSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_UtDU_cplx(Re, Im));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_lsSolve_real
Eigen::VectorXd EigenR_lsSolve_real(const Eigen::MatrixXd& A, const Eigen::VectorXd& b);
RcppExport SEXP _EigenR_EigenR_lsSolve_real(SEXP ASEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_lsSolve_real(A, b));
    return rcpp_result_gen;
END_RCPP
}
// EigenR_lsSolve_cplx
Rcpp::List EigenR_lsSolve_cplx(const Eigen::MatrixXd& ReA, const Eigen::MatrixXd& ImA, const Eigen::VectorXd& Reb, const Eigen::VectorXd& Imb);
RcppExport SEXP _EigenR_EigenR_lsSolve_cplx(SEXP ReASEXP, SEXP ImASEXP, SEXP RebSEXP, SEXP ImbSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type ReA(ReASEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXd& >::type ImA(ImASEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Reb(RebSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Imb(ImbSEXP);
    rcpp_result_gen = Rcpp::wrap(EigenR_lsSolve_cplx(ReA, ImA, Reb, Imb));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_EigenR_EigenR_det_real", (DL_FUNC) &_EigenR_EigenR_det_real, 1},
    {"_EigenR_EigenR_det_cplx", (DL_FUNC) &_EigenR_EigenR_det_cplx, 2},
    {"_EigenR_EigenR_det_sparse_real", (DL_FUNC) &_EigenR_EigenR_det_sparse_real, 5},
    {"_EigenR_EigenR_det_sparse_cplx", (DL_FUNC) &_EigenR_EigenR_det_sparse_cplx, 5},
    {"_EigenR_EigenR_rank_real", (DL_FUNC) &_EigenR_EigenR_rank_real, 1},
    {"_EigenR_EigenR_rank_cplx", (DL_FUNC) &_EigenR_EigenR_rank_cplx, 2},
    {"_EigenR_EigenR_inverse_real", (DL_FUNC) &_EigenR_EigenR_inverse_real, 1},
    {"_EigenR_EigenR_inverse_cplx", (DL_FUNC) &_EigenR_EigenR_inverse_cplx, 2},
    {"_EigenR_EigenR_kernel_COD_real", (DL_FUNC) &_EigenR_EigenR_kernel_COD_real, 1},
    {"_EigenR_EigenR_kernel_COD_cplx", (DL_FUNC) &_EigenR_EigenR_kernel_COD_cplx, 2},
    {"_EigenR_EigenR_kernel_LU_real", (DL_FUNC) &_EigenR_EigenR_kernel_LU_real, 1},
    {"_EigenR_EigenR_kernel_LU_cplx", (DL_FUNC) &_EigenR_EigenR_kernel_LU_cplx, 2},
    {"_EigenR_EigenR_image_LU_real", (DL_FUNC) &_EigenR_EigenR_image_LU_real, 1},
    {"_EigenR_EigenR_image_LU_cplx", (DL_FUNC) &_EigenR_EigenR_image_LU_cplx, 2},
    {"_EigenR_EigenR_image_QR_real", (DL_FUNC) &_EigenR_EigenR_image_QR_real, 1},
    {"_EigenR_EigenR_image_QR_cplx", (DL_FUNC) &_EigenR_EigenR_image_QR_cplx, 2},
    {"_EigenR_EigenR_image_COD_real", (DL_FUNC) &_EigenR_EigenR_image_COD_real, 1},
    {"_EigenR_EigenR_image_COD_cplx", (DL_FUNC) &_EigenR_EigenR_image_COD_cplx, 2},
    {"_EigenR_EigenR_QR_real", (DL_FUNC) &_EigenR_EigenR_QR_real, 1},
    {"_EigenR_EigenR_QR_cplx", (DL_FUNC) &_EigenR_EigenR_QR_cplx, 2},
    {"_EigenR_EigenR_chol_real", (DL_FUNC) &_EigenR_EigenR_chol_real, 1},
    {"_EigenR_EigenR_chol_cplx", (DL_FUNC) &_EigenR_EigenR_chol_cplx, 2},
    {"_EigenR_EigenR_chol_sparse_real", (DL_FUNC) &_EigenR_EigenR_chol_sparse_real, 5},
    {"_EigenR_EigenR_chol_sparse_cplx", (DL_FUNC) &_EigenR_EigenR_chol_sparse_cplx, 5},
    {"_EigenR_EigenR_UtDU_real", (DL_FUNC) &_EigenR_EigenR_UtDU_real, 1},
    {"_EigenR_EigenR_UtDU_cplx", (DL_FUNC) &_EigenR_EigenR_UtDU_cplx, 2},
    {"_EigenR_EigenR_lsSolve_real", (DL_FUNC) &_EigenR_EigenR_lsSolve_real, 2},
    {"_EigenR_EigenR_lsSolve_cplx", (DL_FUNC) &_EigenR_EigenR_lsSolve_cplx, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_EigenR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
