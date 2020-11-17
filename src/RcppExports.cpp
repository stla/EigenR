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

static const R_CallMethodDef CallEntries[] = {
    {"_EigenR_EigenR_det_real", (DL_FUNC) &_EigenR_EigenR_det_real, 1},
    {"_EigenR_EigenR_det_cplx", (DL_FUNC) &_EigenR_EigenR_det_cplx, 2},
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
    {NULL, NULL, 0}
};

RcppExport void R_init_EigenR(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
