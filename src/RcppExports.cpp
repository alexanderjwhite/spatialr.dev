// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// mat_mul_c
SEXP mat_mul_c(const Eigen::Map<Eigen::MatrixXd> A, Eigen::Map<Eigen::MatrixXd> B, int n_cores);
RcppExport SEXP _spatialr_dev_mat_mul_c(SEXP ASEXP, SEXP BSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type A(ASEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(mat_mul_c(A, B, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// exp_uv
SEXP exp_uv(const Eigen::Map<Eigen::MatrixXd> u, Eigen::Map<Eigen::MatrixXd> v, int n_cores);
RcppExport SEXP _spatialr_dev_exp_uv(SEXP uSEXP, SEXP vSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type u(uSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type v(vSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(exp_uv(u, v, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// grad_v
SEXP grad_v(const Eigen::Map<Eigen::MatrixXd> x, Eigen::Map<Eigen::MatrixXd> u, Eigen::Map<Eigen::MatrixXd> v, Eigen::Map<Eigen::MatrixXd> uv_exp, int n_cores);
RcppExport SEXP _spatialr_dev_grad_v(SEXP xSEXP, SEXP uSEXP, SEXP vSEXP, SEXP uv_expSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type u(uSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type v(vSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type uv_exp(uv_expSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_v(x, u, v, uv_exp, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// grad_u
SEXP grad_u(const Eigen::Map<Eigen::MatrixXd> x, Eigen::Map<Eigen::MatrixXd> u, Eigen::Map<Eigen::MatrixXd> v, Eigen::Map<Eigen::MatrixXd> uv_exp, Eigen::Map<Eigen::MatrixXd> w, Eigen::Map<Eigen::MatrixXd> j, Eigen::Map<Eigen::MatrixXd> one, double lambda, int n_cores);
RcppExport SEXP _spatialr_dev_grad_u(SEXP xSEXP, SEXP uSEXP, SEXP vSEXP, SEXP uv_expSEXP, SEXP wSEXP, SEXP jSEXP, SEXP oneSEXP, SEXP lambdaSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type u(uSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type v(vSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type uv_exp(uv_expSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type w(wSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type j(jSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type one(oneSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(grad_u(x, u, v, uv_exp, w, j, one, lambda, n_cores));
    return rcpp_result_gen;
END_RCPP
}
// objective_c
SEXP objective_c(const Eigen::Map<Eigen::MatrixXd> x, Eigen::Map<Eigen::MatrixXd> u, Eigen::Map<Eigen::MatrixXd> v, Eigen::Map<Eigen::MatrixXd> w, Eigen::Map<Eigen::MatrixXd> j, double lambda, int n_cores);
RcppExport SEXP _spatialr_dev_objective_c(SEXP xSEXP, SEXP uSEXP, SEXP vSEXP, SEXP wSEXP, SEXP jSEXP, SEXP lambdaSEXP, SEXP n_coresSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::Map<Eigen::MatrixXd> >::type x(xSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type u(uSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type v(vSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type w(wSEXP);
    Rcpp::traits::input_parameter< Eigen::Map<Eigen::MatrixXd> >::type j(jSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< int >::type n_cores(n_coresSEXP);
    rcpp_result_gen = Rcpp::wrap(objective_c(x, u, v, w, j, lambda, n_cores));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spatialr_dev_mat_mul_c", (DL_FUNC) &_spatialr_dev_mat_mul_c, 3},
    {"_spatialr_dev_exp_uv", (DL_FUNC) &_spatialr_dev_exp_uv, 3},
    {"_spatialr_dev_grad_v", (DL_FUNC) &_spatialr_dev_grad_v, 5},
    {"_spatialr_dev_grad_u", (DL_FUNC) &_spatialr_dev_grad_u, 9},
    {"_spatialr_dev_objective_c", (DL_FUNC) &_spatialr_dev_objective_c, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_spatialr_dev(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
