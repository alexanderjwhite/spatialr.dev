// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// vecnorm_row
arma::rowvec vecnorm_row(arma::mat x);
RcppExport SEXP _spatialr_dev_vecnorm_row(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vecnorm_row(x));
    return rcpp_result_gen;
END_RCPP
}
// vecnorm_diag
arma::mat vecnorm_diag(arma::mat x);
RcppExport SEXP _spatialr_dev_vecnorm_diag(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(vecnorm_diag(x));
    return rcpp_result_gen;
END_RCPP
}
// uv_norm
List uv_norm(arma::mat u, arma::mat v);
RcppExport SEXP _spatialr_dev_uv_norm(SEXP uSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(uv_norm(u, v));
    return rcpp_result_gen;
END_RCPP
}
// fct_c_opt_adam
List fct_c_opt_adam(arma::mat gradients, List state);
RcppExport SEXP _spatialr_dev_fct_c_opt_adam(SEXP gradientsSEXP, SEXP stateSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type gradients(gradientsSEXP);
    Rcpp::traits::input_parameter< List >::type state(stateSEXP);
    rcpp_result_gen = Rcpp::wrap(fct_c_opt_adam(gradients, state));
    return rcpp_result_gen;
END_RCPP
}
// pdist
double pdist(arma::mat x, arma::sp_mat w, arma::mat index);
RcppExport SEXP _spatialr_dev_pdist(SEXP xSEXP, SEXP wSEXP, SEXP indexSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type index(indexSEXP);
    rcpp_result_gen = Rcpp::wrap(pdist(x, w, index));
    return rcpp_result_gen;
END_RCPP
}
// fct_c_optimize
List fct_c_optimize(arma::sp_mat x, arma::mat u, arma::mat v, arma::sp_mat w, arma::mat index, double lambda, Nullable<double> cnorm, double epsilon, int maxiter);
RcppExport SEXP _spatialr_dev_fct_c_optimize(SEXP xSEXP, SEXP uSEXP, SEXP vSEXP, SEXP wSEXP, SEXP indexSEXP, SEXP lambdaSEXP, SEXP cnormSEXP, SEXP epsilonSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type w(wSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type index(indexSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< Nullable<double> >::type cnorm(cnormSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< int >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(fct_c_optimize(x, u, v, w, index, lambda, cnorm, epsilon, maxiter));
    return rcpp_result_gen;
END_RCPP
}
// obs_log_like
arma::rowvec obs_log_like(arma::sp_mat test_x, arma::mat u, arma::mat v, arma::mat test_nn, arma::mat index_map);
RcppExport SEXP _spatialr_dev_obs_log_like(SEXP test_xSEXP, SEXP uSEXP, SEXP vSEXP, SEXP test_nnSEXP, SEXP index_mapSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::sp_mat >::type test_x(test_xSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type u(uSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type v(vSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type test_nn(test_nnSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type index_map(index_mapSEXP);
    rcpp_result_gen = Rcpp::wrap(obs_log_like(test_x, u, v, test_nn, index_map));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_spatialr_dev_vecnorm_row", (DL_FUNC) &_spatialr_dev_vecnorm_row, 1},
    {"_spatialr_dev_vecnorm_diag", (DL_FUNC) &_spatialr_dev_vecnorm_diag, 1},
    {"_spatialr_dev_uv_norm", (DL_FUNC) &_spatialr_dev_uv_norm, 2},
    {"_spatialr_dev_fct_c_opt_adam", (DL_FUNC) &_spatialr_dev_fct_c_opt_adam, 2},
    {"_spatialr_dev_pdist", (DL_FUNC) &_spatialr_dev_pdist, 3},
    {"_spatialr_dev_fct_c_optimize", (DL_FUNC) &_spatialr_dev_fct_c_optimize, 9},
    {"_spatialr_dev_obs_log_like", (DL_FUNC) &_spatialr_dev_obs_log_like, 5},
    {NULL, NULL, 0}
};

RcppExport void R_init_spatialr_dev(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
