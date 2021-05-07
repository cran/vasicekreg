// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// cpp_dvasicekmean
NumericVector cpp_dvasicekmean(const NumericVector x, const NumericVector alpha, const NumericVector theta, const bool logprob);
RcppExport SEXP _vasicekreg_cpp_dvasicekmean(SEXP xSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP logprobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const bool >::type logprob(logprobSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dvasicekmean(x, alpha, theta, logprob));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pvasicekmean
NumericVector cpp_pvasicekmean(const NumericVector x, const NumericVector alpha, const NumericVector theta, const bool lowertail, const bool logprob);
RcppExport SEXP _vasicekreg_cpp_pvasicekmean(SEXP xSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP lowertailSEXP, SEXP logprobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const bool >::type lowertail(lowertailSEXP);
    Rcpp::traits::input_parameter< const bool >::type logprob(logprobSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pvasicekmean(x, alpha, theta, lowertail, logprob));
    return rcpp_result_gen;
END_RCPP
}
// cpp_qvasicekmean
NumericVector cpp_qvasicekmean(const NumericVector x, const NumericVector alpha, const NumericVector theta, const bool lowertail, const bool logprob);
RcppExport SEXP _vasicekreg_cpp_qvasicekmean(SEXP xSEXP, SEXP alphaSEXP, SEXP thetaSEXP, SEXP lowertailSEXP, SEXP logprobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const bool >::type lowertail(lowertailSEXP);
    Rcpp::traits::input_parameter< const bool >::type logprob(logprobSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_qvasicekmean(x, alpha, theta, lowertail, logprob));
    return rcpp_result_gen;
END_RCPP
}
// cpp_dvasicekquant
NumericVector cpp_dvasicekquant(const NumericVector x, const NumericVector mu, const NumericVector theta, const NumericVector tau, const bool logprob);
RcppExport SEXP _vasicekreg_cpp_dvasicekquant(SEXP xSEXP, SEXP muSEXP, SEXP thetaSEXP, SEXP tauSEXP, SEXP logprobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool >::type logprob(logprobSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_dvasicekquant(x, mu, theta, tau, logprob));
    return rcpp_result_gen;
END_RCPP
}
// cpp_pvasicekquant
NumericVector cpp_pvasicekquant(const NumericVector x, const NumericVector mu, const NumericVector theta, const NumericVector tau, const bool lowertail, const bool logprob);
RcppExport SEXP _vasicekreg_cpp_pvasicekquant(SEXP xSEXP, SEXP muSEXP, SEXP thetaSEXP, SEXP tauSEXP, SEXP lowertailSEXP, SEXP logprobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool >::type lowertail(lowertailSEXP);
    Rcpp::traits::input_parameter< const bool >::type logprob(logprobSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_pvasicekquant(x, mu, theta, tau, lowertail, logprob));
    return rcpp_result_gen;
END_RCPP
}
// cpp_qvasicekquant
NumericVector cpp_qvasicekquant(const NumericVector x, const NumericVector mu, const NumericVector theta, const NumericVector tau, const bool lowertail, const bool logprob);
RcppExport SEXP _vasicekreg_cpp_qvasicekquant(SEXP xSEXP, SEXP muSEXP, SEXP thetaSEXP, SEXP tauSEXP, SEXP lowertailSEXP, SEXP logprobSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< const NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const bool >::type lowertail(lowertailSEXP);
    Rcpp::traits::input_parameter< const bool >::type logprob(logprobSEXP);
    rcpp_result_gen = Rcpp::wrap(cpp_qvasicekquant(x, mu, theta, tau, lowertail, logprob));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_vasicekreg_cpp_dvasicekmean", (DL_FUNC) &_vasicekreg_cpp_dvasicekmean, 4},
    {"_vasicekreg_cpp_pvasicekmean", (DL_FUNC) &_vasicekreg_cpp_pvasicekmean, 5},
    {"_vasicekreg_cpp_qvasicekmean", (DL_FUNC) &_vasicekreg_cpp_qvasicekmean, 5},
    {"_vasicekreg_cpp_dvasicekquant", (DL_FUNC) &_vasicekreg_cpp_dvasicekquant, 5},
    {"_vasicekreg_cpp_pvasicekquant", (DL_FUNC) &_vasicekreg_cpp_pvasicekquant, 6},
    {"_vasicekreg_cpp_qvasicekquant", (DL_FUNC) &_vasicekreg_cpp_qvasicekquant, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_vasicekreg(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
