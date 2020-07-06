// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// func_dc
List func_dc(NumericMatrix z);
RcppExport SEXP _BSCOV_func_dc(SEXP zSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    rcpp_result_gen = Rcpp::wrap(func_dc(z));
    return rcpp_result_gen;
END_RCPP
}
// func_dc_by
List func_dc_by(NumericMatrix z, double dmby, double dtby);
RcppExport SEXP _BSCOV_func_dc_by(SEXP zSEXP, SEXP dmbySEXP, SEXP dtbySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type dmby(dmbySEXP);
    Rcpp::traits::input_parameter< double >::type dtby(dtbySEXP);
    rcpp_result_gen = Rcpp::wrap(func_dc_by(z, dmby, dtby));
    return rcpp_result_gen;
END_RCPP
}
// func_coef
NumericMatrix func_coef(NumericMatrix z, int scale);
RcppExport SEXP _BSCOV_func_coef(SEXP zSEXP, SEXP scaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type z(zSEXP);
    Rcpp::traits::input_parameter< int >::type scale(scaleSEXP);
    rcpp_result_gen = Rcpp::wrap(func_coef(z, scale));
    return rcpp_result_gen;
END_RCPP
}
// func_input
NumericMatrix func_input(NumericMatrix coef, NumericMatrix sgn);
RcppExport SEXP _BSCOV_func_input(SEXP coefSEXP, SEXP sgnSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coef(coefSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type sgn(sgnSEXP);
    rcpp_result_gen = Rcpp::wrap(func_input(coef, sgn));
    return rcpp_result_gen;
END_RCPP
}
// func_input_on
NumericMatrix func_input_on(NumericMatrix coef);
RcppExport SEXP _BSCOV_func_input_on(SEXP coefSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type coef(coefSEXP);
    rcpp_result_gen = Rcpp::wrap(func_input_on(coef));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BSCOV_func_dc", (DL_FUNC) &_BSCOV_func_dc, 1},
    {"_BSCOV_func_dc_by", (DL_FUNC) &_BSCOV_func_dc_by, 3},
    {"_BSCOV_func_coef", (DL_FUNC) &_BSCOV_func_coef, 2},
    {"_BSCOV_func_input", (DL_FUNC) &_BSCOV_func_input, 2},
    {"_BSCOV_func_input_on", (DL_FUNC) &_BSCOV_func_input_on, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BSCOV(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
