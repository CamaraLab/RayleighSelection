// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// adjacencyCpp
List adjacencyCpp(List x, IntegerVector& feature_order);
RcppExport SEXP _RayleighSelection_adjacencyCpp(SEXP xSEXP, SEXP feature_orderSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< IntegerVector& >::type feature_order(feature_orderSEXP);
    rcpp_result_gen = Rcpp::wrap(adjacencyCpp(x, feature_order));
    return rcpp_result_gen;
END_RCPP
}
// l1down
arma::mat l1down(arma::mat one_simplices);
RcppExport SEXP _RayleighSelection_l1down(SEXP one_simplicesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type one_simplices(one_simplicesSEXP);
    rcpp_result_gen = Rcpp::wrap(l1down(one_simplices));
    return rcpp_result_gen;
END_RCPP
}
// pushCpp
List pushCpp(arma::vec v, List x, SEXP perm, arma::sp_mat adjacency);
RcppExport SEXP _RayleighSelection_pushCpp(SEXP vSEXP, SEXP xSEXP, SEXP permSEXP, SEXP adjacencySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type v(vSEXP);
    Rcpp::traits::input_parameter< List >::type x(xSEXP);
    Rcpp::traits::input_parameter< SEXP >::type perm(permSEXP);
    Rcpp::traits::input_parameter< arma::sp_mat >::type adjacency(adjacencySEXP);
    rcpp_result_gen = Rcpp::wrap(pushCpp(v, x, perm, adjacency));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_RayleighSelection_adjacencyCpp", (DL_FUNC) &_RayleighSelection_adjacencyCpp, 2},
    {"_RayleighSelection_l1down", (DL_FUNC) &_RayleighSelection_l1down, 1},
    {"_RayleighSelection_pushCpp", (DL_FUNC) &_RayleighSelection_pushCpp, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_RayleighSelection(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
