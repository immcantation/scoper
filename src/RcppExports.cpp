// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// pairwiseMutMatrixRcpp
List pairwiseMutMatrixRcpp(NumericVector informative_pos, StringMatrix mutMtx, NumericMatrix motifMtx);
RcppExport SEXP _scoper_pairwiseMutMatrixRcpp(SEXP informative_posSEXP, SEXP mutMtxSEXP, SEXP motifMtxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type informative_pos(informative_posSEXP);
    Rcpp::traits::input_parameter< StringMatrix >::type mutMtx(mutMtxSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type motifMtx(motifMtxSEXP);
    rcpp_result_gen = Rcpp::wrap(pairwiseMutMatrixRcpp(informative_pos, mutMtx, motifMtx));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_scoper_pairwiseMutMatrixRcpp", (DL_FUNC) &_scoper_pairwiseMutMatrixRcpp, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_scoper(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
