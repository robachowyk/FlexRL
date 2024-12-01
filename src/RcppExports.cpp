// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// callFunction
NumericVector callFunction(NumericVector x, Function f);
RcppExport SEXP _FlexRL_callFunction(SEXP xSEXP, SEXP fSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Function >::type f(fSEXP);
    rcpp_result_gen = Rcpp::wrap(callFunction(x, f));
    return rcpp_result_gen;
END_RCPP
}
// F2
List F2(IntegerVector U, int nvals);
RcppExport SEXP _FlexRL_F2(SEXP USEXP, SEXP nvalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type U(USEXP);
    Rcpp::traits::input_parameter< int >::type nvals(nvalsSEXP);
    rcpp_result_gen = Rcpp::wrap(F2(U, nvals));
    return rcpp_result_gen;
END_RCPP
}
// F33
IntegerMatrix F33(List A, List B, int nvals);
RcppExport SEXP _FlexRL_F33(SEXP ASEXP, SEXP BSEXP, SEXP nvalsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type A(ASEXP);
    Rcpp::traits::input_parameter< List >::type B(BSEXP);
    Rcpp::traits::input_parameter< int >::type nvals(nvalsSEXP);
    rcpp_result_gen = Rcpp::wrap(F33(A, B, nvals));
    return rcpp_result_gen;
END_RCPP
}
// sspaste2
CharacterVector sspaste2(IntegerMatrix A);
RcppExport SEXP _FlexRL_sspaste2(SEXP ASEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerMatrix >::type A(ASEXP);
    rcpp_result_gen = Rcpp::wrap(sspaste2(A));
    return rcpp_result_gen;
END_RCPP
}
// initΔMap
void initΔMap();
RcppExport SEXP _FlexRL_initΔMap() {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    initΔMap();
    return R_NilValue;
END_RCPP
}
// Δfind
IntegerMatrix Δfind();
RcppExport SEXP _FlexRL_Δfind() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(Δfind());
    return rcpp_result_gen;
END_RCPP
}
// sampleD
List sampleD(const Rcpp::IntegerMatrix& S, const Rcpp::NumericVector& LLA, const Rcpp::NumericVector& LLB, const Rcpp::NumericVector& LLL, const Rcpp::NumericVector& gamma, double loglik, int nlinkrec, LogicalVector& sumRowD, LogicalVector& sumColD);
RcppExport SEXP _FlexRL_sampleD(SEXP SSEXP, SEXP LLASEXP, SEXP LLBSEXP, SEXP LLLSEXP, SEXP gammaSEXP, SEXP loglikSEXP, SEXP nlinkrecSEXP, SEXP sumRowDSEXP, SEXP sumColDSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::IntegerMatrix& >::type S(SSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type LLA(LLASEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type LLB(LLBSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type LLL(LLLSEXP);
    Rcpp::traits::input_parameter< const Rcpp::NumericVector& >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< double >::type loglik(loglikSEXP);
    Rcpp::traits::input_parameter< int >::type nlinkrec(nlinkrecSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type sumRowD(sumRowDSEXP);
    Rcpp::traits::input_parameter< LogicalVector& >::type sumColD(sumColDSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleD(S, LLA, LLB, LLL, gamma, loglik, nlinkrec, sumRowD, sumColD));
    return rcpp_result_gen;
END_RCPP
}
// sampleNL
IntegerVector sampleNL(IntegerVector G, NumericVector eta, NumericVector phi);
RcppExport SEXP _FlexRL_sampleNL(SEXP GSEXP, SEXP etaSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type G(GSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleNL(G, eta, phi));
    return rcpp_result_gen;
END_RCPP
}
// sampleL
IntegerVector sampleL(IntegerVector GA, IntegerVector GB, NumericVector survivalpSameH, IntegerMatrix choice_set, IntegerVector choice_equal, int nval, NumericVector phikA, NumericVector phikB, NumericVector eta);
RcppExport SEXP _FlexRL_sampleL(SEXP GASEXP, SEXP GBSEXP, SEXP survivalpSameHSEXP, SEXP choice_setSEXP, SEXP choice_equalSEXP, SEXP nvalSEXP, SEXP phikASEXP, SEXP phikBSEXP, SEXP etaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type GA(GASEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type GB(GBSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type survivalpSameH(survivalpSameHSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type choice_set(choice_setSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type choice_equal(choice_equalSEXP);
    Rcpp::traits::input_parameter< int >::type nval(nvalSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phikA(phikASEXP);
    Rcpp::traits::input_parameter< NumericVector >::type phikB(phikBSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type eta(etaSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleL(GA, GB, survivalpSameH, choice_set, choice_equal, nval, phikA, phikB, eta));
    return rcpp_result_gen;
END_RCPP
}
// cartesianProduct
IntegerMatrix cartesianProduct(IntegerVector vec1, IntegerVector vec2);
RcppExport SEXP _FlexRL_cartesianProduct(SEXP vec1SEXP, SEXP vec2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type vec1(vec1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vec2(vec2SEXP);
    rcpp_result_gen = Rcpp::wrap(cartesianProduct(vec1, vec2));
    return rcpp_result_gen;
END_RCPP
}
// ExpandGrid
IntegerMatrix ExpandGrid(IntegerVector vector1, IntegerVector vector2);
RcppExport SEXP _FlexRL_ExpandGrid(SEXP vector1SEXP, SEXP vector2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type vector1(vector1SEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type vector2(vector2SEXP);
    rcpp_result_gen = Rcpp::wrap(ExpandGrid(vector1, vector2));
    return rcpp_result_gen;
END_RCPP
}
// generateSequence
IntegerVector generateSequence(int n);
RcppExport SEXP _FlexRL_generateSequence(SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(generateSequence(n));
    return rcpp_result_gen;
END_RCPP
}
// sampleH
List sampleH(IntegerVector nA, IntegerVector nB, IntegerMatrix links, NumericMatrix survivalpSameH, LogicalVector pivs_stable, List pivsA, List pivsB, IntegerVector nvalues, LogicalVector nonlinkedA, LogicalVector nonlinkedB, List eta, List phi);
RcppExport SEXP _FlexRL_sampleH(SEXP nASEXP, SEXP nBSEXP, SEXP linksSEXP, SEXP survivalpSameHSEXP, SEXP pivs_stableSEXP, SEXP pivsASEXP, SEXP pivsBSEXP, SEXP nvaluesSEXP, SEXP nonlinkedASEXP, SEXP nonlinkedBSEXP, SEXP etaSEXP, SEXP phiSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type nA(nASEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nB(nBSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type links(linksSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type survivalpSameH(survivalpSameHSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type pivs_stable(pivs_stableSEXP);
    Rcpp::traits::input_parameter< List >::type pivsA(pivsASEXP);
    Rcpp::traits::input_parameter< List >::type pivsB(pivsBSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type nvalues(nvaluesSEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type nonlinkedA(nonlinkedASEXP);
    Rcpp::traits::input_parameter< LogicalVector >::type nonlinkedB(nonlinkedBSEXP);
    Rcpp::traits::input_parameter< List >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< List >::type phi(phiSEXP);
    rcpp_result_gen = Rcpp::wrap(sampleH(nA, nB, links, survivalpSameH, pivs_stable, pivsA, pivsB, nvalues, nonlinkedA, nonlinkedB, eta, phi));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_FlexRL_callFunction", (DL_FUNC) &_FlexRL_callFunction, 2},
    {"_FlexRL_F2", (DL_FUNC) &_FlexRL_F2, 2},
    {"_FlexRL_F33", (DL_FUNC) &_FlexRL_F33, 3},
    {"_FlexRL_sspaste2", (DL_FUNC) &_FlexRL_sspaste2, 1},
    {"_FlexRL_initΔMap", (DL_FUNC) &_FlexRL_initΔMap, 0},
    {"_FlexRL_Δfind", (DL_FUNC) &_FlexRL_Δfind, 0},
    {"_FlexRL_sampleD", (DL_FUNC) &_FlexRL_sampleD, 9},
    {"_FlexRL_sampleNL", (DL_FUNC) &_FlexRL_sampleNL, 3},
    {"_FlexRL_sampleL", (DL_FUNC) &_FlexRL_sampleL, 9},
    {"_FlexRL_cartesianProduct", (DL_FUNC) &_FlexRL_cartesianProduct, 2},
    {"_FlexRL_ExpandGrid", (DL_FUNC) &_FlexRL_ExpandGrid, 2},
    {"_FlexRL_generateSequence", (DL_FUNC) &_FlexRL_generateSequence, 1},
    {"_FlexRL_sampleH", (DL_FUNC) &_FlexRL_sampleH, 12},
    {NULL, NULL, 0}
};

RcppExport void R_init_FlexRL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
