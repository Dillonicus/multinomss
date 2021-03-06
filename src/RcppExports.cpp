// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// spatial_zones
List spatial_zones(const NumericMatrix& locs, const float& spatial_prop);
RcppExport SEXP _multinomss_spatial_zones(SEXP locsSEXP, SEXP spatial_propSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< const float& >::type spatial_prop(spatial_propSEXP);
    rcpp_result_gen = Rcpp::wrap(spatial_zones(locs, spatial_prop));
    return rcpp_result_gen;
END_RCPP
}
// temporal_zones
List temporal_zones(const List& zones, const IntegerVector& dates, const float& time_prop, const int& study_length);
RcppExport SEXP _multinomss_temporal_zones(SEXP zonesSEXP, SEXP datesSEXP, SEXP time_propSEXP, SEXP study_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type zones(zonesSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type dates(datesSEXP);
    Rcpp::traits::input_parameter< const float& >::type time_prop(time_propSEXP);
    Rcpp::traits::input_parameter< const int& >::type study_length(study_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(temporal_zones(zones, dates, time_prop, study_length));
    return rcpp_result_gen;
END_RCPP
}
// zones
List zones(const NumericMatrix& locs, const float& spatial_prop, const IntegerVector& dates, const float& time_prop, const int& study_length);
RcppExport SEXP _multinomss_zones(SEXP locsSEXP, SEXP spatial_propSEXP, SEXP datesSEXP, SEXP time_propSEXP, SEXP study_lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< const float& >::type spatial_prop(spatial_propSEXP);
    Rcpp::traits::input_parameter< const IntegerVector& >::type dates(datesSEXP);
    Rcpp::traits::input_parameter< const float& >::type time_prop(time_propSEXP);
    Rcpp::traits::input_parameter< const int& >::type study_length(study_lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(zones(locs, spatial_prop, dates, time_prop, study_length));
    return rcpp_result_gen;
END_RCPP
}
// max_dist
NumericVector max_dist(List idlist, arma::vec idvec, arma::mat locs);
RcppExport SEXP _multinomss_max_dist(SEXP idlistSEXP, SEXP idvecSEXP, SEXP locsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type idlist(idlistSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type idvec(idvecSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type locs(locsSEXP);
    rcpp_result_gen = Rcpp::wrap(max_dist(idlist, idvec, locs));
    return rcpp_result_gen;
END_RCPP
}
// multinom_mlc
List multinom_mlc(const List& zones, const arma::vec& group, const arma::vec& id);
RcppExport SEXP _multinomss_multinom_mlc(SEXP zonesSEXP, SEXP groupSEXP, SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type zones(zonesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(multinom_mlc(zones, group, id));
    return rcpp_result_gen;
END_RCPP
}
// multinom_permutation
arma::vec multinom_permutation(const List& zones, const arma::vec& group, const int& n_perm);
RcppExport SEXP _multinomss_multinom_permutation(SEXP zonesSEXP, SEXP groupSEXP, SEXP n_permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type zones(zonesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const int& >::type n_perm(n_permSEXP);
    rcpp_result_gen = Rcpp::wrap(multinom_permutation(zones, group, n_perm));
    return rcpp_result_gen;
END_RCPP
}
// p_val
arma::mat p_val(const arma::vec& stats, const arma::vec& perm);
RcppExport SEXP _multinomss_p_val(SEXP statsSEXP, SEXP permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type stats(statsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type perm(permSEXP);
    rcpp_result_gen = Rcpp::wrap(p_val(stats, perm));
    return rcpp_result_gen;
END_RCPP
}
// non_overlap
arma::uvec non_overlap(List ids);
RcppExport SEXP _multinomss_non_overlap(SEXP idsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type ids(idsSEXP);
    rcpp_result_gen = Rcpp::wrap(non_overlap(ids));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_multinomss_spatial_zones", (DL_FUNC) &_multinomss_spatial_zones, 2},
    {"_multinomss_temporal_zones", (DL_FUNC) &_multinomss_temporal_zones, 4},
    {"_multinomss_zones", (DL_FUNC) &_multinomss_zones, 5},
    {"_multinomss_max_dist", (DL_FUNC) &_multinomss_max_dist, 3},
    {"_multinomss_multinom_mlc", (DL_FUNC) &_multinomss_multinom_mlc, 3},
    {"_multinomss_multinom_permutation", (DL_FUNC) &_multinomss_multinom_permutation, 3},
    {"_multinomss_p_val", (DL_FUNC) &_multinomss_p_val, 2},
    {"_multinomss_non_overlap", (DL_FUNC) &_multinomss_non_overlap, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_multinomss(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
