// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>

using namespace Rcpp;

// spatial_zones
Rcpp::List spatial_zones(const Rcpp::NumericMatrix& locs, const int& k);
RcppExport SEXP _multinomss_spatial_zones(SEXP locsSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    rcpp_result_gen = Rcpp::wrap(spatial_zones(locs, k));
    return rcpp_result_gen;
END_RCPP
}
// temporal_zones
Rcpp::List temporal_zones(const Rcpp::List& centroid_zones, const Rcpp::IntegerVector& dates, const float& time_prop, const int& max_time);
RcppExport SEXP _multinomss_temporal_zones(SEXP centroid_zonesSEXP, SEXP datesSEXP, SEXP time_propSEXP, SEXP max_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type centroid_zones(centroid_zonesSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type dates(datesSEXP);
    Rcpp::traits::input_parameter< const float& >::type time_prop(time_propSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_time(max_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(temporal_zones(centroid_zones, dates, time_prop, max_time));
    return rcpp_result_gen;
END_RCPP
}
// zones
Rcpp::List zones(const Rcpp::NumericMatrix& locs, const int& k, const Rcpp::IntegerVector& dates, const float& time_prop, const int& max_time);
RcppExport SEXP _multinomss_zones(SEXP locsSEXP, SEXP kSEXP, SEXP datesSEXP, SEXP time_propSEXP, SEXP max_timeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::NumericMatrix& >::type locs(locsSEXP);
    Rcpp::traits::input_parameter< const int& >::type k(kSEXP);
    Rcpp::traits::input_parameter< const Rcpp::IntegerVector& >::type dates(datesSEXP);
    Rcpp::traits::input_parameter< const float& >::type time_prop(time_propSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_time(max_timeSEXP);
    rcpp_result_gen = Rcpp::wrap(zones(locs, k, dates, time_prop, max_time));
    return rcpp_result_gen;
END_RCPP
}
// multinom_scan
arma::vec multinom_scan(const Rcpp::List& zones, const arma::vec& group, const arma::vec& levels);
RcppExport SEXP _multinomss_multinom_scan(SEXP zonesSEXP, SEXP groupSEXP, SEXP levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type zones(zonesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type levels(levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(multinom_scan(zones, group, levels));
    return rcpp_result_gen;
END_RCPP
}
// multinom_mlc
Rcpp::List multinom_mlc(const Rcpp::List& zones, const arma::vec& group, const arma::vec& id, const arma::vec& levels);
RcppExport SEXP _multinomss_multinom_mlc(SEXP zonesSEXP, SEXP groupSEXP, SEXP idSEXP, SEXP levelsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type zones(zonesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type id(idSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type levels(levelsSEXP);
    rcpp_result_gen = Rcpp::wrap(multinom_mlc(zones, group, id, levels));
    return rcpp_result_gen;
END_RCPP
}
// multinom_permutation
arma::vec multinom_permutation(const Rcpp::List& zones, const arma::vec& group, const arma::vec& levels, int n_perm);
RcppExport SEXP _multinomss_multinom_permutation(SEXP zonesSEXP, SEXP groupSEXP, SEXP levelsSEXP, SEXP n_permSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List& >::type zones(zonesSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type group(groupSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type levels(levelsSEXP);
    Rcpp::traits::input_parameter< int >::type n_perm(n_permSEXP);
    rcpp_result_gen = Rcpp::wrap(multinom_permutation(zones, group, levels, n_perm));
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

static const R_CallMethodDef CallEntries[] = {
    {"_multinomss_spatial_zones", (DL_FUNC) &_multinomss_spatial_zones, 2},
    {"_multinomss_temporal_zones", (DL_FUNC) &_multinomss_temporal_zones, 4},
    {"_multinomss_zones", (DL_FUNC) &_multinomss_zones, 5},
    {"_multinomss_multinom_scan", (DL_FUNC) &_multinomss_multinom_scan, 3},
    {"_multinomss_multinom_mlc", (DL_FUNC) &_multinomss_multinom_mlc, 4},
    {"_multinomss_multinom_permutation", (DL_FUNC) &_multinomss_multinom_permutation, 4},
    {"_multinomss_p_val", (DL_FUNC) &_multinomss_p_val, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_multinomss(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}