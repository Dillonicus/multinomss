#include <iostream>
#include <RcppArmadillo.h>
#include <progress.hpp>
#include <progress_bar.hpp>
using namespace Rcpp;
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppProgress)]]


// Helper and Utility Functions ------------------------------------------------

// The haversine_dist Function calculates the great circle distance between two
// coordinates, which are expressed in latitude/longitude. Based on the
// haversine formula.

NumericVector haversine_dist(const NumericMatrix& locs, NumericVector point) {
  NumericMatrix rads = locs * M_PI / 180;
  
  NumericVector pt = point * M_PI / 180;
  NumericVector d_long = rads(_,0) - pt[0];
  NumericVector d_lat = rads(_,1) - pt[1];
  
  NumericVector a = pow(sin(d_lat / 2.0), 2) + cos(pt[1]) * cos(rads(_,1)) * pow(sin(d_long / 2.0), 2);
  a = ifelse(a < 1, a, 1);
  NumericVector b = 2 * 6378137 * atan(sqrt(a)/sqrt(1-a));
  return b;
}

// k_nearest is used to sort the haversine distances from each point identified
// above in aschending order; used to identify the k nearest points to each
// centroid

IntegerVector k_nearest(NumericVector v, int k) {
  int size = v.size();
  
  IntegerVector idx(size);
  std::iota(idx.begin(), idx.end(), 0);
  
  std::stable_sort(idx.begin(), idx.end(),
                   [&](int i, int j){return v[i] < v[j];});
  
  return idx[Rcpp::Range(0, k-1)];
  
}

// Takes date interval, maximum proportion of study time allowed for each time
// interval, and the total study time as inputs and returns matrix of lower/upper
// bounds for the temporal windows

IntegerMatrix time_windows(const IntegerVector& bounds, const float& time_prop, const int& study_length, const bool& retrospective){
  
  int m = bounds.size();
  float cutoff = time_prop * study_length;
  
  IntegerVector output_l(m*m);
  IntegerVector output_u(m*m);
  
  for(int i = 0; i < m; i++) {
    for(int j = 0; j < m; j++) {
      
      int p = i * m + j;
      output_l[p] = bounds[i];
      output_u[p] = bounds[j];
    }
  }
  
  if(retrospective == TRUE){
    LogicalVector filt = ( output_l < output_u ) & ((output_u - output_l) <= cutoff ) ;
    IntegerVector output_l2 = output_l[ filt ];
    IntegerVector output_u2 = output_u[ filt ];
    IntegerMatrix output = cbind(output_l2, output_u2);
    return output;
  }
  
  else{
    LogicalVector filt = ( output_u == max(output_u) ) & ( output_l < output_u ) & ((output_u - output_l) <= cutoff ) ;
    IntegerVector output_l2 = output_l[ filt ];
    IntegerVector output_u2 = output_u[ filt ];
    IntegerMatrix output = cbind(output_l2, output_u2);
    return output;
  }
  
}

// Similar to the "table" function in R, but requires pre-specified levels

arma::vec counter(const arma::vec& group, const arma::vec& levels){
  arma::vec counts(levels.n_elem, arma::fill::zeros);
  
  int inum = levels.n_rows;
  int jnum = group.n_rows;
  
  for(int i = 0; i < inum; i++){
    int count = 0;
    for(int j = 0; j < jnum; j++){
      if(group(j) == levels(i)){
        count++;
      }
    }
    counts(i) = count;
  }
  return counts;
}

// Spatial and Space-Time Zone Functions ---------------------------------------

//' Calculates the spatial zones, starting with each centroid and increasing size
//' from 1 to pre-specified k points. These represent circular zones in space.
//' 
//' @param locs An n x 2 matrix of centroid point coordinates in latitude/longitude
//' @param k An integer specifying the maximum number of nearest neighbors to include in spatial zones
//' @export
// [[Rcpp::export]]
List spatial_zones(const NumericMatrix& locs, const int& k) {
  int n = locs.nrow();
  List spatial_zones(n*k);
  
  Rcout << "Starting Spatial Zone Calculation..." << endl;
  Progress p(n*k, true);
  
  for(int i = 0; i < n; ++i) {
    NumericVector fullzones = haversine_dist(locs, locs(i, _ ));
    IntegerVector nn = k_nearest(fullzones, k);
    IntegerVector ni = rep(i+1, k);
    IntegerMatrix nns = cbind(ni, nn);
    
    for(int j = 0; j < k; ++j) {
      int m = k * i + j;
      spatial_zones[m] = nns(Rcpp::Range(0, j), _);
      p.increment();
    }
  }
  
  Rcout << "Spatial Zones Done" << endl;
  return spatial_zones;
}

//' Takes spatial zone input and calculates the temporal zones. For each circular
//' zone in space, it effectively varies the height to create cylindrical zones.
//'
//' @param zones A list of spatial zones produced by \code{spatial_zones()} function
//' @param dates A vector containing the dates associated with each observation. If time aggregation is used, this will correspond to time interval that each observation falls into
//' @param time_prop A float/decimal specifying the maximum proportion of the study length to consider as a cluster
//' @param study_length An integer specifying the length of the study in days
//' @param retrospective A logical T/F that specifies whether a retrospective analysis is desired. If \code{retrospective = F} the prospective zones will be calculated
// [[Rcpp::export]]
List temporal_zones(const List& zones, const IntegerVector& dates,
                    const float& time_prop, const int& study_length, const bool& retrospective){
  int m = zones.size();
  List finalzones(m);
  
  Rcout << "Starting Calculation of Temporal Zones..." << endl;
  Progress p(m);
  for(int i = 0; i < m; i++){
    IntegerMatrix subzone = zones[i];
    IntegerVector subzone_ind = subzone(_,1);
    IntegerVector centroid_ind = subzone(_,0);
    IntegerVector subdate = dates[subzone_ind];
    IntegerMatrix time_interval = time_windows(unique(subdate), time_prop, study_length, retrospective);
    
    int n = time_interval.nrow();
    IntegerVector ns = seq(0, subdate.size() - 1);
    List subsubzone(n);
    
    for(int j = 0; j < n; j++){
      IntegerVector time_sub = time_interval(j,_);
      IntegerVector index = ns[ ( subdate > time_sub[0] ) & ( subdate <= time_sub[1] ) ];
      IntegerVector subzone_filt = subzone_ind[index];
      IntegerVector centroid_filt = centroid_ind[index];
      IntegerVector lb(index.size(), time_sub[0]);
      IntegerVector ub(index.size(), time_sub[1]);
      subsubzone[j] = cbind(centroid_filt, lb, ub, subzone_filt);
    }
    finalzones[i] = subsubzone;
    p.increment();
  }
  Rcout << "Temporal Zones Done" << endl;
  return finalzones;
}

//' Performs calculation of both spatial and temporal zones in one step. The arguments are the same as specified in \code{spatial_zones()} and \code{temporal_zones()}
//' @export
// [[Rcpp::export]]
List zones(const NumericMatrix& locs, const int& k, const IntegerVector& dates,
           const float& time_prop, const int& study_length, const bool& retrospective){
  
  List spatialzones = spatial_zones(locs, k);
  List finalzones = temporal_zones(spatialzones, dates, time_prop, study_length, retrospective);
  
  return finalzones;
}

// Multinomial Scan Statistic Functions ----------------------------------------

// multinom_stat_null is used to calculate the total counts and the L0 value for
// the likelihood ratio statistic calculation

List multinom_stat_null(const arma::vec& group, const arma::vec& levels){
  List null(2);
  
  arma::vec totals = counter(group, levels);
  
  arma::mat log_totals = totals % arma::log(totals / sum(totals));
  log_totals.replace(arma::datum::nan,0);
  
  null[0] = totals;
  null[1] = sum(log_totals);
  
  return null;
}

// multinom_stat calculates the test statistic with a matrix of counts and the
// total counts and L0 value (obtained via multinom_stat_null)

arma::mat multinom_stat(const arma::mat& in_zone, const arma::vec& totals,
                        const arma::mat& LL0){
  
  arma::mat logs_in = in_zone % arma::log(in_zone.each_row() / arma::sum(in_zone, 0));
  logs_in.replace(arma::datum::nan, 0);
  arma::mat sum_logs_in = arma::sum(logs_in, 0);
  
  arma::mat out_zone = totals - in_zone.each_col();
  arma::mat logs_out = out_zone % arma::log(out_zone.each_row() / arma::sum(out_zone, 0));
  logs_out.replace(arma::datum::nan,0);
  
  arma::mat LL1_in = arma::sum(logs_in, 0);
  arma::mat LL1_out = arma::sum(logs_out, 0);
  arma::mat LL1 = LL1_in + LL1_out;
  
  return LL1.each_col() - LL0;
}

// [[Rcpp::export]]
arma::vec multinom_scan(const List& zones, const arma::vec& group, const arma::vec& levels){
  int n = zones.size();
  List tot = multinom_stat_null(group, levels)  ;
  arma::vec out(n, arma::fill::zeros);
  
  for(int i = 0; i < n; i++){
    List centroid = zones[i];
    int m = centroid.size();
    arma::mat subout(levels.n_elem, m, arma::fill::zeros);
    
    for(int j = 0; j < m; j++){
      arma::umat group_index = centroid[j];
      arma::vec group_sub = group.elem(group_index.col(3));
      subout.col(j) = counter(group_sub, levels);
    }
    arma::mat test_stat = multinom_stat(subout, tot[0], tot[1]);
    
    if(test_stat.n_elem > 0){
      arma::uword maximum_index = test_stat.index_max();
      double maximum_stat = test_stat(maximum_index);
      out[i] = maximum_stat;
    }
    
    else{
      double maximum_stat = 0;
      out[i] = maximum_stat;
    }
  }
  return out;
}

//' Calculates the most likely cluster for each centroid point by calculating the test statistic for all scanning windows
//' 
//' @param zones list of space-time scanning windows
//' @param group vector corresponding to group membership of each observation
//' @param id vector corresponding to the ID number of the observation
//' @param levels vector denoting the unique levels of the group variable to which an observation can belong
//' @export
// [[Rcpp::export]]
List multinom_mlc(const List& zones, const arma::vec& group, const arma::vec& id, const arma::vec& levels){
  int n = zones.size();
  List tot = multinom_stat_null(group, levels)  ;
  arma::mat mlc(n, 4, arma::fill::zeros);
  List zone_ids(n);
  
  for(int i = 0; i < n; i++){
    List centroid = zones[i];
    int m = centroid.size();
    arma::mat subout(levels.n_elem, m, arma::fill::zeros);
    
    for(int j = 0; j < m; j++){
      arma::umat centroid_sub = centroid[j];
      arma::vec group_sub = group.elem(centroid_sub.col(3));
      subout.col(j) = counter(group_sub, levels);
    }
    
    arma::mat test_stat = multinom_stat(subout, tot[0], tot[1]);
    
    if(test_stat.n_elem > 0){
      arma::uword maximum_index = test_stat.index_max();
      arma::mat max_zone = centroid[maximum_index];
      arma::vec max_zone_ind = max_zone.col(3);
      arma::rowvec sub_max_zone = max_zone(0, arma::span(0, 2));
      arma::vec max_stat(1, arma::fill::zeros);
      max_stat.fill(test_stat(maximum_index));
      
      arma::uvec max_zone_ind_u = arma::conv_to<arma::uvec>::from(max_zone_ind);
      
      sub_max_zone.insert_cols(3, max_stat);
      mlc.row(i) = sub_max_zone;
      arma::vec ids = id.rows(max_zone_ind_u);
      
      zone_ids[i] = ids;
    }
    
    else{
      arma::mat max_zone(1, 4, arma::fill::zeros);
      mlc.row(i) = max_zone.row(0);
      zone_ids[i] = max_zone.col(3);
    }
    
  }
  
  arma::uvec index_out = arma::find(mlc.col(3) > 0);
  arma::mat mlc_out = mlc.rows(index_out);
  zone_ids = zone_ids[as<IntegerVector>(wrap(index_out))];
  
  List out = List::create(mlc_out, zone_ids);
  
  return out;
}

//' Conducts the permutations needed for Monte Carlo hypothesis testing
//' 
//' @param zones list of space-time scanning windows
//' @param group vector corresponding to group membership of each observation
//' @param n_perm number of Monte Carlo permutations to perform
//' @export
// [[Rcpp::export]]
arma::vec multinom_permutation(const List& zones, const arma::vec& group, const arma::vec& levels, const int& n_perm){
  
  arma::vec LLs(n_perm + 1, arma::fill::zeros);
  
  for(int i = 0; i < n_perm; i++){
    Rcout << "Monte Carlo Permutation: " << i << "/" << n_perm << endl;
    arma::vec perm_group = arma::shuffle(group);
    arma::vec perm_stat = multinom_scan(zones, perm_group, levels);
    LLs.row(i) = arma::max(perm_stat);
  }
  return LLs;
}

//' Calculates the p-values from the Monte Carlo permutation-generated test statistics
//' 
//' @param stats vector of test statistics from the data for which p-values are desired
//' @param perm vector of permutation test statistics used to compute p-values
//' @export
// [[Rcpp::export]]
arma::mat p_val(const arma::vec& stats, const arma::vec& perm){
  
  arma::vec a = perm;
  int m = stats.n_elem;
  int n = perm.n_elem;
  
  arma::vec p(m);
  arma::mat out(m, 2);
  
  for(int i = 0; i < m; i ++){
    
    a.row(n-1) = stats.row(i);
    
    arma::vec rank = arma::conv_to< arma::vec>::from(arma::sort_index(arma::sort_index(a, "descend")));
    
    p.row(i) = rank.row(n-1) + 1;
  }
  
  out.col(0) = stats;
  out.col(1) = p/n;
  
  return out;
}