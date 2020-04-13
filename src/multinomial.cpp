#include <RcppArmadillo.h>
#include <iostream>
#include <progress.hpp>
#include <progress_bar.hpp>

// [[Rcpp::plugins("cpp11")]]

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]

// Helper and Utility Functions ------------------------------------------------

// The haversine_dist Function calculates the great circle distance between two
// coordinates, which are expressed in latitude/longitude. Based on the
// haversine formula.

Rcpp::NumericVector haversine_dist(Rcpp::NumericMatrix locs, Rcpp::NumericVector point) {
  Rcpp::NumericMatrix rads   = locs * M_PI / 180;

  Rcpp::NumericVector pt     = point * M_PI / 180;
  Rcpp::NumericVector d_long = rads(Rcpp::_,0) - pt[0];
  Rcpp::NumericVector d_lat  = rads(Rcpp::_,1) - pt[1];

  Rcpp::NumericVector a      = pow(sin(d_lat / 2.0), 2) + cos(pt[1]) * cos(rads(Rcpp::_,1)) * pow(sin(d_long / 2.0), 2);
  a                          = ifelse(a < 1, a, 1);
  
  Rcpp::NumericVector b      = 2 * 6378137 * atan(sqrt(a)/sqrt(1-a));
  
  return b;
}

arma::vec haversine_dist2(arma::mat locs, arma::rowvec point){
  arma::mat rads = locs * arma::datum::pi / 180;
  arma::rowvec pt = point * arma::datum::pi / 180;

  double pt_long = arma::as_scalar(pt.col(0));
  double pt_lat = arma::as_scalar(pt.col(1));

  arma::vec d_long = rads.col(0) - pt_long;
  arma::vec d_lat = rads.col(1) - pt_lat;

  arma::vec a = pow(sin(d_lat / 2), 2) + cos(rads.col(1)) * cos(pt_lat) % pow(sin(d_long / 2), 2);
  a.elem( find(a >= 1) ).ones();

  a = 2 * 6378137 * atan(sqrt(a)/sqrt(1-a));
  return a;
}

// k_nearest is used to sort the haversine distances from each point identified
// above in aschending order; used to identify the k nearest points to each
// centroid

Rcpp::IntegerVector k_nearest(Rcpp::NumericVector v, int k) {
  int size = v.size();

  Rcpp::IntegerVector idx(size);
  std::iota(idx.begin(), idx.end(), 0);

  std::stable_sort(idx.begin(), idx.end(),
                   [&](int i, int j){return v[i] < v[j];});

  return idx[Rcpp::Range(0, k-1)];

}

arma::uvec k_nearest2(arma::vec v, int k){
  arma::uvec ind = arma::sort_index(v);
  return ind(arma::span(0,k-1));
}

// Takes date interval, maximum proportion of study time allowed for each time
// interval, and the total study time as inputs and returns matrix of lower/upper
// bounds for the temporal windows

Rcpp::IntegerMatrix time_windows(Rcpp::IntegerVector bounds, float time_prop, int study_length){

  int m = bounds.size();
  float cutoff = time_prop * study_length;

  Rcpp::IntegerVector output_l(m*m);
  Rcpp::IntegerVector output_u(m*m);

  for(int i = 0; i < m; i++) {
    for(int j = 0; j < m; j++) {

      int p = i * m + j;
      output_l[p] = bounds[i];
      output_u[p] = bounds[j];
    }
  }

  Rcpp::IntegerVector output_l2 = output_l[ (output_l < output_u) & ((output_u - output_l) <= cutoff) ];
  Rcpp::IntegerVector output_u2 = output_u[ (output_l < output_u) & ((output_u - output_l) <= cutoff) ];
  Rcpp::IntegerMatrix output = Rcpp::cbind(output_l2, output_u2);

  return output;
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

// Calculates the spatial zones, starting with each centroid and increasing size
// from 1 to pre-specified k points. These represent circular zones in space.

// [[Rcpp::export]]
Rcpp::List spatial_zones(const Rcpp::NumericMatrix& locs, const int& k) {
  int n = locs.nrow();
  Rcpp::List spatial_zones(n*k);

  std::cout << "Starting Spatial Zone Calculation..." << std::endl;
  Progress p(n*k, true);

  for(int i = 0; i < n; ++i) {
    Rcpp::NumericVector fullzones = haversine_dist(locs, locs(i, Rcpp::_ ));
    Rcpp::IntegerVector nn = k_nearest(fullzones, k);
    Rcpp::IntegerVector ni = Rcpp::rep(i+1, k);
    Rcpp::IntegerMatrix nns = Rcpp::cbind(ni, nn);

    for(int j = 0; j < k; ++j) {
      int m = k * i + j;
      spatial_zones[m] = nns(Rcpp::Range(0, j), Rcpp::_);
      p.increment();
    }
  }

  std::cout << "Spatial Zones Done" << std::endl;
  return spatial_zones;
}

// Takes spatial zone input and calculates the temporal zones. For each circular
// zone in space, it effectively varies the height to create cylindrical zones.

// [[Rcpp::export]]

Rcpp::List temporal_zones(const Rcpp::List& centroid_zones, const Rcpp::IntegerVector& dates,
                    const float& time_prop, const int& max_time){
  int m = centroid_zones.size();
  Rcpp::List finalzones(m);

  std::cout << "Starting Calculation of Temporal Zones..." << std::endl;
  Progress p(m);
  for(int i = 0; i < m; i++){
    Rcpp::IntegerMatrix subzone = centroid_zones[i];
    Rcpp::IntegerVector subzone_ind = subzone(Rcpp::_,1);
    Rcpp::IntegerVector centroid_ind = subzone(Rcpp::_,0);
    Rcpp::IntegerVector subdate = dates[subzone_ind];
    Rcpp::IntegerMatrix time_interval = time_windows(unique(subdate), time_prop, max_time);

    int n = time_interval.nrow();
    Rcpp::IntegerVector ns = Rcpp::seq(0, subdate.size() - 1);
    Rcpp::List subsubzone(n);

    for(int j = 0; j < n; j++){
      Rcpp::IntegerVector time_sub = time_interval(j,Rcpp::_);
      Rcpp::IntegerVector index = ns[ ( subdate > time_sub[0] ) & ( subdate <= time_sub[1] ) ];
      Rcpp::IntegerVector subzone_filt = subzone_ind[index];
      Rcpp::IntegerVector centroid_filt = centroid_ind[index];
      Rcpp::IntegerVector lb(index.size(), time_sub[0]);
      Rcpp::IntegerVector ub(index.size(), time_sub[1]);
      subsubzone[j] = Rcpp::cbind(centroid_filt, lb, ub, subzone_filt);
    }
    finalzones[i] = subsubzone;
    p.increment();
  }
  std::cout << "Temporal Zones Done" << std::endl;
  return finalzones;
}

// One function that ties the spatial and temporal functions together

// [[Rcpp::export]]
Rcpp::List zones(const Rcpp::NumericMatrix& locs, const int& k, const Rcpp::IntegerVector& dates,
           const float& time_prop, const int& max_time){

  Rcpp::List spatialzones = spatial_zones(locs, k);
  Rcpp::List finalzones = temporal_zones(spatialzones, dates, time_prop, max_time);

  return finalzones;
}

// Multinomial Scan Statistic Functions ----------------------------------------

// multinom_stat_null is used to calculate the total counts and the L0 value for
// the likelihood ratio statistic calculation

Rcpp::List multinom_stat_null(const arma::vec& group, const arma::vec& levels){
  Rcpp::List nl(2);

  arma::vec totals = counter(group, levels);

  arma::mat log_totals = totals % arma::log(totals / sum(totals));
  log_totals.replace(arma::datum::nan,0);

  nl[0] = totals;
  nl[1] = sum(log_totals);

  return nl;
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
arma::vec multinom_scan(const Rcpp::List& zones, const arma::vec& group, const arma::vec& levels){
  int n = zones.size();
  Rcpp::List tot = multinom_stat_null(group, levels)  ;
  arma::vec out(n, arma::fill::zeros);

  for(int i = 0; i < n; i++){
    Rcpp::List centroid = zones[i];
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

// [[Rcpp::export]]
Rcpp::List multinom_mlc(const Rcpp::List& zones, const arma::vec& group, const arma::vec& id, const arma::vec& levels){
  int n = zones.size();
  Rcpp::List tot = multinom_stat_null(group, levels);
  arma::mat mlc(n, 4, arma::fill::zeros);
  Rcpp::List zone_ids(n);
  
  for(int i = 0; i < n; i++){
    Rcpp::List centroid = zones[i];
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
  zone_ids = zone_ids[Rcpp::as<Rcpp::IntegerVector>(Rcpp::wrap(index_out))];
  
  Rcpp::List out = Rcpp::List::create(mlc_out, zone_ids);
  
  return out;
}

// [[Rcpp::export]]
arma::vec multinom_permutation(const Rcpp::List& zones, const arma::vec& group, const arma::vec& levels, int n_perm){

  arma::vec LLs(n_perm + 1, arma::fill::zeros);

  //Progress p(n_perm, true);
  for(int i = 0; i < n_perm; i++){
    std::cout.clear();
    std::cout << "Monte Carlo Permutation: " << i << "/" << n_perm << std::endl;
    arma::vec perm_group = arma::shuffle(group);
    arma::vec perm_stat = multinom_scan(zones, perm_group, levels);
    LLs.row(i) = arma::max(perm_stat);
    //p.increment();
  }
  return LLs;
}

// [[Rcpp::export]]

arma::mat p_val(const arma::vec& stats, const arma::vec& perm){
  
  arma::vec a = perm;
  
  int m = stats.n_elem;
  arma::vec p(m);
  arma::mat out(m, 2);
  
  for(int i = 0; i < m; i ++){
    
    a.row(999) = stats.row(i);

    arma::vec rank = arma::conv_to< arma::vec>::from(arma::sort_index(arma::sort_index(a, "descend")));

    p.row(i) = rank.row(999) + 1;
  }
  
  out.col(0) = stats;
  out.col(1) = p/1000;
  
  return out;
}
