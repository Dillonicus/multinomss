#' Used to aggregate time by days, weeks, months or years
#' 
#'  @param date A vector of dates to be aggregated
#'  @param agg_unit Determines what time unit to aggregate across
#'  @param agg_interval Determines how many time units to aggregate across 
#'  @export
time_intervals <- function(date, agg_unit = c('days', 'weeks', 'months', 'years'), agg_interval){
  step <- paste(paste0('-', agg_interval), agg_unit)
  breaks <- rev(seq(max(date), min(date), by = step))
  if (breaks[1] != min(date)) {
    breaks <- c(min(date), breaks)
  }
  
  as.Date(cut(date, breaks, right = F, include.lowest = T))
  
}

#' Used to format results for output
#' 
#' @param mlcs A list produced by the \code{multinom_mlc} function
#' @param perm_p Permutation p-values produced by the \code{multinom_permutation} function
#' @param ids A vector that corresponds to the ID of each observation
#' @param locs A matrix of coordinates for each observation/centroid
#' @export
results <- function (mlcs, perm_p, ids, locs) {
  mlcs[[1]] <- data.table::data.table(mlcs[[1]])
  
  data.table::setnames(
    mlcs[[1]], 
    old = names(mlcs[[1]]), 
    new = c("Centroid", "Time Lower", "Time Upper", "LR Stat"))
  
  
    ps <- data.table::data.table(
      multinomss::p_val(
        unique(mlcs[[1]][["LR Stat"]]),
        perm_p))
    
    data.table::setnames(
      ps, 
      old = names(ps), 
      new = c("LR Stat", "p-value"))
    
    mlcs[[1]] <- data.table::merge.data.table(
      mlcs[[1]], 
      ps,
      all = T, 
      sort = F)[, c(2, 3, 4, 1, 5)]
    
    mlcs <- purrr::reduce(mlcs, cbind)
    
    mlcs <- mlcs[order(`LR Stat`, decreasing = T)
                 ][, .SD[1], by = .(Centroid)]
    
    data.table::setnames(mlcs, old = c("elt"), new = c("IDs"))
    
    mlcs <- cbind(mlcs, locs[mlcs$Centroid + 1,])
    
    mlcs[, `:=`(`Time Lower` = as.Date(`Time Lower`, "1970-01-01"),
                `Time Upper` = as.Date(`Time Upper`, "1970-01-01"),
                Centroid = ids[Centroid+1],
                Radius = multinomss::max_dist(IDs, ids, locs))]
    
    # if(non_overlapping == TRUE){
    #   select <- multinomss::non_overlap(mlcs[, IDs])
    #   mlcs <- mlcs[as.vector(select),]
    #   }
    mlcs[]
    }


#' Produces an interactive map to visualize the results of an analysis
#'
#' @param results data.table containing the formatted results of the analysis.
#' @param mlc_only Boolean value. If true, only map the most likely cluster.
#' @param map_layer Specifies the desired base tile layer for the map. Available layers can be previewed at https://leaflet-extras.github.io/leaflet-providers/preview/.
#' @export
cluster_map <- function(results, mlc_only = TRUE, map_layer = "CartoDB.Voyager"){
  
  if(mlc_only == TRUE){
    dat <- sf::st_as_sf(results[1,], coords = c("long", "lat"))
  }
  else{
    dat <- sf::st_as_sf(results, coords = c("long", "lat"))
  }
  
  
  leaflet::leaflet() %>%
    leaflet::addProviderTiles(provider = leaflet::providers[[map_layer]]) %>%
    leaflet::addCircleMarkers(data = dat, 
                   color = 'red', 
                   fill = 'red', 
                   fillOpacity = 0,
                   radius = 1) %>%
    leaflet::addCircles(data = dat, 
               radius = ~ I(Radius * 1000), 
               weight = 1, 
               color = 'blue', 
               fillOpacity = .2)  
}


