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
#' @export
results <- function (mlcs, perm_p = NULL, ids) 
{
  mlcs[[1]] <- data.table::data.table(mlcs[[1]])
  data.table::setnames(mlcs[[1]], old = names(mlcs[[1]]), new = c("Centroid", 
                                                                  "Time Lower", 
                                                                  "Time Upper", 
                                                                  "LR Stat"))
  
  
  if(!is.null(perm_p)){
    ps <- data.table::data.table(multinomss::p_val(unique(mlcs[[1]][["LR Stat"]]), 
                                                   perm_p))
    
    data.table::setnames(ps, old = names(ps), new = c("LR Stat", 
                                                      "p-value"))
    
    mlcs[[1]] <- data.table::merge.data.table(mlcs[[1]], ps, 
                                              all = T, sort = F)[, c(2, 3, 4, 1, 5)]
    
    mlcs <- purrr::reduce(mlcs, cbind)
  }
  
  else{
    mlcs[[1]] <- mlcs[[1]][, c(2, 3, 4, 1)]
    mlcs <- purrr::reduce(mlcs, cbind)
  }
  
  mlcs <- mlcs[order(`LR Stat`, decreasing = T)][, .SD[1], 
                                                 by = .(Centroid)]
  data.table::setnames(mlcs, old = c("elt"), new = c("IDs"))
  mlcs[, `:=`(Centroid, ids[Centroid])]
  mlcs[]
}


