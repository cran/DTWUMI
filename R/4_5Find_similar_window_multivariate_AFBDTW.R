#' @title Finding similar windows to a query in multivariate signals
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description This function finds similar windows to a query consisting of a multivariate signal.
#' @import dtw
#' @importFrom lsa cosine
#' @keywords internal

.Finding_similar_window_multivariate_AFBDTW <- function(query, database, threshold2, ...){
  # Initialization
  T <- nrow(query)
  Cost_dist <- c()
  pos_i <- threshold2
  
  for (i in pos_i){
    k <- i+T-1
    ref <- database[i:k, ]
    dist_matrix_AFBDTW <- .dist_afbdtw_matrix(query, ref)
    align <- dtw(dist_matrix_AFBDTW, keep.internals=TRUE, ...)
    cost <- align$normalizedDistance
    Cost_dist <- c(Cost_dist, cost)
  }
  
  min_cost <- min(Cost_dist)
  id_similar <- pos_i[which(Cost_dist==min_cost)]
  
  # id_similar <- id_found_start[which(dist_found==min)]
  # k <- id_similar+T-1
  
  return(list(id_similar))
  # return(list(id_similar, id_found_start))
}