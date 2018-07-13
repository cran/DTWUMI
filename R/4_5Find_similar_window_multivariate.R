#' @title Finding similar windows to a query in multivariate signals
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description This function finds similar windows to a query consisting of a multivariate signal.
#' @import dtw
#' @importFrom lsa cosine
#' @keywords internal

.Finding_similar_window_multivariate <- function(query, database, gap_size, threshold2, ...){
  
   # Initialization
   Cost_dist <- c()
   pos_i <- threshold2
  
  for (i in pos_i){
    k <- i+gap_size-1
    ref <- database[i:k, ]
    align <- dtw(query, ref, keep.internals=TRUE, dist.method="Euclidean")
    cost <- align$normalizedDistance
    Cost_dist <- c(Cost_dist, cost)
  }
   
   min_cost <- min(Cost_dist)
   id_similar <- pos_i[which(Cost_dist==min_cost)]

   return(list(id_similar))
}