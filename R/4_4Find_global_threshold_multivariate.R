#' @title Global threshold for missing data imputation in multivariate signals
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Finds a threshold for missing data imputation in multivariate signals.
#' @import dtw
#' @importFrom lsa cosine
#' @keywords internal

.DTW_threshold_global_multivariate <- function(query, database, gap_size, i_start, i_finish, step_threshold, threshold_cos, thresh_cos_stop, ...){
  
  # Estimate global features
  query_FM <- .features_matrix(query)
  
  # Initialization
  Cosine_threshold <- c()
  pos_i <- c()
  threshold_cos_temp <- threshold_cos
  
  while(length(Cosine_threshold)==0&&threshold_cos_temp>thresh_cos_stop){
    i <- i_start
    while(i<=i_finish){ 
      k <- i+gap_size-1
      ref <- database[i:k, ]
      ref_FM <- .features_matrix(ref)
    
      cos_threshold <- cosine(as.numeric(query_FM), as.numeric(ref_FM))
    
      if(cos_threshold[1]>=threshold_cos_temp){
        pos_i <- c(pos_i, i)
        Cosine_threshold <- c(Cosine_threshold, cos_threshold[1])
      }
      i <- i+step_threshold
    }
    threshold_cos_temp <- threshold_cos_temp-0.01
  }
  
  if(length(Cosine_threshold)==0){stop("No similar window looks appropriate for imputation")}

  return(pos_i)
}