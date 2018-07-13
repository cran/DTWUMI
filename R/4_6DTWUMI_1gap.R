#' @title Imputation of a large gap based on DTW for multivariate signals
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Fills a gap of size 'gap_size' begining at the position 'begin_gap' within a multivariate signal using DTW.
#' @param data a multivariate signals containing gaps
#' @param id_sequence id of the sequence containing the gap to fill (corresponding to the column number)
#' @param begin_gap id of the begining of the gap to fill
#' @param gap_size size of the gap to fill
#' @param DTW_method DTW method used for imputation ("DTW", "DDTW", "AFBDTW"). By default "DTW"
#' @param threshold_cos threshold used to define similar sequences to the query
#' @param thresh_cos_stop Define the lowest cosine threshold acceptable to find a similar window to the query
#' @param step_threshold step used within the loops determining the threshold and the most similar sequence to the query
#' @param ... additional arguments from dtw() function
#' @return returns a list containing the following elements:
#' \itemize{
#'  \item{imputed_values: }{output vector containing the imputation proposal}
#'  \item{id_imputation: }{a vector containing the position of the imputed values extracted}
#'  \item{id_sim_win: }{a vector containing the position of the similar window to the query}
#'  \item{id_gap: }{a vector containing the position gap considered}
#'  \item{id_query: }{a vector containing the position of the query}
#' }
#' @import DTWBI
#' @examples
#' data(dataDTWUMI)
#' dataDTWUMI_gap <- dataDTWUMI[["incomplete_signal"]]
#' t <- 207 ; T <- 40
#' imputation <- DTWUMI_1gap_imputation(dataDTWUMI_gap, id_sequence=1, t, T)
#' plot(dataDTWUMI_gap[, 1], type = "l", lwd = 2)
#' lines(y = imputation$imputed_values, x = imputation$id_gap, col = "red")
#' lines(y = dataDTWUMI_gap[imputation$id_query, 1], x = imputation$id_query, col = "green")
#' lines(y = dataDTWUMI_gap[imputation$id_sim_win, 1], x = imputation$id_sim_win, col = "blue")
#' lines(y = dataDTWUMI_gap[imputation$id_imputation, 1], x = imputation$id_imputation, col = "orange")


DTWUMI_1gap_imputation <- function(data, id_sequence, begin_gap, gap_size, DTW_method="DTW", threshold_cos=0.995, thresh_cos_stop=0.8, step_threshold=2, ...){
  
  if(is.null(ncol(data))){stop("For univariate signal use package DTWBI")}
  
method_options <- c("DTW", "DDTW", "AFBDTW",
                    "dtw", "ddtw", "afbdtw")
if(DTW_method %in% method_options){a <- 1} else {stop("Invalid DTW method")}
  
# Check thresh_cos_stop
if(thresh_cos_stop>=threshold_cos){stop("Parameter thresh_cos_stop is higher than threshold_cos")}
  
# Store location and size of gaps and isolated missing data
store_miss <- Indexes_size_missing_multi(data)

# Create a copy of data
data2 <- data

# Fill isolated NA
for (icol in 1:ncol(data)){
  id_gap1 <- which(store_miss[[icol]][, 2]==1)
  if(length(id_gap1)>0){
    pos1 <- store_miss[[icol]][id_gap1, 1]
    data2[, icol] <- imp_1NA(data[, icol], pos1)
  }
}

# Fill other gaps using trapezoid temporarily
store_miss2 <- Indexes_size_missing_multi(data2) # Id and size of remaining gaps
for (icol in 1:ncol(data)){
  pos_dtw <- store_miss2[[icol]][, 1]
  size_dtw <- store_miss2[[icol]][, 2]
  ma <- max(data[, icol], na.rm=T) 
  if(length(pos_dtw)>0){
    for (id in 1:length(pos_dtw)){
      fuzzy_max <- c()
      size_dtw_id <- size_dtw[id]
      pos <- pos_dtw[id]
      alpha_tra <- c()
      t1 <- pos # Begining of the large gap considered
      t2 <- t1+size_dtw_id-1 # End of the large gap considered
      for(x in t1:t2){alpha_tra <- c(alpha_tra, .trapezoid(x, t1, t2))}
      fuzzy_max <- c(fuzzy_max, ma*alpha_tra)
      data2[t1:t2, icol] <- fuzzy_max
    }
  }
}

#### BEGINING OF IMPUTATION FUNCTION ####

N <- nrow(data2)
  
  # Gap before N/2
  if(begin_gap<N/2){
    pos_start <- begin_gap+gap_size
    Researchbase_a <- data2[pos_start:N, ]
    #####
    # Modify dataset following DTW method used
    if(DTW_method=="DDTW"||DTW_method=="ddtw"){
      Researchbase_a <- .local.derivative.ddtw(Researchbase_a)
      if(length(which(is.na(Researchbase_a)))>0){stop("NA in Researchbase_a2!!!")}
    }
    #####
    query_a <- Researchbase_a[1:gap_size, ]
    i_start <- gap_size+1
    i_finish <- nrow(Researchbase_a)-gap_size
    
    ####
    if(DTW_method=="AFBDTW"||DTW_method=="afbdtw"){
      threshold <- .DTW_threshold_global_multivariate(query_a, Researchbase_a, gap_size, i_start, i_finish, step_threshold, threshold_cos, thresh_cos_stop, ...)
      # STEP 3: find similar windows in the research database		
      id_similar_window <- .Finding_similar_window_multivariate_AFBDTW(query_a, Researchbase_a, gap_size, threshold, ...)
    } else {
      threshold <- .DTW_threshold_global_multivariate(query_a, Researchbase_a, gap_size, i_start, i_finish, step_threshold, threshold_cos, thresh_cos_stop, ...)
      # STEP 3: find similar windows in the research database		
      id_similar_window <- .Finding_similar_window_multivariate(query_a, Researchbase_a, gap_size, threshold, ...)
    }

id_sim_win <- id_similar_window[[1]]

# Imputation
id_similar_finish <- id_sim_win-1
id_similar_start <- id_similar_finish-gap_size+1
imp_values <- data2[pos_start:N, ][id_similar_start:id_similar_finish, id_sequence]
id_imputation <- (pos_start:N)[id_similar_start:id_similar_finish]
id_simwin <- (pos_start:N)[id_sim_win:(id_sim_win+gap_size-1)]
id_query <- (begin_gap+gap_size):(begin_gap+2*gap_size-1)
  }
  
  # Gap after N/2
  if(begin_gap>=N /2){
    Researchbase_b <- data2[1:(begin_gap-1), ]
    #####
    # Modify dataset following DTW method used
    if(DTW_method=="DDTW"||DTW_method=="ddtw"){
      Researchbase_b <- .local.derivative.ddtw(Researchbase_b)
      if(length(which(is.na(Researchbase_b)))>0){stop("NA in Researchbase_b!!!")}
    }
    #####
    start_q <- begin_gap-gap_size
    query_b <- Researchbase_b[start_q:(begin_gap-1), ]
    i_start <- 1
    i_finish <- nrow(Researchbase_b)-2*gap_size
    
    ####
    if(DTW_method=="AFBDTW"||DTW_method=="afbdtw"){
      # Finding a threshold
      threshold <- .DTW_threshold_global_multivariate(query_b, Researchbase_b, gap_size, i_start, i_finish, step_threshold, threshold_cos, thresh_cos_stop, ...)
      # STEP 3: find similar windows in the research database
      id_similar_window <- .Finding_similar_window_multivariate_AFBDTW(query_b, Researchbase_b, gap_size, threshold, ...)
    } else {
      threshold <- .DTW_threshold_global_multivariate(query_b, Researchbase_b, gap_size, i_start, i_finish, step_threshold, threshold_cos, thresh_cos_stop, ...)
      # STEP 3: find similar windows in the research database
      id_similar_window <- .Finding_similar_window_multivariate(query_b, Researchbase_b, gap_size, threshold, ...)
    }

    id_sim_win <- id_similar_window[[1]]
    id_similar_start <- id_sim_win+gap_size
    id_similar_finish <- id_similar_start+gap_size-1
    imp_values <- data2[1:begin_gap, ][id_similar_start:id_similar_finish, id_sequence]
    id_imputation <- (1:begin_gap)[id_similar_start:id_similar_finish]
    id_simwin <- id_sim_win:(id_sim_win+gap_size-1)
    id_query <- (begin_gap-gap_size):(begin_gap-1)
  }
        
  return(list("imputed_values" = imp_values,
              "id_imputation" = id_imputation,
              "id_sim_win" = id_simwin,
              "id_gap" = begin_gap:(begin_gap+gap_size-1),
              "id_query" = id_query))
}