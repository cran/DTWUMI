#' @title Large gaps imputation based on DTW for multivariate signals
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Fills all gaps within a multivariate signal. Gaps of size 1 are filled using the average values of nearest neighbours.
#' Gaps of size >1 and <gap_size_threshold are filled using weighted moving average. Larger gaps are filled using DTW.
#' @param data a multivariate signals containing gaps
#' @param gap_size_threshold threshold above which dtw based imputation is computed. Below this threshold, a weighted moving average is calculated
#' @param DTW_method DTW method used for imputation ("DTW", "DDTW", "AFBDTW"). By default "DTW"
#' @param threshold_cos threshold used to define similar sequences to the query
#' @param thresh_cos_stop Define the lowest cosine threshold acceptable to find a similar window to the query
#' @param step_threshold step used within the loops determining the threshold and the most similar sequence to the query
#' @param ... additional arguments from dtw() function
#' @return returns a list containing a dataframe of completed signals
#' @import DTWBI
#' @examples
#' data(dataDTWUMI)
#' dataDTWUMI_gap <- dataDTWUMI[["incomplete_signal"]]
#' imputation <- DTWUMI_imputation(dataDTWUMI_gap, gap_size_threshold = 10)
#' plot(dataDTWUMI_gap[, 1], type = "l", lwd = 2)
#' lines(imputation$output[, 1], col = "red")
#' plot(dataDTWUMI_gap[, 2], type = "l", lwd = 2)
#' lines(imputation$output[, 2], col = "red")
#' plot(dataDTWUMI_gap[, 3], type = "l", lwd = 2)
#' lines(imputation$output[, 3], col = "red")

DTWUMI_imputation <- function(data, gap_size_threshold, DTW_method="DTW", threshold_cos=0.995, thresh_cos_stop=0.8, step_threshold=2, ...){
  gap_size_threshold <- gap_size_threshold
  DTW_method <- DTW_method
  threshold_cos <- threshold_cos
  thresh_cos_stop <- thresh_cos_stop
  step_threshold <- step_threshold

method_options <- c("DTW", "DDTW", "AFBDTW",
                      "dtw", "ddtw", "afbdtw")
if(DTW_method %in% method_options){
  print(DTW_method)} else {stop("Invalid DTW method")}
  
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

# Fill small gaps > 1 and < gap_size_threshold using weighted moving average
for (icol in 1:ncol(data)){
  posMA <- which(store_miss[[icol]][, 2]>1&store_miss[[icol]][, 2]<gap_size_threshold)
  if(length(posMA)>0){
    pMA <- c()
    for(i in posMA){
      count <- 0
      num <- store_miss[[icol]][i, 2]
      while(count<num){ # Create a sequence of missing values which will be imputed using WMA
        (j <- store_miss[[icol]][i, 1] + count)
        (pMA <- c(pMA, j))
        (count <- count+1)
      }
    }
    signal <- data[, icol]
    store_w <- store_miss[[icol]]
    imp_valueMA <- .imp_NA_WMA(data[, icol], posMA, pMA, store_w, lag = round(mean(store_miss[[icol]][, 2])))
    data2[pMA, icol] <- imp_valueMA
  }
}

#### Imputation of large gaps using DTW
store_miss2 <- Indexes_size_missing_multi(data2)

data_final <- data2

  for (icol_general in 1:ncol(data)){
  gaps_number <- nrow(store_miss2[[icol_general]])
      if(gaps_number>0){
      for (id in 1:gaps_number){
        begin <- store_miss[[icol_general]][id, 1]
        size <- store_miss[[icol_general]][id, 2]
        impute_temp <- DTWUMI_1gap_imputation(data=data2, id_sequence=icol_general, begin_gap=begin, gap_size=size,
                                              threshold_cos=threshold_cos, thresh_cos_stop=thresh_cos_stop, step_threshold=step_threshold, ...)
        data_final[begin:(begin+size-1), icol_general] <- impute_temp$imputed_values
        }
    }
  }
return(list("output" = data_final))
}