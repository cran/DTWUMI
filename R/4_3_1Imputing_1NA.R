#' @title Imputing gaps of size 1
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Imputes isolated missing values based on the average of nearest neighbours.
#' @param data a univariate signal
#' @param pos1 the position of the begining of gaps of size 1, as obtained using Indexes_size_missing_multi() function
#' @return returns a new vector of same size with imputed values

imp_1NA<-function(data, pos1){
  N <- length(data)
  begin <- sapply(pos1, function(x){max(1, x-1)})
  end <- sapply(pos1, function(x){min(N, (x+1))})
  imp_value1 <- rowMeans(cbind(data[begin], data[end]), na.rm=TRUE)
  data_out <- data
  data_out[pos1] <- imp_value1
  return(data_out)
}