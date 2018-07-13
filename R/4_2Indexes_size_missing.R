#' @title Indexing gaps size
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Stores the position of the begining of each gap and their respective size within a multivariate signal.
#' @param data multivariate signal
#' @return returns a list with one element per signal.
#' Within each element of this list, the first column gives the position of the begining of each gap and the second column its size.
#' @import rlist
#' @examples
#' data(dataDTWUMI)
#' id_NA <- Indexes_size_missing_multi(dataDTWUMI$incomplete_signal)

Indexes_size_missing_multi <- function(data){
  store_miss <- list()
  for(icol in 1:ncol(data)){
    pos_w <- c()
    num_w <- c()
    count <- 0
    i <- 1
    while(i < nrow(data)){
      if(is.na(data[i, icol])==TRUE){
        count <- count+1
      }
      if(is.na(data[i+1, icol])==FALSE){
        if(count>0){
          num_w <- c(num_w, count)
          pos_w <- c(pos_w, i-count+1)
        }
        count <- 0
      }
      i <- i+1
    }
    store_w <- matrix(0, ncol=2, nrow=length(pos_w))
    store_w[, 1] <- pos_w
    store_w[, 2] <- num_w
    store_miss <- rlist::list.append(store_miss, store_w)
  }
  return(store_miss)
}