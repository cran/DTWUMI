#' @title Imputing small gaps of size 1 < size < gap_size_threshold
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Imputes small gaps by weighted moving average.
#' @param signal signal
#' @param posMA position of the begining of small gaps
#' @param pMA position of all corresponding values of small gaps
#' @param store_w list storing index and size of gaps
#' @param lag lag of the weighted moving average function
#' @return  returns a vector of imputed values
#' @keywords internal

.imp_NA_WMA = function(signal, posMA, pMA, store_w, lag){
  
  denominator <- lag*(lag+1)/2
  y <- signal
  ya <- signal

  for (i in posMA){
    pos <- store_w[i,1]
    num <- store_w[i,2]
    pos_a <- pos+num-1
    count <- 1
    while(count<num)
    {
      fenetre <- max(pos-lag,1)
      fenetre_a <- min(pos_a+lag, length(signal))
      numerator <- 0
      numerator_a <- 0
      for(j in 1:lag){ 
        if(is.na(y[fenetre+j-1])==TRUE){y[fenetre+j-1] <- 0}
        numerator <- numerator+j*y[fenetre+j-1] 
        if (is.na(ya[fenetre_a-j+1])==TRUE){ya[fenetre_a-j+1] <- 0}
        numerator_a <- numerator_a+j*ya[fenetre_a-j+1]
      }
      y[pos] <- numerator/denominator
      ya[pos_a] <- numerator_a/denominator
      
      pos <- pos+1
      pos_a <- pos_a-1
      count <- count+1
    }
    if(count==num){
      if(pos==length(signal)){y[pos] <- y[pos-1]}
      else{y[pos] <- (y[pos-1]+y[pos+1])/2}
      if(pos_a==1){ya[pos_a] <- ya[pos_a+1]}
      else{ya[pos_a] <- (ya[pos_a-1]+ya[pos_a+1])/2}
    }
  }
  out <- c()
  for(t in pMA){
    tmp <- (y[t]+ya[t])/2
    out <- c(out, tmp)
  }
  return(out)
}