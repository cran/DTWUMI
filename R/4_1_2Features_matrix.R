#' @title Feature matrix
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description Store global features of a multivariate signal within a feature matrix.
#' Following .globalfeatures() function, features estimated are:
#'  \itemize{
#'  \item{minx: }{minimum value of the input vector}
#'  \item{maxx: }{maximum value of the input vector}
#'  \item{avg: }{average value of the input vector}
#'  \item{medianx: }{median of the input vector}
#'  \item{std: }{standard deviation of the input vector}
#'  \item{momx: }{skewness}
#'  \item{nop: }{number of peaks of the input vector}
#'  \item{len: }{length of the input vector}
#'  \item{entro: }{measure of entropy}
#'  }
#' @param data multivariate signal
#' @return returns a matrix with one row, each column giving the value of corresponding estimated feature.
#' Values for each signal (column of the dataset storing the multivariate signal) are concatenated.
#' @keywords internal
#' @examples
#' data(dataDTWUMI)
#' gf <- .features_matrix(dataDTWUMI)

.features_matrix <- function(data){
  att <- c()
  for (i in 1:ncol(data)){
    gf <- .globalfeatures(data[, i])
    att <- cbind(att, gf)
  }
  out <- att	
}