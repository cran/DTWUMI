#' @title Adaptive Feature Based Dynamic Time Warping algorithm for multivariate signals
#' @author DEZECACHE Camille, PHAN Thi Thu Hong, POISSON-CAILLAULT Emilie
#' @description 
#' This function estimates a distance matrix which is used as an input in dtw() function (package dtw) to align two multivariate signals following Adaptative Feature Based Dynamic Time Warping algorithm (AFBDTW).
#' @param q query dataframe
#' @param r reference dataframe
#' @param w1 weight of local feature VS global feature.
#'  By default, w1 = 0.5, and by definition, w2 = 1 - w1.

.dist_afbdtw_matrix <- function(q, r, w1=0.5){
  
  w2 <- 1-w1
  
  if(w1<=0){stop("Weights should be positive")}
  
    ql <- .local_feature_multi(q)
    rl <- .local_feature_multi(r)
    qg <- .global_feature_multi(q)
    rg <- .global_feature_multi(r)
    dist_local <- .dist_matrix_multi(ql, rl)
    dist_global <- .dist_matrix_multi(qg, rg)
    dist <- w1*dist_local + w2*dist_global
    
  out <- dist
}