#' @title Trapezoid function
#' @description Used in .SB() function

.trapezoid <- function(t, t1, t2){
  a <- t1 + 0.3*(t2 - t1)
  b <- t1 + 0.7*(t2 - t1)
  if((t<t1)||(t>t2)){out=0}
  else if(t<a){out =(t-t1)/(a-t1)}
  else if(t<=b){out=1}
  else{out=(t-t2)/(b-t2)}
  return(out)
}