#' @title dpareto 
#' @description Continuous Probability density function  of Pereto 
#' @details #Used for modeling noise coverage
#' @param x Outcome
#' @param alpha Model parameter
#' @param x0 Lower limit of distr
#' @param logged Whether logged values should be returned
#' @export 
dpareto = function(x,alpha,x0,logged=FALSE) {
  vals = log(alpha) + alpha*log(x0) - (alpha + 1)*log(x)
  if(!logged) vals = exp(vals)
  return(vals)
}
