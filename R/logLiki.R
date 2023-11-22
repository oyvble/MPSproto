#' @title logLiki
#' @author Oyvind Bleka
#' @description logLiki returns the likelihood for each markers (based on the MLE results)
#'
#' @param mlefit Fitted object using inferEvidence
#' @return A vector with log-likelihood-values for each locus for given model
#' @export

#Get log likelihood value for each marker given model fit.
logLiki <- function(mlefit){
  par = mlefit$fit$par
  c <- mlefit$prepareC #returned from prepareC

  #Obtain  L(E|g,thetahat) for each marker
  logLiki = calcLogLikC_prediction(par,c,returnPerMarker=TRUE) #return logLik per marker
  
  return( logLiki ) 
}



