#' @title logLiki
#' @author Oyvind Bleka
#' @description logLiki returns the likelihood for each markers (based on the MLE results)
#'
#' @param mlefit Fitted object using inferEvidence
#' @return ret A vector with log-likelihood-values for each locus for given model
#' @export

#Get log likelihood value for each marker given model fit.
logLiki <- function(mlefit){
  par = mlefit$fit$par
  c <- mlefit$prepareC #returned from prepareC
  locs = c$locNames #loci to evaluate
  
  logLikv = setNames(rep(NA,length(locs)),locs) #obtain loglik per marker
  if(any(is.na(unlist(par)))) return(logLikv) #return if params not found
  
  #Obtain  L(E|g,thetahat) for each marker
  outcomeList = calcLogLikC_prediction(par,c,TRUE) #returnOutcome (all combinations are there)
  
  #Summing over all outcome to get per-marker logLik
  logLikv = log(sapply(outcomeList,function(x) sum(x[,ncol(x)])))
  return( logLikv ) #simply returns  
}



