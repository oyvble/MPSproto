#' @title prefit_prediction
#' @description Fit parameter start values for MLE
#' @details
#' Notice that the marker efficiency is scaled with the PHs in order to be taken into account.
#' @param c The object returned from prepareC_prediction
#' @return Estimated parameters (mean,coefVar,slopePara) per replicate

#Helpfunction to fit param start-values (Marker amplification taken into account)
prefit_prediction = function(c) {
  
  #Prefitting data based on the model for sum of the peak heights  to find proper startvalues for MLE 
  nLocs = c$nLocs #number of markers
  nReps = c$nReps #obtain number of replicates
  sumY <- meanbp <- matrix(NA,nrow=nReps,ncol=nLocs) #sum PH, bp per marker per replicates
  for(m in seq_len(nLocs)) { #for each marker
    markerEff = c$markerEfficiency[m]
    indA = (0:(c$nAlleles[m]-1)) + c$startIndMarker_nAlleles[m] #obtain indices for alleles (not replicates)
    for(r in seq_len(nReps)) { #for each replicate
      indR = (0:(c$nAlleles[m]-1))*c$nSamples[m] + r + c$startIndMarker_nAllelesReps[m] #obtain indices for replicate r
      heights = c$peakHeights[indR] #obtain PHs
      sumY[r,m] <- sum(c$peakHeights[indR])/markerEff #take sum of the peak heights (and scale with marker efficiency)
      meanbp[r,m] <- mean(c$basepairs[indA + 1])  #note adding of index here!
    }
  }
  useDEG = !all(meanbp==0)
  sumY <- colMeans(sumY) #take average
  meanbp <- colMeans(meanbp) #take average
  if(!useDEG) meanbp = NULL
  
  #CALL FUNCTION TO PREFIT MODEL PARAMEtERS:
  th0 <- fitSUMmodel(y=sumY,x=meanbp,niter=10,delta=2,restrictDeg=FALSE,model=c$model)
  
  #Post handling of start parameters:
  th0[2] = th0[2]*0.9 #shrinking start value of omega (ALWAYS)
  return(th0)
}