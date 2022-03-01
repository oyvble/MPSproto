#' @title createDiagnosticPlots 
#' @description A function to plot the diagnostic of the MarkerEfficienc
#' @param obj An returned object from the inferMarkerEfficiency function
#' @export 

createDiagnosticPlots = function(obj) {

  locNames = obj$locNames
  sampleNames = obj$sampleNames
  MLEobj=obj$MLE
  MCMCobj=obj$mcmc
  
  nLocs = length(locNames) #number of loci
  nSamples = length(sampleNames) #number of samples
  
  MEAN=MLEobj$theta2 # (approx posterior mean)
  SD=sqrt(diag(MLEobj$Sigma_theta2)) # (approx posterior SD)
  
  #1) Params: Am Marker efficiceny (Am)
  indUse_Am = 1:nLocs
  post_Am = MCMCobj$posttheta[,indUse_Am,drop=F]
  
  pdf("diagnoseMCMC_Am.pdf",height=2*nLocs,width=10)
  diagnoseMCMC(post_Am,xlim=range(post_Am),MEAN=MEAN[indUse_Am],SD=SD[indUse_Am])
  dev.off()
  
  #2) Params: Coverage expectation per sample (Mu)
  indUse_mu = nLocs + 1:nSamples
  post_mu = MCMCobj$posttheta[,indUse_mu,drop=F]
  colnames(post_mu) = sampleNames
  
  pdf("diagnoseMCMC_mu.pdf",height=2*nSamples,width=10)
  diagnoseMCMC(post_mu,xlim=range(post_mu),MEAN=MEAN[indUse_mu],SD=SD[indUse_mu])
  dev.off()
  
  #3) Params: Coverage variation per sample (Omega)
  indUse_omega = nLocs + 1:nSamples + nSamples
  post_omega =  MCMCobj$posttheta[,indUse_omega,drop=F]
  colnames(post_omega) = sampleNames
  
  pdf("diagnoseMCMC_omega.pdf",height=2*nSamples,width=10)
  diagnoseMCMC(post_omega,xlim=range(post_omega),MEAN=MEAN[indUse_omega],SD=SD[indUse_omega])
  dev.off()
  
    
}

