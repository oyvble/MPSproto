#' @title inferMarkerEfficiency 
#' @description Infer marker efficiency 
#' @details Using the sum of the coverage per marker as data (per sample), 
#' the function estimates marker efficiency and sample specific EXP(mu) and CV(omega). 
#' Using Maximum likelihood estimation for a gamma distribution. 
#' Requires specification of AT to account for dropout markers.
#' 
#' The function also performs MCMC for diagnostic: the mleObj argument allows for possibility for doing MCMC directly.
#' @param dat Datatable with columns (SampleName,Locus,Coverage)
#' @param AT Vector with analytical threshold per marker
#' @param locNames The order of markers can be specified (otherwise extracted as unique from dataset)
#' @param runMCMC Whether to run MCMC for inferring the posterior distribution
#' @param mleOptions Option input for optimization
#' @param mcmcOptions Option input for MCMC simulation
#' @param verbose Whether to print out progress
#' @param mleObj An output from inferMarkerEfficiency with mle fit (used for directly performing MCMC)
#' @export 

#locNames=NULL;runMCMC=FALSE;mleOptions=list(maxIter=1000,delta=.1,seed=1000); mcmcOptions=list(niter=10000,delta=.025,seed=1);verbose=TRUE
inferMarkerEfficiency = function(dat,AT,locNames=NULL, runMCMC=FALSE, mleOptions=list(maxIter=1000,delta=.1,steptol=1e-6,seed=1000), mcmcOptions=list(niter=10000,delta=.025,seed=1),verbose=FALSE,mleObj=NULL) {

  if( !all(c("SampleName","Locus","Coverage")%in%colnames(dat)) ) stop("Wrong or missing column names!")

  #Likelihood function to fit marker efficiency
  neglogLik = function(th,transform=FALSE,useNegative=FALSE) {
    if(transform) th = exp(th)
    Avec = th[1:(nLocs-1)]
    Am = nLocs - sum(Avec) #obtain last param
    Avec = c(Avec,Am) #insert last
    mu = th[ 1:nSamples + nLocs - 1]
    omega = th[ (nLocs + nSamples):length(th)]
    
    if(any(Avec<0)) {
      if(useNegative) return(Inf)
      if(!useNegative) return(-Inf)
    }
    shapemat = outer(2/omega^2,Avec)
    scalemat = replicate(nLocs,mu*omega^2)
    shapeVector = c(shapemat) #vectorize with number of sampleNames for each locis
    scaleVector = c(scalemat)  #vectorize with number of sampleNames for each locis
    
    #MODELING BOTH OBSERVED AND DROPOUT:
    val = sum(dgamma(sumCovObs[!markerDropout],shapeVector[!markerDropout],scale=scaleVector[!markerDropout],log=TRUE))
    if(any(markerDropout)) val = val + sum(pgamma( AT[ locusNamesVectorized[markerDropout] ] ,shapeVector[markerDropout],scale=scaleVector[markerDropout],log.p=TRUE))
    #curve(dgamma(x,shapeVector[markerDropout],scale=scaleVector[markerDropout]),from=0,to=1000);abline(v=AT)
    if(useNegative) {
      return(-val)
    } else {
      return(val)
    }
  }

  #Extracting unique sampleNames and loci
  sampleNames = unique(dat$SampleName) #obtain unique samples
  if(is.null(locNames)) locNames = unique(dat$Locus) #obtain unique loci
  nLocs = length(locNames) #number of loci
  nSamples = length(sampleNames) #obtain number of sampleNames
  if(verbose) print(paste0("Dataset consists of ",nSamples," samples with ",nLocs," markers."))
  
  #Obtain sum marker of each data
  sumCov = aggregate( dat$Coverage,by=list(dat$SampleName,dat$Locus),sum,drop=FALSE) #obtain per-marker
  colnames(sumCov) = c("Sample","Locus","value")
  
  if(length(AT)==1) AT = setNames(rep(AT,nLocs),locNames)
  names(AT) = toupper(names(AT)) #get as upper case
  
  #SET VALUES BELOW AT TO ZERO
  for(loc in locNames) sumCov$value[sumCov$Locus==loc & sumCov$value<AT[loc]] = 0
  
  ###########################
  #Step 1: obtain startvalue#
  ###########################
  
  #obtain locus specific parameters
  agg = aggregate( sumCov$value,by=list(sumCov$Locus),function(x) c(mean(x/2)),drop=FALSE)#,sd(x/2))) #obtain per-marker
  meanPHperLoc = setNames( agg$x,agg[,1] )[locNames] 
  mu0 = mean(meanPHperLoc) #mean(dat$Coverage)/2 #obtain expectead PH
  Avec0 = meanPHperLoc/mu0 #scale variable wrt mu0
  
  agg = aggregate(  sumCov$value,by=list(sumCov$Sample),function(x) c(mean(x/2),sd(x/2)),drop=FALSE)
  mu_samples = setNames(agg$x[,1],agg[,1])[sampleNames]
  sd_samples = setNames(agg$x[,2],agg[,1])[sampleNames]
  omega_samples = sd_samples/mu_samples*0.5 #obtain start value for omega variable
  
  #Obtain start value of params
  theta0 = c(Avec0[-length(Avec0)],mu_samples,omega_samples) #start value of param
  nparam = length(theta0) #number of param
  #  th=theta0
  
  #Obtain info about which markers that drops out based on "expected" data vector (filled)  
  expectedMarkers = c(t(replicate(nSamples,locNames))) #obtain expected markers
  expectedSamples = rep(sampleNames,nLocs) #obtain expected markers
  expectedMarkerSamples = paste0(expectedMarkers,expectedSamples) #obtain unique vector for expectation
  #  sumCov = sumCov[sample(1:nrow(sumCov)),]
  dataMarkerSamples = paste0(sumCov$Locus,sumCov$Sample) #obtain unique vector for data
  
#  all(expectedMarkerSamples==dataMarkerSamples): THIS WILL ONLY MATCH IF NO DROPOUT
  sumCovObs = sumCov$value[match(expectedMarkerSamples,dataMarkerSamples)] #change order to make it follow expectedMarkerSamples
  sumCovObs[is.na(sumCovObs)] = 0 #set as missing if not found
  #sum(sumCovObs==0)
  inferEFF = sum(sumCovObs>0)/nparam #inference Efficiency (number of datapoints per params)
  if(verbose) print(paste("Number of datapoints per parameter: ",round(inferEFF,1)))
  
  #  tail(sumCov,100)
  markerDropout = sumCovObs==0 #!expectedMarkerSamples%in%dataMarkerSamples #obtain index of dropout markers
  if(verbose) print( paste0("Number of marker dropouts=", sum(markerDropout)))
  locusNamesVectorized = toupper(expectedMarkers) #vectorize markers (needed for accessing locus specific AT )
  
    #  expectedMarkerSamples[which(markerDropout)]
    # tail( cbind(sumCov[,1:3],expectedMarkerSamples[!markerDropout]),10) #CHECK IF SAME ORDER
    #  logLik(th=log(theta0))
    
  if(is.null(mleObj)) {
    
    #FITTING MARKER EFFICIENCY MODEL
    if(verbose) print("Fitting Marker efficiency parameters using MLE...")
    set.seed(mleOptions$seed)
    #neglogLik(theta0,transform=FALSE)==neglogLik(log(theta0),transform=TRUE)
    #fit = nlm(neglogLik,log(theta0),iterlim = mleOptions$maxIter,transform=TRUE)
    fit = optimizeLik(neglogLik,log(theta0),maxIter=mleOptions$maxIter ,delta=mleOptions$delta, steptol=mleOptions$steptol, verbose=verbose)
    fit$phi = fit$estimate #obtain estimates (transformed scale)
    fit$theta = exp(fit$phi)  #transform back
    
    #obtain param names
    paramNames2 = c(locNames, paste0("mu.",sampleNames),paste0("omega.",sampleNames))
    paramNames = paramNames2[-nLocs]
    names(fit$theta) <- names(fit$estimate) <- paramNames
    
    #Also including last marker to Marker efficiency:
    indUse_Am = 1:(nLocs-1)
    indUse_Am2 = 1:nLocs #include full vector of Amhat
    Amhat = fit$theta[indUse_Am] #obtain Marker efficiency params
    A_M = setNames( nLocs - sum(Amhat), locNames[nLocs]) #obtain the restricted last param (deduced from the others)
    fit$theta2 = c(Amhat, A_M, fit$theta[nLocs:length(fit$theta)]) #include all params
  
  #  plot(fit$theta2[grep("omega",names(fit$theta2))],omega_samples);abline(a=0,b=1)
    
    #OBTAIN COVARIANCE MATRIX for both phi and theta (based on delta-method)
    Sigma_phi = solve(fit$hessian)#+diag(0.1,nparam)) #Calculating covariance matrix of proposal fujnction
    Jmat = diag(fit$theta) #this is Jacobian matrix
    Sigma_theta <- Sigma_theta2 <- t(Jmat)%*%Sigma_phi%*%Jmat #this is delta-method
    Sigma_theta2 <- cbind(Sigma_theta2[,indUse_Am],0,Sigma_theta2[,-indUse_Am]) #put in space for AmL (horisontal)
    Sigma_theta2 <- rbind(Sigma_theta2[indUse_Am,],0,Sigma_theta2[-indUse_Am,]) #put in space for AmL (vertical)
   #dim(Sigma_theta2)
    #isSymmetric(Sigma_theta2)
    COV_Am = Sigma_theta[indUse_Am,indUse_Am] #covariance matrix of Avec
    Sigma_theta2[nLocs,nLocs] = sum(COV_Am) #insert variance (sum of covariances)
    Sigma_theta2[indUse_Am,nLocs] <- Sigma_theta2[nLocs,indUse_Am] <-  -colSums(COV_Am) #this is covariance between Ai,AL
    #sum(Sigma_theta2==0) #remaining is zero
    #image(log(abs(Sigma_theta2)))
    #image(log(abs(Sigma_phi)))
    #Amhat_VAR = c(diag(thetaCOV_Am), sum(thetaCOV_Am)) #Obtain 
    colnames(Sigma_theta) <- rownames(Sigma_theta) <- paramNames
    colnames(Sigma_theta2) <- rownames(Sigma_theta2) <- paramNames2
    
    #Storing covariances:
    fit$Sigma_phi = Sigma_phi
    fit$Sigma_theta = Sigma_theta
    fit$Sigma_theta2 = Sigma_theta2
    
  } else { #end if mleObj not specified
    fit=mleObj$MLE
    paramNames2 = names(fit$theta2)
  }
  #RUNNNING MCMC TO OBTAIN POSTERIOR DISTRIBUTION
  mcmc = NULL
  if(runMCMC) {
    set.seed(mcmcOptions$seed)
    #Assuming flat prior on phi param (correspond to exponential in theta domain):
    if(verbose) print("Running MCMC...")
    theta0 = fit$theta #start values
    Sigma0 = fit$Sigma_theta
#    logLik=neglogLik;th0=fit$theta;Sigma0=Sigma_theta;niter=mcmcOptions$niter; delta=1
#    doMCMC(logLik=neglogLik, th0=fit$theta, Sigma0=Sigma_theta, niter=mcmcOptions$niter, delta=1)$accrat
    mcmc = doMCMC(logLik=neglogLik, th0=theta0, Sigma0=Sigma0, niter=mcmcOptions$niter, delta=mcmcOptions$delta)
#   print(mcmc$accrat)
    tmp = mcmc$posttheta #keep a copy
    mcmc$posttheta = cbind(tmp[,1:(nLocs-1)],0,tmp[,nLocs:ncol(tmp)])
    mcmc$posttheta[,nLocs] = nLocs - rowSums(tmp[,1:(nLocs-1),drop=FALSE]) #obtain last Marker efficiency variable
    colnames(mcmc$posttheta) = paramNames2
  }
  return(list(MLE=fit,mcmc=mcmc,locNames=locNames,sampleNames=sampleNames,AT=AT,sumCovObs=sumCovObs)) 
}