#' @title genMPSevidence
#' @author Oyvind Bleka 
#' @description Function for generating replicated mixture samples given same contributors and model parameters.
#' @details genDataset samples random mixture peak heights given as gamma(rho*sum(h_k),tau), with h_k as peak height of k-te contributor.
#' genData conditions on alleles given by refData. Empty references are generated with population frequencies.

#' @param NOC Number of contributors
#' @param calibration An object indicating the fitted calibration model (marker efficiency,stutters,noise) (MLE based)
#' @param popFreq A list of allele frequencies for a given population.
#' @param mu Expected peak heights for a het. single contributor allele
#' @param omega Coeffecient of variance of peak heights.
#' @param beta Coeffecient related to degree of degradation.
#' @param sorted Boolean for wheter sorting the contributors with respect to decreasingly mixture proportions.
#' @param refData A list of given reference profiles given as refData[[i]][[s]]. Default is random from population. 
#' @param mx A vector of known mixture proportions. Default is random uniform.
#' @param nrep Number of replicates (same contributors) to generate. Default is 1.
#' @param kit Kitname for considered kit (shortname). Required for degradation model
#' @param platform The platform name that is used
#' @return List with elements theta,samples,refData where theta is the true parameters of the model. samples is a list with samples which for each samples has locus-list elements with list elements adata and hdata
#' @export
#' @examples
#' \dontrun{ 
#' kit = "ForenSeq"
#' pkg = path.package("MPSproto")
#' calib =  readRDS(paste0(pkg,"/paper_stutterChar/calibrated_MPSproto.RDS"))
#' popFreq =  importMPSfreqs(paste0(pkg,"/paper_stutterChar/freqFile_ForenSeqFWbrack_Norway.csv"))[[1]]
#' gen = genMPSevidence(calib,2,popFreq,mu=1000,omega=0.2,beta=1,kit=kit )
#' plotMPS(gen$samples,gen$refData,AT=10)
#' }

#calibration=calib
#NOC=2;mu=1000;omega=0.1;beta=1;sorted=FALSE;mx=NULL;refData=NULL;nrep=1; kit="ForenSeq";platform="MPS"
genMPSevidence = function(calibration, NOC, popFreq, mu=1000,omega=0.1, beta=1, sorted=FALSE,mx=NULL,refData=NULL,nrep=1, kit=NULL,platform="MPS") {
  if(mu<0 || omega<0) stop("Parameters must be positive!")
  noiseCap = 1000 #assigning a noise cap (notice same default in inferNoiseModel function)    

  #Indicate whether to model degradation or not:
  modelDegradation = FALSE
  if(beta!=1) {
    if(is.null(kit)) stop("kit must be specified in order to apply degradation model.") 
    modelDegradation = TRUE
    kitinfo = getMPSkit(kit)  #Obtain kitinfo for selected kit (only short name of kit name given)
  }
  
  #Preparet mixture proportions (mx)
  if(is.null(mx)) {
   mx <- rgamma(NOC,1) #generate a flat distribution of mx if not provided
   mx = mx/sum(mx) #rdirichlet(1,rep(1,nC))  #simulate mx for contributors
  }
  if(sorted)  mx  <- sort(mx,decreasing=TRUE)
  mx <- c(mx)
  
  #obtain marker info:  
  locNames = names(popFreq) #obtain locus names to consider
  nLocs = length(locNames) #number of loci
  
  #convert (mu,sigma) to gamma-parameters
  shape0 <- 1/(omega^2) #constant part for all contributors
  scale0 <- mu*omega^2 #constant part for all contributors

  sampleNames = paste0("Sample",1:nrep) #create sample names
  mixH = c(t(replicate(2,mx))) #obtain mixture proportions for each alleles
  
  nR = length(refData)
  refNames = names(refData)
  if(NOC>nR) refNames = c(refNames, paste0("genRef",1:(NOC-nR))) #include generate names
  
  ##init sample evid and refs
  samples <- refData2 <- list() 
  for(sample in sampleNames) samples[[sample]] = list()
  for(ref in refNames) refData2[[ref]] = list()
  
  
  for(loc in locNames) { #traverse through each markers
#  loc=locNames[1]
    calibLoc = calibration[[loc]]
    if(is.null(calibLoc)) warning(paste0("Could not find calibration for locus ",loc))
    
    AT0 = calibLoc$AT #analytical threshold
    if(length(AT0)==0) stop(paste0("Analytical threshold (AT) not found for ",loc))
    
    markerEff = calibLoc$markerEff #marker efficiency
    if(length(markerEff)==0) stop(paste0("Marker efficiency (markerEff) not found for ",loc))
    
    markerEff_exp = markerEff[1] #get expectation (first element)
    #A second element may also be presented indicating sample variation
    
    #Draw random samples for marker efficiency (assume different for different samples)
    #markerEffRND = abs(rnorm(1,markerEffNstat_exp,markerEffNstat_sd)); #hist(markerEffRND)
    markerEffRND = markerEff_exp
    
    #Obtain param of noise model:
    nNoiseParam = calibLoc$nNoiseParam
    noiseSizeParam = calibLoc$noiseSizeParam 
    
    #Stutter types (these are always last elements:
    stuttTypes = NULL #no stutter by defaul
    if(length(calibLoc)>4) stuttTypes = names(calibLoc[-c(1:4)])
    
    #BEGIN GENERATING PROFILE
    #Obtaining number of known references to condition on.  
    refDataLoc = lapply(refData,function(x) unlist(x[[loc]]))
    nR = length(refDataLoc)
    nU = NOC-nR #number of contributors (unknowns) to generate
    refNames2 = names(refDataLoc) #obtain relevant reference names
    isKnown = refNames%in%refNames2 #indicate which referenes are known
    
    #Simulate TRUE contributors (if any unknowns to generate)
    if(nU>0) simAllelesRefs <-  matrix(sample(names(popFreq[[loc]]),size=2*nU,prob=popFreq[[loc]],replace=TRUE),ncol=2)
    #dose = as.integer(simAllelesRefs[,1]==simAllelesRefs[,2]) + 1 #obtain dose per reference
    
    contrib = NULL #obtain contribution vector (across all contributors)
    counter = 1
    for(k in seq_len(NOC)) { #for each contributor
      if(isKnown[k]) {
        refAllele = unlist(refDataLoc[[refNames[k]]]) #obtain known alleles
        if( length(refAllele)!=2) stop(paste0("Wrong allele vector length for reference ",refNames[k]," at marker ",loc))
      } else { #not known (need to impute )
        refAllele <- simAllelesRefs[counter,] #obtain alleles from simulations
        counter = counter + 1 #add counter of unknowns
      }
      refData2[[refNames[k]]][[loc]] = refAllele
      contribTmp = setNames(rep(mx[k],2),refAllele)
      contrib = c(contrib,contribTmp)
    }    

    #Obtain unique alleles (true alleles)
    agg=aggregate(contrib,by=list(names(contrib)),sum) #aggregate contributions for each alleles
    trueAllele = agg[,1] #obtain unique alleles
    nTrueAllele = length(trueAllele)
    shapev0 = agg[,2]*markerEffRND*shape0 #must multiply with marker efficency
    names(shapev0) = trueAllele
    
    #Scale shape vector with degradation model (euroformix type)
    if(modelDegradation) {
   	  subkitinfo = kitinfo[toupper(kitinfo$Marker)==toupper(loc),,drop=FALSE]
  	  kitinfo_size =  subkitinfo$Size 
  	  kitinfo_allele = subkitinfo$Allele #assume name as allele names
      
      bp <- numeric() #fragment length/size in bp
      for(allel in trueAllele) { #for each (true) allele
        bind <- which( kitinfo_allele==allel ) #obtain index of position
        if(length(bind)==0) { #if missing alleles in kitinfo
          
          if(platform=="CE") { #alleles are numeric
            bind <- which.min(abs(as.numeric(allel) - as.numeric(kitinfo_allele))) 
          } else if(platform=="MPS") { #alleles are strings
            #dist = adist(allel,kitinfo_allele) #calculate number of string-missmatches
            #bind <- which.min(dist)[1]
			stop("Degradation not implemented yet for MPS")
          } 
        } 
        bp <- c(bp,kitinfo_size[bind]) #get fragment length
      } #end for each alleles
      degscale = beta^((bp-125)/100) #note: must be same as in prepareC-function
      shapev0 = degscale*shapev0 #scale shape param
    } #end if degradation
    
    #obtain stutters (STR or MPS):
    shapev = shapev0 #copy shape param
    alleles = trueAllele #used to append new alleles
    
	  for(type in stuttTypes) { #for each stutter type
	   # type = stuttTypes[1]
	    stuttInfo = getStutteredSequence(trueAllele,type,platform=platform)  #obtain stutters (with BLMM)
  		stuttAllele = stuttInfo$stutter
  		stuttBLMM = stuttInfo$BLMM
  		
  		isNA = is.na(stuttAllele) #get stutter alleles as NA (i.e. no stutter returned)
  		if(sum(!isNA)==0) next #skip if no obtained stutters
  		
  		stuttAllele = stuttAllele[!isNA]
  		stuttBLMM = stuttBLMM[!isNA]

  		#OBtain BLMM variables
  		stuttAlleleBLMM = strsplit(stuttBLMM ,split="/") #also extract correspondning BLMM
  		stuttAlleleBLMM1 = as.numeric(sapply(stuttAlleleBLMM,function(x) x[1]))
  		stuttAlleleBLMM2 = as.numeric(sapply(stuttAlleleBLMM,function(x) x[2]))
  		
  		#Constructing liear predictor to obtain expected stutter proportion
  		stuttCoef = calibLoc[[type]] #obtain coefficients for particular tpe
  		linPred = rep(stuttCoef[1],length(stuttAlleleBLMM)) #intercept (always used)
  		if(length(stuttCoef)<=2 && !is.na(stuttCoef[2])) linPred = linPred + stuttCoef[2]*stuttAlleleBLMM1 #Inlcude slope 1
  		if(length(stuttCoef)<=3 && !is.na(stuttCoef[3])) linPred = linPred + stuttCoef[3]*stuttAlleleBLMM2 #Inlcude slope 2
  		stuttProps = 1/(1+exp(-linPred)) #obtaining expected stutter proportions
  		
  		newAllele = stuttAllele[!stuttAllele%in%alleles] #obtain new alleles (because of stutter)
  		if(length(newAllele)>0) {
  		  alleles = c(alleles,newAllele) #add to vector
  		  shapev = c(shapev,rep(0,length(newAllele))) #extend shape vector
  		}
  		
  		stuttIndFrom = match(trueAllele[!isNA],alleles)
  		stuttIndTo = match(stuttAllele,alleles)
  		
  		shapev[stuttIndTo] = shapev[stuttIndTo] + stuttProps*shapev0[!isNA] #add
  		shapev[stuttIndFrom] = shapev[stuttIndFrom] - stuttProps*shapev0[!isNA]  #subtract
	  }
    names(shapev) = alleles
    
    #Draw stutters with rgamma:
    heights <- rgamma(nrep*length(alleles),shape=rep(shapev,nrep),scale=scale0) #shape/scale given. We round to integer later.
    heights = matrix(heights,ncol=nrep)
    colnames(heights) = sampleNames
    rownames(heights) = alleles
    
    #always simulate noise if params are given
    simNoise =!is.null(nNoiseParam) && !is.null(noiseSizeParam)
    if(simNoise) {
      discrete_range = AT0:noiseCap #outcome
      pdf_noiseCov = dpareto(discrete_range, noiseSizeParam,AT0)
      pdf_noiseCov = pdf_noiseCov/sum(pdf_noiseCov) #normalize
      #plot(discrete_range,pdf_noiseCov)
    }
      
    #insert samples for each repeat
    for(repind in 1:nrep) {
#        repind=1
      isObserved = heights[,repind]>=AT0
      allelesObs = alleles[isObserved] #obtain observations
      heightsObs = round(heights[isObserved,repind]) #obtain observations
      
      #SIMULATE Noise alleles in addition:  sequence errors of trueAllele (1 base errors etc)

      #NOTE: DROP-IN ALLELES SHOULD NOT BE EXPLAINED AS MODELED STUTTERS!
      nNoiseAlleles = rgeom(1,nNoiseParam)
      if(nNoiseAlleles>0) {
        noiseReads = sample(discrete_range, nNoiseAlleles, prob=pdf_noiseCov,replace=TRUE)

        #if(platform=="CE") {
        other_alleles = setdiff(names(popFreq[[loc]]),allelesObs) #get pool of other alleles
        noiseAlleles = sample(other_alleles,min(nNoiseAlleles, length(other_alleles)),replace=FALSE)
        #} else {
          
        #trueAlleleSeq = convertBracket2seq(trueAllele)
        #mutatedSeq = trueAlleleSeq[mutate]
        #noiseAllele = LUSstrR::convert(mutatedSeq) #convert back to bracket format
        
        allelesObs = c(allelesObs,noiseAlleles) #add noise allele
        heightsObs = c(heightsObs,noiseReads) #add noise coverage
      }
      
      #Insert to table:
      samples[[repind]][[loc]] = list(adata=allelesObs,hdata=as.numeric(heightsObs))
    } #end for each repeats
  } #end for each marker
  return( list(calibration=calibration, NOC=NOC,samples=samples,refData=refData2, mu=mu,omega=omega, mx=mx) )
}


