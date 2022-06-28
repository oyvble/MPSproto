#' @title calcLogLikC_prediction
#' @description Likelihood calculation of (evidence) prediction model (uses calibration info)
#' @param par parameters A list(mx,mu,omega,beta), The ordinary scale (ttheta)
#' @param c data object (returned from prepareC_calibration)
#' @param returnOutcome Whether to return all joint genotype outcome values (can be used for Deconvolution and marker-specific values)
#' @export 

calcLogLikC_prediction = function(par,c, returnOutcome=FALSE) {
  #keep_threshold A threshold used to 
   #stutter proportions should not exceed this
  
  #MARKER EFFICIENCY (Am):
  theta_Am = c$markerEfficiency #fixed and known
  #theta_Am = c(theta_Am, nLocs - sum(theta_Am))
  
  #Hypothesis:
  #as.numeric(bigsumVEC) vector of likelihood per combination, over all combinations ( P(E|g)*P(g) ) 
  #as.integer(combUseVEC), index of which genotypes to traverse 
  #c$nJointCombs #join genotype combinations per marker
  #c$NOC #number of contirbutors (constant)
  #c$NOK #Number of known contributors (per marker)
  
  #Other params: 
  #as.numeric(c$AT) #analytical/detection threshold (per locus)
  #as.numeric(c$fst) #theta-correction (per locus)
  #as.numeric(c$nNoiseParam) #number of noise (per locus)
  #as.numeric(c$noiseSizeWeight) #noise size weighting (per allele)

  #names(c)
  #data structure and observations
  #as.integer(c$nLocs) #number of loci
  #as.integer(c$nSamples) #number of samples per locus
  #as.integer(c$nAlleles) #number of observed alleles (inlcuding Q)
  #as.integer(c$nGenos) #number of genotype combination (1 contributor)
  #as.numeric(c$peakHeights) #signal intensity vector (including Q)
  #as.numeric(c$peakHeights2) #Transformation of the signals (used in pdf calculation)
  #as.numeric(c$freq) #frequency vector
  #as.numeric(c$nTyped) #Number of previously typed alleles. Used for theta-correction
  #as.numeric(c$maTyped) #Number of previously typed alleles (per allele). Used for theta-correction
  #as.numeric(c$basepairs) #adjusted fragmentlength in basepairs (used for degradation)
    
  #stutter part:  
  #as.integer(c$nStutters) #Number of stutters
  #as.integer(c$stuttFromInd) #index stutters from (starts from 0)
  #as.integer(c$stuttToInd) #index stutters to  (starts from 0)
  #as.numeric(c$stuttExp) #stutter expectations to use
  
  #Allele outcome/contribution  
  #as.integer(c$outG1allele) #Allele index for genotype outcome
  #as.integer(c$outG1contr) #Allele contriubtion for genotype outcome
  
  #Cumulative outcome   (start index for markers)
  #as.integer(c$startIndMarker_nAlleles)  #Startindex for number of alleles
  #as.integer(c$startIndMarker_nAllelesReps)  #Startindex for number of alleles (including replicates)
  #as.integer(c$startIndMarker_nStutters) #Startindex for number of stutters
  #as.integer(c$startIndMarker_nGenos) #Startindex for number of genotypes
  #as.integer(c$startIndMarker_outG1allele) #Startindex for number of alleles
  #as.integer(c$startIndMarker_outG1contr) #Startindex for number of alleles
  #as.integer(c$startIndMarker_nJointCombs) #Startindex for number of joint genotypes

  nJointCombs = c$nGenos^c$NOU #obtain number of combinations per markers
  startIndMarker_nJointCombs = cumsum(c(0,nJointCombs)) 
  nCombs = sum(nJointCombs) #total number of combinations
  
  #Select C script to run depending on chosen model:
  model = c$model #chosen model is indicated in c object
  cfun_allcomb = paste0("loglikPrediction_allcomb_",model)
  cfun_allcomb2 = paste0("loglikPrediction_allcomb2_",model)
  if(!returnOutcome) {
    calc = .C(cfun_allcomb,as.numeric(0), as.integer(c$nJointCombs), c$NOC, c$NOK,
              as.numeric(par$mx),  as.numeric(par$mu), as.numeric(par$omega), as.numeric(par$beta), as.numeric(theta_Am),
              c$AT,c$fst,c$nNoiseParam,c$noiseSizeWeight,
              c$nLocs, c$nSamples, c$nAlleles, c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps,
              c$peakHeights,c$peakHeights2, c$freq, c$nTyped, c$maTyped, c$basepairs,
              c$nGenos, c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele, c$startIndMarker_outG1contr, 
              c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttExp , c$startIndMarker_nStutters,
              c$knownGind) 
    ret = calc[[1]] #returning loglik
  } else {
    bigsumVEC = rep(0.0,nCombs)
    calc = .C(cfun_allcomb2,as.numeric(bigsumVEC), as.integer(c$nJointCombs), c$NOC, c$NOK,
              as.numeric(par$mx),  as.numeric(par$mu), as.numeric(par$omega), as.numeric(par$beta), as.numeric(theta_Am),
              c$AT,c$fst,c$nNoiseParam,c$noiseSizeWeight,
              c$nLocs, c$nSamples, c$nAlleles, c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps,
              c$peakHeights,c$peakHeights2, c$freq, c$nTyped, c$maTyped, c$basepairs,
              c$nGenos, c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele, c$startIndMarker_outG1contr, 
              c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttExp , c$startIndMarker_nStutters,
              c$knownGind, as.integer(startIndMarker_nJointCombs)) 
  
    #INSERT JOINT GENOTYPE NAMES FOR EACH CALCULATED COMBINATIONS
    bigsumVEC = calc[[1]]
    #Whether to obtain the likelihood of each specific genotypes
    likOutcome = list() #index with indexes to keep (updated after calculations)
    for(locind in seq_len(c$nLocs)) {       #Running through each marker
      #locind=1
      loc = c$locNames[locind] #obtain locus names
      nGjoint = c$nJointCombs[locind] #number of joint combinations
      
      indJointGenotypes = c$startIndMarker_nJointCombs[locind] + seq_len(nGjoint) #obtain index of unknown genotpyes
      probOutcome = bigsumVEC[indJointGenotypes] #take subset of values (all for specific marker)
      
      G = c$Gset[[loc]] #obtain genotype outcome
      G2 = paste0(G[,1],"/",G[,2]) #join into a vector
      nG = c$nGenos[locind] #nrow(G) #number of genotypes
      nU = c$NOU[locind] #numbers of unknowns for that marker
      
      #CONSTRUCT INFO ABOUT JOINT GENOTYPES FOR ALL JOINT COMBINATIONS:
      genOutcome = matrix(NA,nrow=nGjoint,ncol=c$NOC)
      knownGind = c$knownGind[ c$NOC*(locind-1) + 1:c$NOC] #obtain index of known genotypes
      posKnownInd = which(knownGind!=-1) #position of knowns
      posUnknownInd = which(knownGind==-1) #position of unknowns
      
      #INSERT GENOTYPE OF KNOWN CONTRIBUTORS:      
      if(length(posKnownInd)>0) { 
        knownGind = knownGind[posKnownInd]  #obtain indices
        for(pind in seq_along(posKnownInd)) genOutcome[,posKnownInd[pind]] = G2[knownGind[pind]+1]
      }
      
      #Index of joint genotype is ([0,0,..],[1,0,..],[2,0,..],...,[0,1,..],[0,2,..],....)
      if(nU>0) {
        GindJoint = list()
        for(uind in seq_len(nU)) GindJoint[[uind]] = seq_len(nG)
        GindJoint = expand.grid(GindJoint) #obtain full span of genotypes (correspond to algorithm in C++)
        for(uind in seq_len(nU)) genOutcome[,posUnknownInd[uind]] = G2[GindJoint[,uind]]     
      }
      
      #INSERT GENOTYPE OF UNKNOWN CONTRIBUTORS:   
      colnames(genOutcome) = paste0("C",seq_len(c$NOC))
      df = data.frame(genOutcome,likProb=probOutcome)
      likOutcome[[loc]] =  df
    }#end for each markers
    ret = likOutcome
  }

#calc[[1]]
  #return log likelihood funciton:
  return(ret)
  
} #end function