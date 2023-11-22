#' @title prepareC_prediction
#' @author Oyvind Bleka
#' @description Helpfunction to prepare data for calling C (used in calcLogLikC_prediction)
#' @details Assumes that all PH below threshold are removed. The function builds the data input to the C-code
#' @param dat A list returned from prepareData_prediction function (includes evididence,reference,frequency data)
#' @param hypothesis A list which defines the hypothesis to evaluate. Must contain the following elements:
#' NOC Number of contributors in model (Mandatory).
#' condOrder (conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' fst The co-ancestry coefficient. Can be a vector (must contain the marker names). Default is 0.
#' knownRef Specify known non-contributing references from refData (indices). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param calibration A list from calibration results (per marker). Must contain the following elements (in following order):
#' AT Analytical threshold 
#' markerEff Marker efficiency param
#' nNoiseParam Parameter for modeling distr of number of noise (Noise model)
#' noiseSizeParam  Parameter for modeling size of noise alleles (Noise model)
#' Reg-coeffs (b0,b1,b2) for different stutter types (each element considers a particular stutter type)
#' @param kit shortname of kit. Used to model degradation model. Obtained from getMPSkit(). 
#' @param platform Platform type {"MPS","CE"}. NOTE: CE has only consider only some stutter types (BW1,DBW1,FW1)
#' @param model Selected model for the signals (read/peak heights): {"GA"=gamma,"NB"=negative binomial}
#' @return ret A list of data input to call the C-code (used in calcLogLikC_prediction)
#' @export 

#kit=NULL;platform="MPS"
prepareC_prediction = function(dat, hypothesis, calibration, kit=NULL, platform="MPS",model="GA"){
 Qallele="99" #Name of allele given if missing in evidence (drop-out allele). Defualt is 99. This is important when considering the degradation model since 99 is closest to maximum allelein a locus. 
 #LUSsymbol="_" #a character symbol used to separate repeatunit and LUS.
 MPSsymbol = ":" #Used for extracting CE for MPS strings in evidence data. Example is "10:[ATCG]10".
 splitsymb = "/" #split symbol for serveral SI and BLMM values (returned from getStutterIndex)
 noiseCap = 1000 #assigning a noise cap (notice same default in inferNoiseModel function)
 
 #CHECK THAT SAME MARKERS ARE INCLUDED
 names(dat) = toupper(names(dat))
 names(calibration) = toupper(names(calibration))
 locs_dat = names(dat)
 locs_cal = names(calibration)
 if(!all(sort(locs_dat)==sort(locs_cal))) {
   locs_missing = setdiff(locs_cal,locs_dat)
   txt = paste0("Missing markers not found in evidence: ",paste0(locs_missing,collapse="/"),". Using markers of calibration instead of evidence.")
   print(txt)
   locs = locs_cal
 } else {
   locs = locs_dat #This is the order that we will use (same as data)
 } 
 nLocs = length(locs)
 
 sampleNames=names(dat[[1]]$samples) #obtaining sample names
 nReps = length(sampleNames) #number of replicates
 
 modelDegradation = FALSE
 if(!is.null(kit)) {
   kitinfo = getMPSkit(kit)  #Obtain kitinfo for selected kit (only short name of kit name given)
   if(is.na(kitinfo)[1]) stop("Specified kit not found. Needs to be specified in order to model degradation.")
   modelDegradation = TRUE
 }

 #obtain variables regarding hypothesis:
 NOC = hypothesis$NOC #obtain number of contributors
 condOrder = hypothesis$cond #obtain contribution order
 fst0 = hypothesis$fst
 if(is.null(fst0)) fst0 = 0 #insert as zero if missing
 knownRef = hypothesis$knownRef

 #Get list of genotypes for each markers (MUST BE SAME ORDER AS THOSE CREATED IN C++ code!!)
 Gset <- list() 
 for(loc in locs) {
    freq = dat[[loc]]$freq
    if(length(freq)==0) stop(paste0("Marker ",loc," not found frequency data!"))
    alleles <- names(freq)   
    nn = length(alleles)
    
    #G matrix is the vectorized upper triangular (1,1),(1,2),...,(1,n),(2,2),(2,3),...,(2,n),....,(n,n)
    G = numeric()
    for(i in 1:nn) {
       G = rbind(G, cbind( alleles[rep(i,nn - i + 1)], alleles[i:nn] ))
    }
    Gset[[loc]] <- G #store genotype
 }
 nGenos = sapply(Gset,nrow) #number of genotypes
 
 #Fix references as known contributors: Assign genotypes of known references to knownGind-matrix
 NOK = rep(0,nLocs) #number of known contributors per loci
 knownGind <- matrix(-1,ncol=nLocs,nrow=NOC) #default is no references (=-1)
 #assign references to knownGind-matrix by values of Glist
 if(!is.null(condOrder) && any(condOrder>0)) {
   for(loc in locs) {
     locind = which(loc==locs)
     subRef <- dat[[loc]]$refs #take out relevant reference
     if(length(subRef)==0)  stop(paste('Missing locus (',loc,') in refData.',sep=''))
     for(k in seq_along(subRef)) { #for each reference
       if(length(subRef[[k]])==0) next #Allowing unknown contributors for missing markers.
       if(length(subRef[[k]])!=2) stop("References need to have exactly two alleles.") #Must have 2 alleles
       if(condOrder[k]>0) { #If reference should be a contributor
         NOK[locind] = NOK[locind] + 1 #add known contributor
         Gind1 <- subRef[[k]][1]==Gset[[loc]][,1] & subRef[[k]][2]==Gset[[loc]][,2]
         Gind2 <- subRef[[k]][2]==Gset[[loc]][,1] & subRef[[k]][1]==Gset[[loc]][,2]
         knownGind[condOrder[k],locind] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
       }
     }
   }
 }

 #TRAVERSE EACH MARKER:
 #PREPARE FREQUENCY AND COVERAGE INFO
 datList = list() 
 stuttParamIndCounter = 0 #counting index for stutters
 for(m in seq_along(locs)) { 
# m=27
   datList[[m]] = list()
   
   loc = locs[m] #for selected loc (already upper)
   freq = dat[[loc]]$freq #get allele freqs
   allelesALL = names(freq) #allele also including Q-alele
   alleles = setdiff(allelesALL,Qallele) #assumed order as in freqs (assume alleles already added!)
      
   datList[[m]]$freq = freq #STORE FREQS (Q-allele already included)
   datList[[m]]$alleles = allelesALL #alleles #store allele names (Include possible Q-alleles)
   datList[[m]]$fst <- fst0
   
   #PREPARE OBSERVATIONS AND DROPIN WEIGHTS
   datList[[m]]$nAlleles <- nAlleles <- length(allelesALL) #number of allees to travers through (Including allele 99)
   datList[[m]]$nSamples <- length(dat[[loc]]$samples) #obtain number of samples
   
   #LAST: Obtain contribution matrices (per locus)
   numGenos1p = nGenos[loc] #obtain number of genotypes (1 contr)
   outG1contr = matrix(0,nrow=numGenos1p, ncol = nAlleles, 0.0); #init nG1xnA matrix (contribution matrix). Indicating what alleles that are contributoed
   outG1allele = matrix(0,nrow=numGenos1p,ncol=2) #init nG1x2 matrix (allele names as indices 0,...,nA-1)
   cc = 1 #0; #counter oveer all
   for (i  in 0:(nAlleles-1)) {
     for (j in i:(nAlleles-1)) { 
       outG1allele[cc,1] = i; #include index
       outG1allele[cc,2] = j; #include index
       outG1contr[cc,i+1] = outG1contr[cc,i+1] + 1.0; #insert contr at allele i
       outG1contr[cc,j+1] = outG1contr[cc,j+1] + 1.0; #insert contr at allele j
       cc = cc + 1; #iterate to next genotype outcome
     }
   }
   datList[[m]]$outG1allele = outG1allele
   datList[[m]]$outG1contr = outG1contr
   
   #Get number of typed allelse (must count alleles)
   tmp <- rep(0, length(allelesALL))
   if(!is.null(dat[[loc]]$refs)) {
     typedRefs = unique( c(which(condOrder>0),knownRef) ) #get unique referneces
     for(k in typedRefs) { #for each typed refs
       ind <- which( allelesALL%in%dat[[loc]]$refs[[k]] )
       tmp[ind] = tmp[ind] + (length(ind)==1) + 1 #add twice sampled if homozygote 
     }
   }
   datList[[m]]$nTyped <- sum(tmp) #number of total sampled (for each loci)
   datList[[m]]$maTyped <- tmp  #add vector of typed
   
   ###############################  
   #Obtain models from calibration
   ###############################  
   datList[[m]]$AT <- AT0 <- calibration[[loc]]$AT
   
   #Insert marker efficiency
   datList[[m]]$markerEff = calibration[[loc]]$markerEff[1] 

   #Obtain Noise model param:
   datList[[m]]$nNoiseParam <- calibration[[loc]]$nNoiseParam #obtain noise param (of geom.)
   datList[[m]]$noiseSizeParam <- noiseSizeParam0 <- calibration[[loc]]$noiseSizeParam[1] #obtain noise param (of pareto)

   #Obtain maximum obtainable noise:
   #upper = min(ceiling(AT0*prob0^(-1/noiseSizeParam0)),noiseCap) #formula for corresponding quantile
   discrete_range = seq(AT0,noiseCap)
   
   pdf_noiseCov = dpareto(discrete_range, noiseSizeParam0,AT0)
   pdf_noiseCov = pdf_noiseCov/sum(pdf_noiseCov) #normalize
   pdf_noiseCov_logged  = log(pdf_noiseCov) #obtain on logged scale
   cdf_noiseCov_logged = log(cumsum(pdf_noiseCov)) #also obtain cumulative (logged)
   #barplot(setNames(pdf_noiseCov,discrete_range))
   
   #Inesrt peak height vector   
   yv <- fv <- zv <- matrix(0,ncol=nAlleles,nrow=datList[[m]]$nSamples) #create Read-matrix and noize weight matrix (pdf and cdf)
   samples = dat[[loc]]$samples #obtain samples
   
   for(r in seq_along(samples)) { #for each replicates (following handle unordered loci under each sample)
     alleles_sample = names(samples[[r]]) #obtain alleles in sample 
     
     #INCLUDE THE POSSIBILITY THAT ALLELES CAN HAVE FORMAT "CE:STRING" (KEEP ONLY STRING)
     indChange = grep(MPSsymbol,alleles_sample)
     for(ind in indChange) alleles_sample[ind] = strsplit(alleles_sample[ind],MPSsymbol)[[1]][2] 
     
     yv[r, match(alleles_sample,alleles)] = samples[[r]] #insert COVERAGE

     #Modeling noise size as Discrete pareto: 
     #NOTE WE COULD SCALE THE WEIGHT HERE to make it discrete
     weight1 = rep( -100, length(yv[r,]) ) #init weights with very small weight (logged pdf)
     weight2 = rep( 0, length(yv[r,]) ) #this is logged cumulativ weights (default is 0)
     insVal = yv[r,]<=noiseCap & yv[r,]>=AT0
     if(any(insVal)) { #insert weight value if relevant to include
       vals = round(yv[r,insVal]) #be sure that it is whole integer
       matchind = match(vals,discrete_range)
       weight1[insVal] = pdf_noiseCov_logged[matchind] 
       weight2[insVal] = cdf_noiseCov_logged[matchind] 
     }
     fv[r,] = weight1 #insert
     zv[r,] = weight2 #insert
   }
   datList[[m]]$Yvec = c(yv) #insert PHs for obs alleles
   datList[[m]]$noiseSizeWeightPDF = c(fv) #insert noise weights, pdf (for each observations)
   datList[[m]]$noiseSizeWeightCDF = c(zv) #insert noise weights cdf (for each observations)
   
   #PREPARE DEGRADATION INFO
   if(modelDegradation) {
      subkitinfo = kitinfo[toupper(kitinfo$Marker)==toupper(loc),,drop=FALSE]
  	  if(nrow(subkitinfo)==0) stop("Marker not found in kitinfo!") 
  	  kitinfo_size =  subkitinfo$Size 
  	  kitinfo_allele = subkitinfo$Allele #assume name as allele names

  	  alleles0 = alleles #alleles to traverse (CE by default)
  	  if(platform=="MPS") {
  	    alleles_CE = unlist(lapply(samples,function(x) sapply(strsplit(names(x),MPSsymbol),function(x) x[1])))
  	    alleles_String = unlist(lapply(samples,function(x) sapply(strsplit(names(x),MPSsymbol),function(x) x[2])))
  	    
  	    indmatch = match(alleles,alleles_String) 
  	    alleles0 = alleles_CE[indmatch] #obtain CE for each corresponding string
  	    if(length(alleles0)!=length(alleles)) stop("Something went wrong when extracting CE from allele names (evidence)")
  	  }
  	  
      bp <- numeric() #fragment length/size in bp
      for(allel in alleles0) { #for each alleles (not Q-allele
        bind <- which( kitinfo_allele==allel ) #obtain index of position
        if(length(bind)==0) { #if missing alleles in kitinfo
            bind <- which.min(abs(as.numeric(allel) - as.numeric(kitinfo_allele))) 
        } 
        bp <- c(bp,kitinfo_size[bind]) #get fragment length
      } #end for each alleles
	  if(Qallele%in%allelesALL) bp = c(bp,max(kitinfo_size)) #set max size for Q-allele if considered
	  bp = (bp-125)/100 #rescale before
  } else {
	  bp = rep(0,nAlleles)
	}	   
  datList[[m]]$basepairs = bp #insert info
        
  #PREPARE STUTTER INFO	  
  stuttModel = NULL
  ntmp = length(calibration[[loc]]) #number of calibration objects
  if(ntmp>4) stuttModel = calibration[[loc]][5:ntmp]  #obtain stutter model
  
  if(!is.null(stuttModel)) {
    stuttTypes = names(stuttModel) #obtain stutter types
    stuttTypeLong <- stuttFromIndLong <- stuttToIndLong <- stuttExpLong <- NULL #store stuttered info
    
    for(type in stuttTypes) { #for each stutter type
#type = stuttTypes[1]
  	 stuttIndex  <- getStutterIndex(alleles,type,platform,locus=loc) #get stutter index for alleles
	 if( is.null(stuttIndex) ) next #skip if none of the alleles can provide stutter type
  	 stuttFromInd = stuttIndex$SI
  	 
  	 #Handle situations where from allele can have several "from canidates":
  	 #USING ALLELE WITH MOST READS
  	 hasMoreInd = grep(splitsymb,stuttFromInd) 
  	 for(ind in hasMoreInd) {
  	   stuttFromIndCands = as.integer(strsplit(stuttFromInd[ind], splitsymb)[[1]]) #obtain allele indices of candidates
  	   sumReps = colSums(yv[,stuttFromIndCands,drop=FALSE]) #obtain sum of each replicates
  	   indIns = which.max(sumReps) #select which to insert (the one with largest coverage)
  	   stuttFromInd[ind] = stuttFromIndCands[indIns] #use largest (most realistic)
  	   
  	   #Also update the BLMM vector (only include selected canidate)
  	   BLMMsplitted = strsplit(stuttIndex$BLMM[ind],splitsymb)[[1]]
  	   stuttIndex$BLMM[ind] = BLMMsplitted[indIns]
  	 } 
  	 stuttFromInd = as.integer(stuttFromInd) #avoid string on this
  	 
  	 isOK = stuttFromInd>0
  	 nStutt = sum(isOK) #number of stutters
  	 if(nStutt==0) next #skip if none observed
  	 stuttToInd = which(isOK) - 1 #subtract with 1 (c++)
  	 stuttFromInd = stuttFromInd[isOK] - 1 #subtract with 1 (c++)
  	 
  	 bhat = stuttModel[[type]] #obtain coefficients
  	 BLMM = stuttIndex$BLMM[isOK] #obtain BL of parental
  	 
  	 linPred = rep(bhat[1],nStutt) #obtain linear predictor for each stutters
  	 if(length(bhat)==2 && !is.na(bhat[2])) {
  	   linPred <- linPred + bhat[2]*as.integer(BLMM)  #include linear effect
  	 } else if(length(bhat)==3 && !is.na(bhat[3])) {
  	   BLMM2 = strsplit(BLMM,"-") #BLMM is multivariate
  	   bhat2 = bhat[-1] #consider var1+var2
  	   for(b in seq_len(nStutt) ) {
  	     linPred[b] <- linPred[b] + sum(bhat2*as.integer(BLMM2[[b]]))  #Need to include both types
  	   }
  	 }
  	 
  	 expStuttProp = 1/(1+exp(-linPred)) #obtain expected stutter proportion for each 
  	 
  	 stuttFromIndLong <- append(stuttFromIndLong, stuttFromInd)
  	 stuttToIndLong <- append(stuttToIndLong, stuttToInd)
  	 stuttExpLong <- append(stuttExpLong, expStuttProp)
  	 stuttTypeLong <- append(stuttTypeLong, rep(type,nStutt))
   } #loop through each stutter type

   datList[[m]]$nStutters = length(stuttFromIndLong) #total number of stutters to traverse
   datList[[m]]$stuttFromInd = as.integer(stuttFromIndLong)
   datList[[m]]$stuttToInd = as.integer(stuttToIndLong)
   datList[[m]]$stuttExp = as.numeric(stuttExpLong )
   datList[[m]]$stuttType = stuttTypeLong
   #datList[[m]]$alleles
  } else { #end if stutter model is included
    datList[[m]]$nStutters <- 0
    datList[[m]]$stuttFromInd <- datList[[m]]$stuttToInd <- as.integer()
    datList[[m]]$stuttExp = as.numeric()
    datList[[m]]$stuttType = as.character()
  }
 } #end for each marker
      
 
 #Vectorize list
 #Variables to use in the C code:
#  names(datList[[1]])
 nLocs #number of loci
 nSamples = as.integer(sapply(datList,function(x) x$nSamples)) #number of samples per locus
 nAlleles =  as.integer(sapply(datList,function(x) x$nAlleles)) #number of observed alleles (including Q) per locus
 #nAlleles2 =  as.integer(sapply(datList,function(x) x$nAlleles2)) #number of observed alleles (including Q) per locus
 freq =   as.numeric(unlist(sapply(datList,function(x) x$freq)))

 peakHeights = as.numeric( unlist(sapply(datList,function(x) x$Yvec)) ) #NOte the vectorization makes each replicates after each other for each allele
 nStutters = as.integer(sapply(datList,function(x) x$nStutters)) #number of stutters per locus
 stuttFromInd = as.integer( unlist(sapply(datList,function(x) x$stuttFromInd)) ) #- 1 #Note the -1 subtraction
 stuttToInd = as.integer( unlist(sapply(datList,function(x) x$stuttToInd)) ) #- 1 #Note the -1 subtraction
 stuttExp = as.numeric( unlist(sapply(datList,function(x) x$stuttExp)) ) #Obtain expectations
 stuttType = unlist(sapply(datList,function(x) x$stuttType))  #Obtain types
 
 markerEfficiency =  as.numeric( unlist(sapply(datList,function(x) x$markerEff)) ) #Obtain marker efficiencies expectations

 #Genotype traversing:
 nGenos =  as.integer(nGenos) #number of genotypes to traverse
 outG1allele =  as.integer( unlist(sapply(datList,function(x) t(x$outG1allele))) ) #obtain genotype outcome
 outG1contr = as.integer( unlist(sapply(datList,function(x) t(x$outG1contr))) ) #obtain genotype outcome
 
 #Diverse info:
 maTyped = as.numeric( unlist(sapply(datList,function(x) x$maTyped)) )
 nTyped = as.numeric( sapply(datList,function(x) x$nTyped) )
 AT = as.numeric(sapply(datList,function(x) x$AT) )
 
 fst = as.numeric(sapply(datList,function(x) x$fst) )
 nNoiseParam = as.numeric(sapply(datList,function(x) x$nNoiseParam) )
 noiseSizeParam  = as.numeric(sapply(datList,function(x) x$noiseSizeParam) ) #also keep this
 noiseSizeWeightPDF = as.numeric( unlist(sapply(datList,function(x) x$noiseSizeWeightPDF)) ) #NOte the vectorization makes each replicates after each other for each allele
 noiseSizeWeightCDF = as.numeric( unlist(sapply(datList,function(x) x$noiseSizeWeightCDF)) ) #NOte the vectorization makes each replicates after each other for each allele
 basepairs =  as.numeric( unlist(sapply(datList,function(x) x$basepairs)) ) 
 
 #CREATE CUMULATIVE vectors to make simple lookup in C++ code (per marker)
 startIndMarker_nAlleles = as.integer(c(0,cumsum(nAlleles))) #used to loop peak heights
 #startIndMarker_nAlleles2 = c(0,cumsum(nAlleles2)) #used to init shape vector
 startIndMarker_nStutters = as.integer(c(0,cumsum(nStutters))) #used for looping stutter indices
 startIndMarker_nGenos = as.integer(c(0,cumsum(nGenos)))
 startIndMarker_nAllelesReps = as.integer(c(0,cumsum(nAlleles*nSamples))) #used to loop peak heights (Replicates)
 
 #Vectorize outG1 matrices and calculate start index (per locus) for these
 startIndMarker_outG1allele = as.integer(2*startIndMarker_nGenos) #start index for 1p contributor matrix
 startIndMarker_outG1contr = as.integer(c(0,cumsum(nAlleles*nGenos))) #start index for per-allele contribution matrix

 #REMAINING variables:
 NOU = as.integer(NOC - NOK) #number of unknowns
 nJointCombs = as.integer(nGenos^NOU) #joint combinations
 startIndMarker_nJointCombs = as.integer(c(0,cumsum(nJointCombs)))
 if(any(is.na(nJointCombs))) stop("The number of combinations to traverse exceeded the maximum. Please reduce the number of alleles or unknowns contributors in the hypothesis.")
 
 #relatedness etc
 NOC=as.integer(NOC)
 NOK=as.integer(NOK)
 knownGind=as.integer(knownGind)

 #t(matrix(outG1contr[1:(nAlleles[1]*nGenos[1])],nrow=nAlleles[1]))
 #datList[[1]]$outG1contr
 alleleNames = lapply(datList,function(x) x$alleles) #obtain allele names for each marker (only observed)
 names(alleleNames) = locs #insert locus names

 #Special handle peak heights depending on model (pre-transform)
 if(model=="NB") {
   peakHeights2 = lgamma(peakHeights+1) 
 } else if(model=="GA") {
   peakHeights2 = log(peakHeights)
   peakHeights2[is.infinite(peakHeights2)] = 0 #will not be used
 }
 
 retlist = list( NOC=NOC, NOK=NOK, NOU=NOU,knownGind=knownGind, #hypotheseis
              nLocs=nLocs,nSamples=nSamples,nAlleles=nAlleles,nGenos=nGenos, nJointCombs=nJointCombs,#dimensions
              freq=freq,peakHeights=peakHeights,peakHeights2=peakHeights2, basepairs=basepairs, markerEfficiency=markerEfficiency,#data
              AT = AT,fst=fst, nTyped=nTyped,maTyped=maTyped, #settings data
              nNoiseParam=nNoiseParam,noiseSizeWeight=noiseSizeWeightPDF,noiseSizeWeightCDF=noiseSizeWeightCDF, #noise model
              nStutters=nStutters,stuttFromInd=stuttFromInd,stuttToInd=stuttToInd,stuttExp=stuttExp,stuttType=stuttType, #stutter info
              startIndMarker_nAlleles=startIndMarker_nAlleles,startIndMarker_nAllelesReps=startIndMarker_nAllelesReps,
              startIndMarker_nStutters=startIndMarker_nStutters,startIndMarker_nGenos=startIndMarker_nGenos, startIndMarker_nJointCombs=startIndMarker_nJointCombs,#cumulative
              outG1allele=outG1allele, outG1contr=outG1contr,startIndMarker_outG1allele=startIndMarker_outG1allele, startIndMarker_outG1contr=startIndMarker_outG1contr,  #contributor matrix
              locNames=locs,sampleNames=sampleNames,nReps=nReps,Gset=Gset,alleleNames=alleleNames,model=model)  #meta info
 
 return(retlist)
} #end function

