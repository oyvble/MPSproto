#' @title calibrateModel
#' @author Oyvind Bleka 
#' @description A wrapper function for performing calibration 
#' @details Helpfunction to obtain calibration objects ready for evidence predictions
#' 
#' @param dat Datatable with columns (RefName,SampleName,Locus,Allele,Coverage,Dose)
#' @param AT Defined analytical threshold used for the analysis (may be defined per marker)
#' @param platform Type of platform. Supports "MPS" or "CE".  
#' @param locNames The order of markers can be specified (otherwise extracted as unique from dataset)
#' @param minStuttOccurence Required number of occurence of a particular stutter type (avoiding overfit)
#' @param stutterTypes Stutter types. See getStutterTypeRule for supported types
#' @param createPlots Whether to create PDF files showing fitted noise model
#' @param runMCMCvalidation Boolean of whether to run MCMCvalidation for checking the posterior distributions
#' @param mleOptions Option input for optimization
#' @param mcmcOptions Option input for MCMC simulation
#' @param verbose Whether to print out progress
#' @return Calibration object used for interpreting a new evidence
#' @export


calibrateModel = function(dat,AT=10,platform="MPS",locNames=NULL,minStuttOccurence=3,stutterTypes=NULL,createPlots=FALSE,runMCMCvalidation=FALSE, mleOptions=list(maxIter=1000,delta=.1,seed=1,printLevel=0,penalty=1,omegaAdj=0.5),mcmcOptions=list(niter=100000,delta=0.01),verbose=TRUE) {
  
  #############################################################################
  #Step 0: Insert missing alleles with dose>0 (true alleles) (attached to end)#
  #############################################################################
  dat = insertMissingTrueAlleles(dat) 
  
  #Possible to specify order of markers with locNames argument (CASE SENSITIVE)
  if(is.null(locNames)) locNames = unique(dat[,3]) #obtain locus names from data if not specified (this will be the order)
  if(length(AT)==1) AT = setNames( rep(AT,length(locNames)),locNames)
    
  ############################################
  #Step 1: OBTAINING MARKER EFFICIENCY PARAMS#
  ############################################
  markerEfficiency = inferMarkerEfficiency(dat,AT,verbose = verbose,runMCMC = FALSE) #fitting model MLE model
  Amhat = markerEfficiency$MLE$theta2[locNames] #estimates
  #  markerEff = getMarkerEffStat(markerEfficiency)
  #plot(markerEff$MEAN,Amhat)
  
  ######################################
  #Step 2: OBTAINING STUTTER PARAMETERS#
  ######################################
  #SEARCH AFTER STUTTER DATA (including index of potential stutters of all types={BW1,DBW1,FW2})
  #dat2 = getFilteredData(dat,platform=platform) #obtain subset of data suitable for extracting stutters
  stutterData = getStutterData(dat,platform=platform,verbose=verbose,minStuttOccurence=minStuttOccurence,stutterTypes=stutterTypes)

  #Infer STUTTER MODELS BASED ON INFO FROM STEP 1 (one marker at the time):
  pdffile = "StutterModel" #name of pdf to be printed
  if(!createPlots) pdffile=NULL
  stuttMod = inferStutterModel(stutterData$stutterdata,print2pdf=pdffile)
  stuttModTab = stuttMod$betaTable #obtain coefficient table
  
  stuttModList = list()
  for(i in seq_len(nrow(stuttModTab))) {
    row = unlist(stuttModTab[i,])
    loc = row[1] #obtain locus name
    type = row[2] #obtain stutter type
    est =  as.numeric( na.omit(row[3:5]))    #obtain estimats
    if( !loc%in%names(stuttModList)) stuttModList[[loc]] = list()
    stuttModList[[loc]][[type]] = est
  }
  
  
	#########################################################
	#STEP 3: INFER NOISE MODEL AFTER ASSIGNING STUTTER MODEL#
	#########################################################
  pdffile = "NoiseModel" #name of pdf to be printed
  if(!createPlots) pdffile=NULL
	noiseMod = inferNoiseModel(sampleNames=unique(dat$SampleName), locNames=unique(dat$Locus), noisedata=stutterData$noisedata, AT=AT,print2pdf=pdffile)

  #FINALLY INSERT EVERYTHING TO ONE OBJECT
  calibration = list()
  for(loc in locNames) {
    markerEff= Amhat[loc] #setNames(c(,AmSD[loc],AmSD[loc]),c("est","SD","SD2"))
    noiseMods = lapply(noiseMod,function(x) x[loc])
    stuttMods = stuttModList[[loc]] 
  
    calibration[[loc]] = list(markerEff=markerEff)  
    calibration[[loc]] = append(calibration[[loc]],noiseMods)
  
    #put stutter model last (requirement!!)
    calibration[[loc]] = append(calibration[[loc]],stuttMods)
  }

  #RETURN PREPARED CALIBRATION OBJECT FOR PREDICTION (based on fitted calibration model);
  return(calibration)
}
