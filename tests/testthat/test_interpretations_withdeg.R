#library(MPSproto);library(testthat)
expect = function(x,y,dig=3) { #helpfunction with specified tolerance
  expect_equal(round(as.numeric(x),dig),round(as.numeric(y),dig))
}

pkg = path.package("MPSproto") #get package install folder
path = paste0(pkg,"/examples") #path of data
calibratedModelFile = paste0(path,"/tutorialDataCalibrated.RDS") #calibrated model
calibration = readRDS(calibratedModelFile) #Or which is already stored in MPSproto installation

#DATA IMPORT (in addition to calibration)
freqFile = paste0(path,"/ForenSeq_NIST1036cauc_FWbrack.csv") #frequency file to use
popFreq =  importMPSfreqs(freqFile)[[1]]
#Note: Ensure that the alleles are all in forward direction for all markers (not UAS)

#Obtained generated examples:
evidData = importMPSsample(paste0(path,"/MPS_evid2p.csv") ) #LOAD EVIDENCE DATA
refData = importMPSsample(paste0(path,"/MPS_evid2p_refs.csv") ) #LOAD REFERENCE DATA

#OBTAIN CE FOR SEQUENCES IN EVIDENCE PROFILE (NECESSARY FOR DEGRADATION): 
#USING LUSstrR for this (already in forward direction)
if(0) {
  require(LUSstrR,warn.conflicts = FALSE)
  df = NULL
  for(loc in names(evidData[[1]]) ) {
#    loc="D1S1656"
    from = evidData[[1]][[loc]]$adata
    seq = unlist(LUSstrR::convertBracket2seq(from)) #obtain sequence
    rows = data.frame(Locus=loc,Sequence=seq,From=from)    
    df = rbind(df,rows)
  }
  df2 = LUSstrR::convert(df,hasFlanks = FALSE,format = "STRaitRazor") #dont reverse complement (already forwarded)
  df$CE = df2$Traditional_STR_Allele #obtain CE
  #View(subset(df,select=-Sequence))
  
  evidDataCE = evidData
  for(loc in names(evidData[[1]]) ) {
    df2 = subset(df,Locus==loc)
    from = evidDataCE[[1]][[loc]]$adata #get 
    if(!all(from==df2$From)) stop("wrong")
    to = paste0(df2$CE,":",from) #attach CE in front
    evidDataCE[[1]][[loc]]$adata = to #replace
  }
  plotMPS(evidDataCE,grpsymbol = ":")
  
  write.csv2(sample_listToTable(evidDataCE),file="MPS_evid2p_CE.csv",row.names = FALSE)
}

evidDataCE = importMPSsample(paste0(path,"/MPS_evid2p_CE.csv") ) #LOAD EVIDENCE DATA (with CE)

#Hypothesis set: minor as POI
#Hp: ref1 + 1unknown
#Hd: 2unknowns
kit = "ForenSeq" #specify kit to include degradation (with param beta<1)
NOC = 2 #Number of contributors
fst0 = 0.01

hypHp = list(NOC=NOC,cond=c(1,0),fst=fst0) #(also include fst here)
hypHd = list(NOC=NOC,cond=c(0,0),fst=fst0,knownRef=1) #(also include fst here)

seed0 = 1

#HELPFUNCTION TO COMPARE WITH INTERNAL R FUNCTION ("MANUAL")
checkLogLikLocs = function(mle) {
  # par=mle$fit$par;c=mle$prepareC
  likList = calcLogLikR_prediction(mle$fit$par,mle$prepareC)  #Obtain likelihood calculations
  #sapply(likList$lik_outcomeList,sum)==likList$logLik_marker
  logLikLocs = logLiki(mle)
  expect(likList$logLik_marker,logLikLocs)
  #list(R=likList$logLik_marker,C=logLikLocs)
}

#############
#GAMMA MODEL#
#############
model = "GA"
mleHp = inferEvidence(evidDataCE,popFreq,refData,hypHp , calibration, verbose=FALSE,model=model,seed = seed0,kit=kit)
mleHd = inferEvidence(evidDataCE,popFreq,refData,hypHd , calibration, verbose=FALSE,model=model,seed = seed0,kit=kit)

#Calculated LR:
par = getPar(mleHp,mleHd)
log10LR = getlog10LR(mleHp,mleHd) #log10 scale

#CHECK VALS
expect(log10LR,28.342)
expect(par[,1],c(0.072,0.928,1030.756,0.097,0.799))
expect(par[,2],c(0.930,0.070,1030.450,0.098,0.799))

#Get model fit:
validhp = validMLE(mleHp,"Hp",verbose = FALSE,createplot = FALSE)
validhd = validMLE(mleHd,"Hd",verbose = FALSE,createplot = FALSE)

expect(validhp$cumProb[1:5],c(0.836,0.79,0.76,0.252,0.371))
expect(validhd$cumProb[1:5],c(0.881,0.861,0.752,0.482,0.446))

#Provide Deconvolution results:
DCtableHp = deconvolve(mleHp)$table2
DCtableHd = deconvolve(mleHd)$table2
expect(as.numeric(DCtableHd[1:5,5]), c(0.655,1,1,0.925,1))

#check loglik results based on R implementation
checkLogLikLocs(mleHp)
checkLogLikLocs(mleHd)

##############
#NegBin MODEL#
##############
model = "NB"
mleHp = inferEvidence(evidDataCE,popFreq,refData,hypHp , calibration, verbose=FALSE,model=model,seed = seed0,kit=kit)
mleHd = inferEvidence(evidDataCE,popFreq,refData,hypHd , calibration, verbose=FALSE,model=model,seed = seed0,kit=kit)

#Calculated LR:
par = getPar(mleHp,mleHd)
log10LR = getlog10LR(mleHp,mleHd) #log10 scale

#CHECK VALS
expect(log10LR,22.208)
expect(par[,1],c(0.072,0.928,1005.112,0.221,0.811))
expect(par[,2],c(0.932,0.068,1020.330,0.172,0.787))

#Get model fit:
validhp = validMLE(mleHp,"Hp",verbose = FALSE,createplot = FALSE)
validhd = validMLE(mleHd,"Hd",verbose = FALSE,createplot = FALSE)

expect(validhp$cumProb[1:5],c(0.994,0.796,0.78,0.138,0.288))
expect(validhd$cumProb[1:5],c(0.971,0.917,0.818,0.409,0.361))

#Provide Deconvolution results:
DCtableHp = deconvolve(mleHp)$table2
DCtableHd = deconvolve(mleHd)$table2
expect(as.numeric(DCtableHd[1:5,5]), c(0.857,1,1,0.957,1))

#check loglik results based on R implementation
checkLogLikLocs(mleHp)
checkLogLikLocs(mleHd)
#paste0(round(validhp$cumProb[1:5],3), collapse=",")

#paste0(round(as.numeric(DCtableHd[1:5,5]),3), collapse=",")

