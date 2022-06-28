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
mleHp = inferEvidence(evidData,popFreq,refData,hypHp , calibration, verbose=FALSE,model=model,seed = seed0)
mleHd = inferEvidence(evidData,popFreq,refData,hypHd , calibration, verbose=FALSE,model=model,seed = seed0)

#Calculated LR:
par = getPar(mleHp,mleHd)
log10LR = getlog10LR(mleHp,mleHd) #log10 scale

#CHECK VALS
expect(log10LR,27.878)
expect(par[,1],c(0.074,0.926,960.177,0.132,1))
expect(par[,2],c(0.928,0.072,959.821,0.137,1))

#Get model fit:
validhp = validMLE(mleHp,"Hp",verbose = FALSE,createplot = FALSE)
validhd = validMLE(mleHd,"Hd",verbose = FALSE,createplot = FALSE)

expect(validhp$cumProb[1:5],c(0.742,0.698,0.644,0.412,0.507))
expect(validhd$cumProb[1:5],c(0.839,0.769,0.639,0.697,0.605))

#Provide Deconvolution results:
DCtableHp = deconvolve(mleHp)$table2
DCtableHd = deconvolve(mleHd)$table2
expect(as.numeric(DCtableHd[1:5,5]), c(0.712,1,1,0.915,1))

#check loglik results based on R implementation
checkLogLikLocs(mleHp)
checkLogLikLocs(mleHd)


##############
#NegBin MODEL#
##############
model = "NB"
mleHp = inferEvidence(evidData,popFreq,refData,hypHp , calibration, verbose=FALSE,model=model,seed = seed0)
mleHd = inferEvidence(evidData,popFreq,refData,hypHd , calibration, verbose=FALSE,model=model,seed = seed0)

#Calculated LR:
par = getPar(mleHp,mleHd)
log10LR = getlog10LR(mleHp,mleHd) #log10 scale

#CHECK VALS
expect(log10LR,24.280)
expect(par[,1],c(0.073,0.927,912.500,0.269,1))
expect(par[,2],c(0.929,0.071,914.247,0.236,1))

#Get model fit:
validhp = validMLE(mleHp,"Hp",verbose = FALSE,createplot = FALSE)
validhd = validMLE(mleHd,"Hd",verbose = FALSE,createplot = FALSE)

expect(validhp$cumProb[1:5],c(0.991,0.764,0.724,0.355,0.534))
expect(validhd$cumProb[1:5],c(0.961,0.853,0.737,0.6,0.608))

#Provide Deconvolution results:
DCtableHp = deconvolve(mleHp)$table2
DCtableHd = deconvolve(mleHd)$table2
expect(as.numeric(DCtableHd[1:5,5]), c(0.689,1,1,0.946,1))

#check loglik results based on R implementation
checkLogLikLocs(mleHp)
checkLogLikLocs(mleHd)



