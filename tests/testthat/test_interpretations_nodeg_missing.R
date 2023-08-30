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
evidData = importMPSsample(paste0(path,"/MPS_evid2p_missing.csv") ) #LOAD EVIDENCE DATA
evidData2 = importMPSsample(paste0(path,"/MPS_evid2p_missing2.csv") ) #LOAD EVIDENCE DATA
refData = importMPSsample(paste0(path,"/MPS_evid2p_refsMissing.csv") ) #LOAD REFERENCE DATA

#evidData[[1]][["PENTA E"]]
#refData[[1]][["PENTA D"]] #missing marker
#Hypothesis set: minor as POI
#Hp: ref1 + 1unknown
#Hd: 2unknowns
kit = "ForenSeq" #specify kit to include degradation (with param beta<1)
NOC = 2 #Number of contributors
fst0 = 0.01

hypHp = list(NOC=NOC,cond=c(1,0),fst=fst0) #(also include fst here)
hypHd1 = list(NOC=NOC,cond=c(0,0),fst=fst0,knownRef=1) 
hypHd2 = list(NOC=NOC,cond=c(0,0),fst=fst0,knownRef=c(1,2)) 
seed0 = 1
verbose = FALSE

#Framework 1: Missing markers not in file
model = "GA"
mleHp = inferEvidence(evidData,popFreq,refData,hypHp , calibration, verbose=verbose,model=model,seed = seed0)
mleHd1 = inferEvidence(evidData,popFreq,refData,hypHd1 , calibration, verbose=verbose,model=model,seed = seed0)
mleHd2 = inferEvidence(evidData,popFreq,refData,hypHd2 , calibration, verbose=verbose,model=model,seed = seed0)
expect(getUpperLR(mleHd1),32.81)
expect_equal(getUpperLR(mleHd2),NULL)

#Framework 2: Missing markers in file
mleHpB = inferEvidence(evidData2,popFreq,refData,hypHp , calibration, verbose=verbose,model=model,seed = seed0)
expect(mleHp$fit$loglik,mleHpB$fit$loglik)