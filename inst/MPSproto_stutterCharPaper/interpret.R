#THIS SCRIPT SHOWS EXAMPLES OF MIXTURES INTERPRETATION BASED ON CALIBRATED MODEL
#We simulate and interpret fictive mixtures
library(MPSproto)
projname = "MPSproto_stutterCharPaper" 
path = path.package("MPSproto") #get package install folder
datfold = paste0(path,"/",projname,"/data")

#IMPORT CALIBRATION FILE
calibfile = paste0(datfold,"/calibrated_",projname) #file name of the calibrated object
fn_calibration = paste0(calibfile,".RDS")
calibration = readRDS(fn_calibration) #import calibration object

#IMPORT FREQENCY DATA
freqFile = paste0(datfold,"/ForenSeq_Norway_FWbrack.csv") #frequency file to use
popFreq =  importMPSfreqs(freqFile)[[1]]

#DATA IMPORTED  


fst0 = 0.01 #specified theta/fst (same for all hypotheses: Affects situations with unknown contibutors)

#EVALUATION OF A GENERATED 2p-MIXTURE
NOC = 2 #Number of contributors
set.seed(1) #make reproducible
gen = genMPSevidence(calibration,NOC,popFreq,mu=1000,omega=0.1)
gen$mx
evidData = gen$samples #this is evidence data
refData = gen$refData #this is reference data
plotMPS(evidData,refData) #Show in graphical interface

#############################
#SPECIFYING HYPOTHESIS SET 1# 
#############################
#Hp: 'ref1 + ref2' Hd: 'ref2 + 1 unknown'
hypHp = list(NOC=NOC,cond=c(1,2),fst=fst0) 
hypHd = list(NOC=NOC,cond=c(0,1),fst=fst0,knownRef=c(1)) 

#Fit model:
mleHp = inferEvidence(evidData,popFreq,refData,hypHp , calibration, verbose=TRUE)
mleHd = inferEvidence(evidData,popFreq,refData,hypHd , calibration, verbose=TRUE)

#Obtain calculated LR
par = getPar(mleHp,mleHd) #get estimated params
LR = getLR(mleHp,mleHd) #get calculated LR (unlogged)
log10LR = getlog10LR(mleHp,mleHd) #log10 scale
upperLR = getUpperLR(mleHd) #get upper LR-limit 1/RMP(theta)

#Perform model validation
validhp = validMLE(mleHp,"Hp")
validhd = validMLE(mleHd,"Hd")

#Show model fit:
plotTopMPS(mleHp)
plotTopMPS(mleHd)

#############################
#SPECIFYING HYPOTHESIS SET 2# 
#############################
#Hp: 'ref1 + 1 unknown' vs Hd: '2 unknowns'
hypHp = list(NOC=NOC,cond=c(1,0),fst=fst0) #(also include fst here)
hypHd = list(NOC=NOC,cond=c(0,0),fst=fst0,knownRef=c(1)) #(also include fst here)

#Fit model:
mleHp = inferEvidence(evidData,popFreq,refData,hypHp , calibration, verbose=TRUE)
mleHd = inferEvidence(evidData,popFreq,refData,hypHd , calibration, verbose=TRUE)

#Obtain calculated LR
par = getPar(mleHp,mleHd) #get estimated params
LR = getLR(mleHp,mleHd) #get calculated LR (unlogged)
log10LR = getlog10LR(mleHp,mleHd) #log10 scale
upperLR = getUpperLR(mleHd) #get upper LR-limit 1/RMP(theta)

#Perform model validation
validhp = validMLE(mleHp,"Hp")
validhd = validMLE(mleHd,"Hd")

#Show model fit:
plotTopMPS(mleHp)
plotTopMPS(mleHd)


#EVALUATION OF A GENERATED 3p-MIXTURE
NOC = 3 #Number of contributors
set.seed(1) #make reproducible
gen = genMPSevidence(calibration,NOC,popFreq,mu=1000,omega=0.1)
gen$mx
evidData = gen$samples #this is evidence data
refData = gen$refData #this is reference data
plotMPS(evidData,refData) #Show in graphical interface

#############################
#SPECIFYING HYPOTHESIS SET 1# 
#############################
#Hp: 'ref1 + ref2 + 1 unknown' Hd: 'ref2 + 2 unknown'
hypHp = list(NOC=NOC,cond=c(1,2,0),fst=fst0) 
hypHd = list(NOC=NOC,cond=c(0,1,0),fst=fst0,knownRef=c(1)) 

#Fit model:
mleHp = inferEvidence(evidData,popFreq,refData,hypHp , calibration, verbose=TRUE)
mleHd = inferEvidence(evidData,popFreq,refData,hypHd , calibration, verbose=TRUE)

#Obtain calculated LR
par = getPar(mleHp,mleHd) #get estimated params
LR = getLR(mleHp,mleHd) #get calculated LR (unlogged)
log10LR = getlog10LR(mleHp,mleHd) #log10 scale
upperLR = getUpperLR(mleHd) #get upper LR-limit 1/RMP(theta)

#Perform model validation
validhp = validMLE(mleHp,"Hp")
validhd = validMLE(mleHd,"Hd")

#Show model fit:
plotTopMPS(mleHp)
plotTopMPS(mleHd)

#############################
#SPECIFYING HYPOTHESIS SET 2# 
#############################
#Hp: 'ref1 + 2 unknowns' Hd: '3 unknowns'
hypHp = list(NOC=NOC,cond=c(1,0,0),fst=fst0) #(also include fst here)
hypHd = list(NOC=NOC,cond=c(0,0,0),fst=fst0,knownRef=c(1)) #(also include fst here)

#Fit model:
mleHp = inferEvidence(evidData,popFreq,refData,hypHp , calibration, verbose=TRUE)
mleHd = inferEvidence(evidData,popFreq,refData,hypHd , calibration, verbose=TRUE)

#Obtain calculated LR
par = getPar(mleHp,mleHd) #get estimated params
LR = getLR(mleHp,mleHd) #get calculated LR (unlogged)
log10LR = getlog10LR(mleHp,mleHd) #log10 scale
upperLR = getUpperLR(mleHd) #get upper LR-limit 1/RMP(theta)

#Perform model validation
validhp = validMLE(mleHp,"Hp")
validhd = validMLE(mleHd,"Hd")

#Show model fit:
plotTopMPS(mleHp)
plotTopMPS(mleHd)

