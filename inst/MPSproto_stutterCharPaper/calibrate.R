#CALIBRATION OF MPSproto based on paper data
#1) Marker efficiency
#2) Stutter model (already done)
#3) Noise model

#Note: This script differs from the tutorial since we import aggregated data instead of a full data table

library(MPSproto)
projname = "MPSproto_stutterCharPaper" 
path = path.package("MPSproto") #get package install folder
datfold = paste0(path,"/",projname,"/data")
calibfile = paste0(datfold,"/calibrated_",projname) #file name of the calibrated object
datfn_read = paste0(datfold,"/data_sumPerMarker.csv") #used for marker efficiency calibration
datfn_stuttermodel = paste0(datfold,"/calibrated_StutterCoef.csv") #the stutter model is already calibrated
datfn_noise = paste0(datfold,"/data_noise.csv") #used for noise model calibration (part 1)
datfn_stutterFew = paste0(datfold,"/data_excludedStutters.csv") #used for noise model calibration (part 2)

#SETTINGS:
AT=10 #ANALYTICAL THREHOLD (same for all markers)

#####################################
#1) CALIBRATION OF MARKER EFFICIENCY#
#####################################
#IMPORTANT: THESE DATA SHOULD NOT BE DEGRADED!
datSum = tableReader(datfn_read)
colnames(datSum) = c("SampleName","Locus","Coverage") #ensure correct header names

#Following is needed for noise modeling later on
sampleNames = unique(datSum$SampleName)
locNames = unique(datSum$Locus)

#Find MLE of Marker efficency (takes ~15 min):
markerEfficiency = inferMarkerEfficiency(dat=datSum,AT,verbose = TRUE)#,runMCMC = FALSE) #fitting model MLE model
Amhat = markerEfficiency$MLE$theta2[locNames] #estimates

#show estimated params
print(round(Amhat,3))


if(0) { #Example of showing emprical calculated vs model fit:
  markerEff = getMarkerEffStat(markerEfficiency)
  markerEff = markerEff[match(locNames,markerEff$Marker),]
  plot(markerEff$MEAN,Amhat,xlab="Empirical mean",ylab="Estimated",main="Comparison of amplification efficiency estimates");abline(a=0,b=1)
  
  #Compare Param with empirical observations:
  Amhat = markerEfficiency$MLE$theta2[locNames]
  AmSD = sqrt(diag(markerEfficiency$MLE$Sigma_theta2)[locNames])
  AmSD2 = markerEff$SD[ match(markerEff$Marker,locNames) ]
  
  expectedMarkers = c(t(replicate(length(markerEfficiency$sampleNames),markerEfficiency$locNames))) #obtain expected markers
  expectedSamples = rep(markerEfficiency$sampleNames,length(markerEfficiency$locNames)) #obtain expected markers
  df = data.frame(Sample=expectedSamples,Marker=factor(expectedMarkers,levels=locNames),Read=markerEfficiency$sumCovObs)
  library(dplyr);library(ggplot2)
  df <- df  %>% group_by(Sample)  %>% mutate( markerEff = Read/mean(Read) )
  
  gg <- ggplot(df, aes(x=Marker, y=markerEff)) +  geom_boxplot()
  gg <- gg + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  gg <- gg + ggtitle("Marker Efficiency (Empirical vs calibration)")
  for(i in seq_along(Amhat))  {
    gg <- gg + geom_segment(x=i, y=Amhat[i]-2*AmSD2[i], xend=i, yend=Amhat[i]+2*AmSD2[i],col=3,lwd=2)
    gg <- gg + geom_segment(x=i, y=Amhat[i]-2*AmSD[i], xend=i, yend=Amhat[i]+2*AmSD[i],col=2,lwd=2)
  }
  gg <- gg + geom_hline(yintercept=1,linetype="dashed") + ylim(0,6)
  print(gg)
}


##############################
#2) OBTAIN STUTTER MODEL DATA#
##############################

stutterModTab = read.table(datfn_stuttermodel, header=TRUE,sep=";")
stutterModList = list()
for(i in seq_len(nrow(stutterModTab))) {
  row = stutterModTab[i,]
  loc = row$Marker #obtain locus name
  type = row$Shortname #obtain stutter type
  est =  as.numeric(row[5:7])    #obtain estimats
  if( !loc%in%names(stutterModList)) stutterModList[[loc]] = list()
  stutterModList[[loc]][[type]] = as.numeric(na.omit(est)) #nicer to remove the NAs
}


###############################
#3) CALIBRATION OF NOISE MODEL#
###############################
#Infer noise model based on inferred stutter types
#sampleNames=unique(dat$SampleName); locNames=unique(dat$Locus); noisedata=stutterData$noisedata
datNoise1 = tableReader(datfn_noise) #non-stutters
datNoise1 = subset(datNoise1,select=c(Sample,Locus,Reads))
datNoise2 = tableReader(datfn_stutterFew) #stutters, but too few
datNoise2 = subset(datNoise2,select=c(Sample,Locus,Reads))
colnames(datNoise1) <-colnames(datNoise2) <- c("SampleName","Locus","Coverage")  #need correct header names
datNoise = rbind(datNoise1,datNoise2) #combine noise and rare stutters

noiseMod = inferNoiseModel(sampleNames=sampleNames, locNames=locNames, noisedata=datNoise, AT=AT)#,createPlots=TRUE)

#show estimated params
noisetab = list2DF(noiseMod)
rownames(noisetab) = locNames
print(round(noisetab,3))


#PREPARE CALIBRATION OBJECT (DEFINED PER MARKER)
calibration = list()
for(loc in locNames) {
  markerEff= Amhat[loc] #setNames(c(,AmSD[loc],AmSD[loc]),c("est","SD","SD2"))
  noiseMods = lapply(noiseMod,function(x) x[loc])
  stuttMods = stutterModList[[loc]] 

  calibration[[loc]] = list(markerEff=markerEff)  
  calibration[[loc]] = append(calibration[[loc]],noiseMods)

  #put stutter model last (requirement!!)
  calibration[[loc]] = append(calibration[[loc]],stuttMods)
}

fn_calibration = paste0(calibfile,".RDS")
saveRDS(calibration,file=fn_calibration)

#Calibration complete

