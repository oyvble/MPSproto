#' @title getMarkerEffStat 
#' @description Obtain marker efficiency statistics (empirical)
#' @param markerEfficiency Object returned from inferMarkerEfficiency
#' @param dfAsInput If markerEfficiency is a dataframe input (optional)
#' @param locNames Locus name can be included to override order (optional)
#' @export 
getMarkerEffStat = function(markerEfficiency,dfAsInput=FALSE,locNames=NULL) {
  library(dplyr)
  if(!dfAsInput) {
    expectedMarkers = c(t(replicate(length(markerEfficiency$sampleNames),markerEfficiency$locNames))) #obtain expected markers
    expectedSamples = rep(markerEfficiency$sampleNames,length(markerEfficiency$locNames)) #obtain expected markers
    if(is.null(locNames)) {
      df = data.frame(Sample=expectedSamples,Marker=factor(expectedMarkers),Read=markerEfficiency$sumCovObs)
    } else {
      df = data.frame(Sample=expectedSamples,Marker=factor(expectedMarkers,levels=locNames),Read=markerEfficiency$sumCovObs)
    }
  } else {
	df = markerEfficiency
    colnames(df) = c("Sample","Marker","Read") #get correct header names
	if(!is.null(locNames)) df$Marker = factor(df$Marker,levels=locNames)
  }
  df$Read[df$Read==0] <- 1 #impute missing as 1 (avoid logged zero)
  df <- df  %>% dplyr::group_by(Sample)  %>% dplyr::mutate( markerEff = Read/mean(Read) )
  stat = df  %>% dplyr::group_by(Marker) %>% dplyr::summarise(MEAN=mean(markerEff),SD=sd(markerEff), MEAN2=mean(log(markerEff)),SD2=sd(log(markerEff)) )
  return(stat)
}