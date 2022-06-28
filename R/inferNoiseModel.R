#' @title inferNoiseModel 
#' @description Fitting noise model in data
#' @details Fitting two parts of the noise model: 
#' 1) Fitting the number of noise-observations (per sample) to a geometric distr (1 param) 
#' 2) Fitting the coverage/reads to the pareto distribution
#' @param sampleNames The samples names used for calibrating the noise model (should also count situations without noise)
#' @param locNames The marker names used in calibration.
#' @param noisedata A table with potential noise data (not typical stutters). Can be obtained from the ObtainStutterData function (noisetable)
#' @param AT A detection threshold value or thresholds per marker. NULL ignores filtering.
#' @param print2pdf Name of pdf file for plots
#' @export 
inferNoiseModel = function(sampleNames, locNames, noisedata, AT, print2pdf=NULL) {
  sig0 = 3 #signif number in plot
  noiseCap = 1000 #this is theoretical maximum of noise
  
  #Check that all column names are included  
  if(!all(c("Locus","Sample","Coverage")%in%colnames(noisedata) )) stop("Missing column name!")
  
  #Go through each marker and fit model:
  createPlots = !is.null(print2pdf)
  if(createPlots) pdf(paste0(print2pdf,".pdf"))
  
  noiseFit = list()
  for(loc in locNames) {
#    loc = locNames[2]
#    print(loc)
    sub = subset(noisedata,noisedata$Locus==loc) #obtain noise subset data
    
    AT0 = AT[loc] #obtain Analytical threshold
    if(is.na(AT0)) {
      AT0 = AT[1] #obtain Analytical threshold
      if(length(AT)>1)  stop(paste0("AT was not found for marker ",loc))
    } 
    
    if( nrow(sub)>0 ) { #if observations 
       #Obtain number of noise alleles per sample:
      noise_samples = factor(sub$Sample,levels = sampleNames)
      nNoise = table(noise_samples) #obtain number of reads per sample
      #nNoise = table(sub$SampleNames, sub$Noise)[,2]
    } else {
      nNoise = rep(0,length(sampleNames)) #assign zero-vector for all samples
    }
    nNoise = c(nNoise,max(nNoise)+1)#add an extra noise
    phat = 1/(mean(nNoise)+1) #estimate paramer for geometrical parmam
        
    #Model # occurences
    nNoiseF = factor(nNoise,levels=seq(0,max(nNoise)))
    noiseCounts = table(nNoiseF) #aggregate number of noises per sample (for marker)
    freq = noiseCounts/sum(noiseCounts) #obtain frequency
    #xi0 = freq[1] #no drop-in prob
    #xi1 = 1-xi0 #more than one drop-in prob
    #C = xi1/(1+xi1) #obtain prob when more than one param
    
    #exp = sum(outf*freq) #obtain expectation 
    #outf1 = outf[-1]
    
    if(createPlots) {
      outf = as.integer(names(freq))
      outf = c(outf,max(outf)+1)
      freq2 = c(freq,0) 
      names(freq2) = outf

      dgeomVal =  dgeom(outf, prob = phat)
      barplot(freq2,space = 0,main=loc,xlab="Number of noise alleles per sample (k)",ylim=c(0,max(dgeomVal,freq)))
      lines(outf+0.5,dgeomVal ,col=1,ty="o",lwd=2)
      mtext("Distribution of noise number")
      
      phat0term = signif(phat,sig0)
      #phat1term = signif(1-phat,sig0)
      #txt = paste0("Pr(#noise=k)=",phat0term,"x",phat1term,"^k")
      txt = paste0("phat=",phat0term)
      legend("topright",txt,pch=15,col=1)
    }
    
    #Modeling distr of reads (assume continuous Pareto):
    sizeNoise = sub$Coverage #obtain reads
    sizeNoise = sizeNoise[sizeNoise>=AT0]
    sizeNoise = c(sizeNoise,AT0) #add one noise at threshold (always)
    
    #Exception when no noise is included: Adding two more noise levels (this is default model)
    needAdd = length(unique(sizeNoise))==1 #whether more noise must be added to fit model
    if(needAdd) sizeNoise = c(sizeNoise,AT0,AT0+1) #add more noise
    
    #Obtain MLE of continous pareto (use as start value for discrete prob)
    alphahat0 = 1/(mean(log(sizeNoise)) - log(AT0)) 

    #Discrete     
    negLogLik_paretoDiscrete = function(logalpha) {
      zv = dpareto(discrete_range,exp(logalpha),AT0)#,logged=TRUE
      zv = zv/sum(zv)
      return( -sum(log(zv[match(sizeNoise,discrete_range)])))
    }
    
    #maxCap = max(noiseCap,sizeNoise) #obtain max cap
    discrete_range = seq(AT0,noiseCap)
    
    suppressWarnings({
      foo = nlm(negLogLik_paretoDiscrete,log(alphahat0))
    })
    alphahat = exp(foo$estimate) #transfer back
    #exp( log(alpha) + alpha*log(AT) - (alpha+1)*log(x) )
    
    if(createPlots) {
      sizeNoise = factor(sizeNoise,levels=seq(AT0,max(sizeNoise)))
      freq = table(sizeNoise)/length(sizeNoise)
      
      outf = as.integer(names(freq))
      outf = c(outf,max(outf)+1)
      freq2 = c(freq,0) 
      names(freq2) = outf
      
      dparetoDiscrete = function(x,alpha,AT,discrete_range) {
        zv = dpareto(discrete_range,alpha,AT)#,logged=TRUE
        return( zv[match(x,discrete_range)] ) #return prob for selected vals
      } 
      zv = dparetoDiscrete(outf, alphahat,AT0,discrete_range)
      names(zv) = outf #insert  range
      
      barplot(freq2,space=0,ylim = c(0,max(zv,freq2)),main=loc,xlab="Coverage of noise alleles (all samples)")
      #plot(h$mids,h$density,main=loc,ty="h",ylim=c(0,max(h$density,zv)))
      lines(outf-AT0+0.5,zv,col=1,lty=1,lwd=2,ty="o")
      mtext("Distribution of noise coverage")
      alphahat0 = signif(alphahat,sig0) #estimated alpha
      txt1 = paste0("alphahat=",alphahat0)
      legend("topright",txt1,pch=15,col=1)
    }      
  
    noiseFit[[loc]] = list(noiseCounts=noiseCounts,thGeom=phat,thPareto=alphahat,AT=AT0)
  }
  if(createPlots) dev.off()
  
  #Rearrange the data storage:
  retList = list()
  retList$nNoiseParam = sapply(noiseFit,function(x)  as.numeric(x$thGeom))
  retList$noiseSizeParam = sapply(noiseFit,function(x)  x$thPareto)
  retList$AT = sapply(noiseFit,function(x) as.numeric(x$AT))
  
  return(retList)
}
