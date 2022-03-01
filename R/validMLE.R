#' @title validMLE
#' @description Performing model validation of whether the observed peak heights fits the maximum likelihood fitted gamma distribution.
#' @details The cumulative probability of the observed allele peaks are calculated and compared with a uniform distribution.
#' Function calls a 'cumulative procedure' in C++ 
#'
#' @param mlefit Fitted object using inferEvidence function.
#' @param plottitle Maintitle text used in the PP-plot
#' @param alpha The significance level used for the envelope test. Default is 0.01
#' @param createplot Boolean of whether plot should be created
#' @param verbose Boolean of whether printing out information about significant alleles
#' @return retinfo A dataframe with information from the calculated validation model
#' @export

#plottitle="PP-plot";alpha=0.01;createplot=TRUE;verbose=TRUE
validMLE <- function(mlefit,plottitle="PP-plot",alpha=0.01,createplot=TRUE,verbose=TRUE) {
  cols = c("black","#df462a","#466791","#60bf37","#953ada","#5a51dc","#d49f36","#507f2d","#db37aa","#5b83db","#c76c2d","#552095","#82702d","#55baad","#62aad3","#8c3025","#417d61")
  xlab="Expected probabilities"#: Unif(0,1)"
  ylab="Observed probabilities"#: (Pr(Yj<=yj|Y_{-j}<=y_{-j},Yj>=thresh,model))"
  sz <- 1.5

  pchmax = 26 #number of possible point types
  c <- mlefit$prepareC #Object already stored in mlefit. returned from prepareC function
  par = mlefit$fit$par #estimated params
  theta_Am = c$markerEfficiency
  
  if(any(is.na(unlist(par)))) return(NULL) #return if params not found
  
  #Obtaining number of alleles to obtain an upper peak height boundary to integrate up to (maxY)
  nAtotObs = sum(c$peakHeights>0) #number of observations
  alphaQ <- 0.001 #ensure very far out in quantile (used for estimating probs in gamma-distribution).
  alpha2 <- alphaQ/nAtotObs #"bonferroni outlier"
  suppressWarnings({ #don't show missing allele warning
    maxYobs <- 1.5*max(c$peakHeights) #max observation ++
    maxYexp <- qgamma(1-alpha2,2*theta_Am/par$PHvar^2,scale=par$PHexp*par$PHvar^2) #max observation in theory (taking into account marker efficiency)
    maxY <- ceiling(max( na.omit(c(maxYobs,maxYexp)))) #get max observed
  }) 
  
  #PArt 1/2: Consider evaluation at PH-levels  
  noiseSizeWeight1 <- noiseSizeWeight2 <- 0*c$noiseSizeWeight #modify dropin weights to account for cdf instead of pdf
  for(locind in seq_len(c$nLocs)) { #traverse each marker
    AT0 = c$AT[ locind ]
    
    for(rind in seq_len(c$nSamples[locind]) ) { #for eac replicate
      for(aind in seq_len(c$nAlleles[locind]) ) {
        cind = c$startIndMarker_nAllelesReps[locind] + (aind-1)*c$nSamples[locind] + rind
        peak = c$peakHeights[cind] #obtain replicate
        #if(peak==AT0) peak = AT0 + 1 #put PH one increment above (avoid zero p-value)
        if(peak < AT0) next #skip if less than zero
        #freq = c$freq[ c$startIndMarker_nAlleles[locind] + aind] #obtain frequency (not used)
        #cumPareto = function(x,AT,alpha) 1-exp( alpha*(log(AT)-log(x)) ) #not used
        
        #obtain cumulative weight for observed PH (and maxY, which is assumed to be 1)
        noiseSizeWeight1[cind] = c$noiseSizeWeightCDF[cind] #log(cdf_noiseCov[getClosest(peak)]) #pexp(peak-AT0,noiseSizeParam0,log.p=TRUE) #= log( 1 - exp(- (lambda0)*(peak-AT0)))
        #noiseSizeWeight2[cind] = 0 #+ log(cdf_noiseCov[getClosest(maxY)]) #pexp(maxY-AT0,noiseSizeParam0,log.p=TRUE) #= log(freq) + 1
	   }
    }
  }
  #max(abs(noiseSizeWeight1-noiseSizeWeight2))
  
  #obtain cumulative vals int_0^y 
  pvalVEC = c$peakHeights #must have a copy
  UaPH = .C("loglikPrediction_cumprob",as.numeric(pvalVEC), as.numeric(0), as.integer(c$nJointCombs), c$NOC, c$NOK,
            as.numeric(par$mx),  as.numeric(par$mu), as.numeric(par$omega), as.numeric(par$beta), as.numeric(theta_Am),
            c$AT,c$fst,c$nNoiseParam,c$noiseSizeWeight, as.numeric(noiseSizeWeight1),
            c$nLocs, c$nSamples, c$nAlleles, c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps,
            c$peakHeights, c$freq, c$nTyped, c$maTyped, c$basepairs,
            c$nGenos, c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele, c$startIndMarker_outG1contr, 
            c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttExp , c$startIndMarker_nStutters,
            c$knownGind )[[1]] #obtain 
  #print(UaPH)
  
  #PArt 2/2: Consider evaluation at PH-levels  
  #obtain cumulative vals int_0^inf 
  pvalVEC = c$peakHeights #must have a copy
  UaMAX = .C("loglikPrediction_cumprob",as.numeric(pvalVEC), as.numeric(maxYobs), as.integer(c$nJointCombs), c$NOC, c$NOK,
             as.numeric(par$mx),  as.numeric(par$mu), as.numeric(par$omega), as.numeric(par$beta), as.numeric(theta_Am),
             c$AT,c$fst,c$nNoiseParam,c$noiseSizeWeight, as.numeric(noiseSizeWeight2),
             c$nLocs, c$nSamples, c$nAlleles, c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps,
             c$peakHeights, c$freq, c$nTyped, c$maTyped, c$basepairs,
             c$nGenos, c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele, c$startIndMarker_outG1contr, 
             c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttExp , c$startIndMarker_nStutters,
             c$knownGind )[[1]] #obtain 
  
  #print(UaMAX)
  #Combing the product from the two parts:  
  cumprobi = UaPH/UaMAX #obtaining cumulative probs
  #print(cumprobi)
  #hist(cumprobi)
  #Obtain names:
  locNames = c$locNames #obtain locus names
  alleleNames = c$alleleNames
  repNames = c$sampleNames
  
  tab = NULL
  for(locind in seq_len(c$nLocs)) { #traverse each marker
    loc0 = locNames[locind] #obtain locus name
    for(rind in seq_len(c$nSamples[locind])) { 
      rep0 =  repNames[rind] #obtain replicate name (NB: THIS WOULD LEAD TO WRONG NAME IF MISSING REPLICATE FOR A MARKER!)
      for(aind in 1:c$nAlleles[locind]) {
        cind = c$startIndMarker_nAllelesReps[locind] + (aind-1)*c$nSamples[locind] + rind
        peak = c$peakHeights[cind]
        if(peak==0) next 
        allele = c$alleleNames[[loc0]][aind] #obtain allele name
        new = c( loc0, rep0, allele, peak, cumprobi[cind])
        tab = rbind(tab, new)
      }
    }
  }
  colnames(tab) = c("Marker","Replicate","Allele","Coverage","cumProb")
  #nrow(tab)==sum(!is.nan(cumprobi))
  tab = as.data.frame(tab)
  tab$cumProb = as.numeric(tab$cumProb)
  
  cumprobi = as.numeric(tab[,ncol(tab)]) #obtain values again
  N <- nrow(tab) #length(cumprobi) #number of peak heights
  alpha2 <- alpha/N #0.05 significanse level with bonferroni correction
  cumunif <- ((1:N)-0.5)/N #=punif((1:N)-0.5,0,N)
  
  #hist(cumprobi)
  #sum(cumprobi<=0.5)/N
  ord <- order(cumprobi,decreasing = FALSE)
  ord2 <- match(1:length(ord),ord) #get reverse index
  pval <- 1-pbeta(cumprobi[ord],N*cumunif,N-N*cumunif+1) #one sided p-value
  
  #Must indicate the values below the line (Symmetry)
  ind <- cumprobi[ord]<cumunif #those below the line (two-sided p-value)
  pval[ind] <- pbeta(cumprobi[ord],N*cumunif,N-N*cumunif+1)[ind] #calculate pval in opposite direction
  #cumprobi<qbeta(alpha/2,N*cumunif,N-N*cumunif+1) | cumprobi<qbeta(1-alpha/2,N*cumunif,N-N*cumunif+1)
  outside <- pval[ord2] < (alpha2/2) #criterion outside region (divide by 2 to get two-sided)
  tab = cbind(tab,pvalue=pval[ord2],Significant=outside) #update table
  
  if(verbose) { #print info:
    #print points outside the bonferroni-adjusted envolopment
    print(paste0("Total number of peak height observation=",N))
    print(paste0("Significance level=",signif(alpha*100,digits=2),"%"))
    print(paste0("Bonferroni-adjusted significance level=",signif(alpha2*100,digits=2),"%"))
    print("List of observations outside the envelope (with the Bonferroni-level):") #two sided check
    print(tab[tab$Significant,],drop=FALSE)  #print list of outliers (outside envelop)
  }
  
  #plot   
  if(createplot) {
    zones <- matrix(c(1,1,1, 2,5,4, 0,3,0), ncol = 3, byrow = TRUE)
    layout(zones, widths=c(0.3,7,1), heights = c(1,7,.75))
    par(xaxt="n", yaxt="n",bty="n",  mar = c(0,0,0,0))   # for all three titles: 
    plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))   # fig 1 from the layout
    text(0,0,paste(plottitle), cex=2)  
    plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))  # fig 2
    text(0,0,paste(ylab), cex=2, srt=90)   
    plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1)) # fig 3
    text(0,0,paste(xlab), cex=2)  
    
    # fig 4: Rigth margin shows distr per replicate
    par(mar = c(1,0,0,0))
    
    plot(0,0,xlim=0:1,ylim=0:1,ty="n")
    for(rep in repNames) {
#rep=repNames[1]
      rind = which(repNames==rep)
      ind = tab$Replicate==rep #obtain index of certain replicate
      subDat = tab[ind,]
      locind = match(subDat$Marker,locNames) #obtain loci to plot
      
      points( rep(rind/(length(repNames)+1),sum(ind)), subDat$cumProb ,pch=(locind-1)%%pchmax,col=cols[rind],cex=sz)
    }
    rect(0,0,1,1)
    
    # fig 5, finally, the scatterplot-- needs regular axes, different margins
    par(mar = c(1,2,0,.5), xaxt="s", yaxt="s", bty="n")
    
    plot(0,0,ty="n",xlim=0:1,ylim=c(0,1),cex.axis=sz,asp=1)
    segments(x0=0,y0=0,x1=1,y1=1,lwd=1.2)
    
    #Goodness of fit test  #REMOVED: pval <- ks.test(cumprobi, "punif")$p.value
    #Updated block in 0.6.2: DRAW ENVELOPE LINES FOR ORDER STATISTICS UNDER RANDOMNESS:
    xsq <- seq(0,1,l=1000) #Draw envolope lines
    ysq <- c(alpha,alpha2) #quantiles to consider
    for(qq in ysq) { #for each quanitle
#      qq=ysq[1]
      lines(xsq,qbeta(qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
      lines(xsq,qbeta(1-qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
    }
    legend("bottomright",legend=paste0("1-Envelope-coverage=",c("",paste0(alpha,"/",N,"=")),signif(ysq,2)),col=1:length(ysq),lty=2,cex=1.3)
    
    #  abline(0,1)
    locind = match(tab$Marker,locNames) #obtain loci to plot
    repind = match(tab$Replicate,repNames) #get replicate indices
    points(cumunif,tab$cumProb[ord],pch=(locind[ord]-1)%%pchmax,cex=sz,col=cols[repind[ord]])
    
    legend("topleft",legend=repNames,pch=19,cex=sz,col=cols[1:length(repNames)])
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    par(op)
  } #end createplot
  
  #return information
  return(tab)
  
}

