#' @title diagnoseMCMC 
#' @description Helpfunction to obtain MCMC diagnostic
#' @param posttheta Matrix with posterior samples (n x p)
#' @param showTrace Whether a trace plot of samples should be shown
#' @param xlim x-axis limits
#' @param MEAN Mean of superimposed distribution 
#' @param SD SD of superimposed distribution 
#' @param fitBeta Whether a Beta instead of normal distribution should be superimposed
#' @export 

diagnoseMCMC = function(posttheta,showTrace=FALSE,xlim=NULL,MEAN=NULL,SD=NULL,fitBeta=FALSE) {
  npar = ncol(posttheta)
  parname = colnames(posttheta)
  par(mfrow=c(npar,1+sum(showTrace)))
  npoints = 10000
  for(i in 1:npar) {
    
    if(!is.null(xlim)) {
      d = density(posttheta[,i],from=xlim[1],to = xlim[2],n = npoints)
    } else {
      ispos = all(posttheta[,i]>=0)
      if(ispos) {
        d = density(posttheta[,i],from=0,n = npoints)
      } else {
        d = density(posttheta[,i],n = npoints)
      }
    }
      
    plot(d,main=parname[i],xlab="Param",ylab="Density",lwd=2)
    if(!is.null(MEAN) && !is.null(SD)) {
      lines(d$x,dnorm(d$x,MEAN[i],SD[i]),col=2,lty=2,lwd=1)
      s0 = 2 #signif
      
      leg = paste0("N(",signif(MEAN[i],s0),",",signif(SD[i],s0),")")
      col = 2
        
      if(fitBeta) { #wheter to fit Beta
        tmp = MEAN[i]*(1-MEAN[i])/SD[i]^2 - 1
        sh1 = MEAN[i]*tmp
        sh2 = (1-MEAN[i])*tmp
        lines(d$x,dbeta(d$x,sh1,sh2),col=3,lty=2,lwd=1)
        
        leg = c(leg,paste0("B(",signif(sh1,s0),",",signif(sh2,s0),")")) #add to legend
        col = c(col,3)
      }
      legend("topright",leg,lty=2,col=col)
      
    }
    if(showTrace) plot(posttheta[,i],ty="l",main=parname[i],ylab="Param")
  }
  #par(mfrow=c(npar))
}