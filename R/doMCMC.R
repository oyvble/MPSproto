#' @title doMCMC 
#' @description Helpfunction to do MCMC
#' @param logLik Objective function to perform MCMC on
#' @param th0 Start value of parameters
#' @param Sigma0 Assumed covariance matrix of sampler
#' @param niter Number of samples to draw
#' @param delta Scaling covariance of the sampler (tweaking acceptance rate)
#' @export 

doMCMC = function(logLik,th0,Sigma0,niter=NULL,delta=NULL) {
  if(is.null(niter)) niter=1000
  if(is.null(delta)) delta=1

  np = length(th0) #number of param
  loglik0 = logLik(th0) #obtain lik val
  
  C <- chol(delta*Sigma0) #scale variance with a factor 2: ensures broad posterior
  X <- matrix(rnorm(np*niter),ncol=np,nrow=niter)%*%C #proposal values
  
  logdmvnorm <- function(X2,mean,cholC) { #function taken from mvtnorm-package
    p <- nrow(cholC)
    tmp <- backsolve(cholC,X2-mean,p,transpose=TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(cholC))) - 0.5*p*log(2*pi) - 0.5*rss
    return(logretval)
  }
  
  #Performing MCMC Metropolis Hastings by Gelfand and Dey (1994), using h() = Normal(th0,delta*Sigma)
  postth <- matrix(NA,ncol=np,nrow=niter) #accepted th
  postlogL <- rep(NA,niter) #accepted th
  postth[1,] <- th0 #init with proper value X[1,] #add noise into start point
  postlogL[1] <- loglik0 #init with proper value logliktheta(postth[1,]) #get start-likelihood  
  U <- runif(niter) #random numbers
  m <- 2 #counter for samples
  nacc <- 0  
  while(m<=niter) {
    postth[m,] <-  X[m,] + postth[m-1,] #proposed theta
    postlogL[m] <- logLik(th=postth[m,])
    pr <- exp(postlogL[m]- postlogL[m-1]) #acceptance rate
    if(U[m]>pr) { #if not accepted, i.e. random prob too large (above pr)
      postth[m,] <-  postth[m-1,]
      postlogL[m] <- postlogL[m-1 ]
    } else {
      nacc <- nacc + 1
    }
    m <- m + 1 #update counter
  } #end while not done
  accrat <- nacc/(niter-1) #acceptance ratio (dont count first)
  #  logpX <- logdmvnorm(postth,mean=th0,cholC=chol(Sigma0)) #insert with Normal-approx of post-th
  logpX <- logdmvnorm(t(postth),mean=th0,cholC=C) #insert with Normal-approx of post-th
  
  #calculate marginalization
  logVals = logpX - postlogL #this is logged values intended to be "exped"
  offset = max(logVals) #find max value
  #margL <- 1/mean(exp(logVals - offset)) #estimated marginal likelihood
  logMargL <- log(length(logVals)) - offset - log(sum(exp(logVals-offset))) #estimated marginal likelihood (logged)
  colnames(postth) = names(th0) #insert param names
  
  return(list(logMargL=logMargL,posttheta=postth,accrat=accrat,postlogL=postlogL,logpX=logpX,th0=th0,Sigma0=Sigma0))
}