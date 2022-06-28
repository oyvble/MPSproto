#' @title fitSUMmodel
#' @description A function to fit sumPH of the observed signal intensities
#' @details This function is used for fitting a regression model (Gamma or Negative binomial)
#'
#' @param y Vector with sum peak heights (per loci)
#' @param x Vector with fragment lengths (base pair) (only provided for degradation)
#' @param niter Number of random samples 
#' @param delta Standard deviation of normal distribution when drawing random startpoints. Default is 1.
#' @param restrictDeg Whether to restrict degradation param to be less than 1
#' @param model Selected model for the signals (read/peak heights): {"GA"=gamma,"NB"=negative binomial}
#' @return Estimated parameters (mean,coefVar,slopePara)
#' @examples
#'\dontrun{
#' th = c(1000,0.3,0.6)
#' n = 100 #number of samples
#' x = seq(10,300,l=n)
#' x = (x-125)/100
#' y = rgamma(n,shape=2/th[2]^2*th[3]^x,scale=th[1]*th[2]^2)
#' plot(x,y)
#' fitSUMmodel(y,x,model="GA")
#' fitSUMmodel(y,x,model="NB")
#'}
#' @export

#x=NULL;niter=10;delta=1;restrictDeg=FALSE
fitSUMmodel <- function(y,x=NULL,niter=10,delta=1,restrictDeg=FALSE,model="GA") {
   DEG = TRUE #degradation used by default
   if(is.null(x)) DEG=FALSE #degradation set to FALSE if x not given
   isOK = !is.na(y) #get non-NA values
   y = y[isOK] #remove NAs
   if(DEG) x = x[isOK] #remove NAs
   
   #plot(x,y)
   #Impute small Y-vals for zero observations (more robust for low-template profiles)
   minY = min(y[y>0]) #get minimum observed (non-zero)
   y[y==0] = minY/2
   
   #define likelihood function depending on model:
   if(model=="GA") {
     loglikfun = function(th1) {
       shapev = 2/th1[2]^2 #expression without degradation
       if(DEG) shapev = shapev*th1[3]^x #include DEG part
       dgamma(y,shape=shapev,scale=th1[1]*th1[2]^2,log=TRUE)
     }
   } else if(model=="NB") {
     loglikfun = function(th1) {
       muv = 2*th1[1] #expression without degradation
       if(DEG) muv = muv*th1[3]^x #include DEG part
       size = th1[1]/(th1[1]*th1[2]^2-1) #obtain size from (mu,cv)
       dnbinom(as.integer(y),size=size,mu=muv,log=TRUE)
     }
   } else {
     stop("Model was not recognized!")
   }
   
   #FUNCTION TO BE MINIMIZED
   negloglik <- function(th) {
    th <- exp(th) #assume log-transformed
    val <- -sum(loglikfun(th))
    if( is.infinite(val)) val <- Inf #.Machine$integer.max 
    return(val)
   }
   
   #Helpfunction for optmization (may throw error)
   helpOptimize = function(th) {
      opt = list(min=Inf)
      suppressWarnings({
         tryCatch({ opt = nlm(negloglik, th )}
         , error = function(e) e )
      })
      return(opt)
   }
   
   #Obtain good start values (depending if degradation model or not)
   if(DEG) { #if degradation 
     coefs = exp(coef(lm(log(y)~x))) #prefit using logged values (convert params back)
     mu0 = coefs[1]/2 #heterozugout allele
     deg0 = coefs[2] #degrad slope
     sd0 = sqrt(mean((2*mu0*deg0^x - y)^2)) #recalculating sigma
   } else {
     mu0 = mean(y/2)
     sd0 = sd(y/2)
   }
   omega0 = sd0/mu0*0.9
    #"NB" = mu0/( sd0^2/mu0 -1 )*0.9
   
   th0 <- c(mu0,omega0)
   if(DEG) th0 <- c(th0,deg0)
   #if(verbose) cat(paste0("theta0=",paste0(th0,collapse=",")))  
   
   #Provide repeated optimizations
   largeVal = 1e300 #a large value (below inf)
   t0 = log(th0) #log transform pre-estimates
   bestFoo = helpOptimize(t0)
   if(niter==1 && bestFoo$minimum>largeVal) niter=30 #if 1st optimization failed
     cc <- 1 #
     while(cc<niter) { #repeat until accepted fitted
      t1 <- rnorm(length(t0),mean=t0,sd=delta) #generate new
      foo <- helpOptimize(t1)
      if(foo$min<bestFoo$min) bestFoo <- foo #if better    
      cc <- cc + 1 #add
     }
   if(bestFoo$minimum>largeVal) return(NULL) #return NULL if optimim failed
   
  th = exp(bestFoo$est) #get estimated parameters
  if(restrictDeg && DEG && th[3]>1 ) th[3] = 0.999 #set start value of degrad close to 1 (IMPORTANT TO AVOID CRASH IN contLikMLE)
  return(th)
} 

