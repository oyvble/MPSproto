#' @title inferEvidence2
#' @description Function to perform inference of new evidence (Penalized MLE based) 
#' @details Optimizes sample specific params (mx,mu,omega, Am). Using MLE as initial params for optimization. Marker efficiency parameter is penalized on prior.

#' @param mlefit Output object from inferEvidence function
#' @param calibration A list from calibration results (per marker). Must contain the following elements (in following order):
#' AT Analytical threshold 
#' markerEff Marker efficiency param (possibly also SD if running inferEvidence2)
#' @param priorAsLN Whether model prior of marker efficiency as a log-normal (SD from logged values)
#' @param reOpt Whether to re-run optimization with random start values
#' @param nDone Number of optimizations required providing equivalent results (same logLik value obtained)
#' @param delta Scaling of variation of normal distribution when drawing random startpoints. Default is 1.
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @param seed The user can set seed if wanted
#' @param verbose Whether to print out progress
#' @export 

inferEvidence2 = function(mlefit, calibration=NULL, priorAsLN=FALSE, reOpt=FALSE, nDone=3, delta=1,steptol=1e-4, seed=NULL, verbose=FALSE) {
  if(!is.null(seed)) set.seed(seed) #set seed if provided
  
  #Obtain already prepared C-dataobject 
  c_obj = mlefit$prepareC
  hypothesis = mlefit$hypothesis
  
  if(is.null(calibration)) calibration = mlefit$calibration
  
  #construct prior-object for marker efficiency
  prior_MEAN =  sapply(calibration,function(x) as.numeric(x$markerEff[1]) )
  prior_SD =  sapply(calibration,function(x) as.numeric(x$markerEff[2]) )

  #Be sure to get correct order:
  prior_MEAN = prior_MEAN[c_obj$locNames]
  if(priorAsLN)  prior_MEAN =  log(prior_MEAN) #take log if want log-normal assumption
  prior_SD = prior_SD[c_obj$locNames]
  prior_SD1 = head(prior_SD,-1) #remove last 
  
  #Obtain fitted params from previous model fit (model 1) as initial value 
  phihat0 = mlefit$fit$phihat #obtain previous fitted ordinary vals
  phihat1 = head(log(prior_MEAN),-1) #
  
  #PREPEARING THE LIKELIHOOD OPTIMZATION (what parameters are provided?):
  NOC = hypothesis$NOC
  condOrder = hypothesis$cond
  nparam = length(phihat0) + length(phihat1)  #number of params to estimate

  modelDeg = length(phihat0)==(NOC+2) #wheter modeling degradation (beta is last param)
      
  #INNER HELPFUNCTIONS:
  logit = function(x) log(x/(1-x))
  drawMixprop = function(NOC) {
    #CONSIDER Mixture porpotions    
    ncond = sum(condOrder>0)
    nU = NOC - sum(condOrder>0)
    mxrnd = rgamma(NOC,1) #Draw simplex (flat)
    mxrnd = mxrnd/sum(mxrnd)
    if(nU>1) { #sort if more than 1 unknown
      ind = (ncond+1):NOC #sort Mix-prop for the unknowns 
      mxrnd[ind] = sort(mxrnd[ind],decreasing = TRUE) #sort Mx in decreasing order
    }
    
    #convert Mx values to real domain (nu:
    nurnd = numeric()
    if(NOC>1) {
      cs = 0 #c( 0,cumsum(mxrnd)) #cumulative sum of mixture proportins
      for(cc in 1:(NOC-1)) { #traverse contributors (Restricted)
        nurnd = c(nurnd, logit( mxrnd[cc]/(1-cs))) 
        cs = cs + mxrnd[cc] #update sum
      }
    }
    return(nurnd)
  }
  
  drawPHvars = function() {   #CONSIDER PH prop variables
    th0 = exp(na.omit(phihat0[NOC + 0:2])) #remove possible NAs
    sdPH = delta*0.15*th0 #obtain considered SD of PH props 
    return( log( abs( rnorm(length(th0),th0,sd=sdPH)) ) )  #Obtain random start for mu/sigma/tau, Note using the delta here (should be small)
  }

  drawMarkerEff = function() {   #Draw from defined prior distr
    if(priorAsLN) {
      return( rnorm(length(phihat1),phihat1, prior_SD1) )  #already transformed
    } else {
      return( log(abs(rnorm(length(phihat1),exp(phihat1),prior_SD1) ) ) ) # transform
    }
  }

  #Helpfunction to convert transfered param back  (also includes marker efficiency Am)
  convParamBack = function(phi) {
    indMxParam = numeric()
    if(NOC>1) indMxParam = 1:(NOC - 1) # number of param
    
    indMuOmegaBetaParam = (NOC - 1)  + 1:3 #assume common mu,omega param (each replicates)
    if(!modelDeg) indMuOmegaBetaParam = head(indMuOmegaBetaParam,-1) #remove last index
    
    #obtain params:
    theta_mx = phi[indMxParam]
    theta_mu = phi[indMuOmegaBetaParam[1]]
    theta_omega = phi[indMuOmegaBetaParam[2]]
    if(modelDeg) {
      theta_beta = phi[indMuOmegaBetaParam[3]]
    } else {
      theta_beta = 0 #set a number if not modeling degradation      
    }
    
    #Mixture proportions
    mx_vec = c(theta_mx,1) #obtain mixture proportions (default)
    theta_mu = exp(theta_mu)  
    theta_omega = exp(theta_omega)
    theta_beta = exp(theta_beta)

    theta_Am = exp( phi[-c(indMxParam,indMuOmegaBetaParam)] ) #invert back
    theta_Am = c(theta_Am, length(theta_Am)+ 1 - sum(theta_Am) )
    #mean(theta_Am)==1
    
    names(theta_Am) = c_obj$locNames
    #Transform mix prop:
    if(NOC>1) {
      cs = 0.0; #init cumulative sum
      for(i in 1:(NOC-1)) { #for each mix props
        mx_vec[i] = (1-cs)/(1+exp(-theta_mx[i])) #convert back (ilogit)
        cs = cs +  mx_vec[i]; #add mixture proportion to cumulative sum 
      } 
      mx_vec[NOC] = 1-sum(mx_vec[-NOC]) #last contributor is 1- sum of mix-prop
    }
    return(list(mx=mx_vec,mu=theta_mu,omega=theta_omega,beta=theta_beta,Am=theta_Am))
  }
  
  #function for calling on C-function: Must convert "real domain" values back to model params 
  negloglik_phi <- function(phi) { #assumed order: mixprop(1:C-1),mu,sigma,beta,xi
    parList = convParamBack(phi) #obtain list with parameters (on original scale)
    
    if(any(parList$Am<0)) return(Inf)
    
    #Update C-object with marker efficiencies:
    c_obj$markerEfficiency = parList$Am
    
    logLik = calcLogLikC_prediction(parList,c_obj)
	
	  if(priorAsLN) parList$Am = log(parList$Am) #simply log the values (notice that prior_MEAN,prior_SD must also be on logged values
	
	  log_prior = sum(dnorm(parList$Am, prior_MEAN,prior_SD,log = TRUE)) #obtain prior weight (normal assumption)	
	
    #print(calcObj$logLik)
    return( -(logLik+log_prior) ) #return negative log-likelihood
  }

  
  #Iterate with several startvalue (mx is random)
  nOK <- 0 #number of times for reaching largest previously seen optimum
  maxITERS <- 30 #number of possible times to be INF or not valid optimum before any acceptance
  nITER <- 0 #number of times beeing INF (invalid value)
  
  maxL <- -Inf #value of maximum obtained loglik
  maxPhi <- rep(NA,nparam) #Set as NA
  maxSigma <- matrix(NA,nparam,nparam)#Set as NA
  
  logLik_tolerance = 0.01 #tolerance of accepting similar liklihood optimization
  if(verbose) print(paste0("Number of parameters to optimize: ",nparam))
  suppressWarnings({
    while(nOK<nDone) {
      
      #FIRST: generate start values (real domain): 
      if(reOpt) {
        phi0a = drawMixprop(NOC) #draw mixture proportions
        phi0b = drawPHvars()
        phi0c = drawMarkerEff() #propose new
        phi0 = c(phi0a,phi0b,phi0c  )#append with the other params
        # convParamBack(phi0)
      } else {
        phi0 = c(phihat0,phihat1) #use already fitted values as startvals
      }
      
      #Check if proposal is OK
      timeOneCall = system.time({ #estimate the time for calling the likelihood function one time 
        likval <- -negloglik_phi(phi=phi0)   #check if start value was accepted
        if(verbose && !reOpt) print(paste0("Prev. maximum at loglik=",likval))
      })[3] #obtain time in seconds
      
      if( is.infinite(likval) ) { #if it was infinite (invalid)
        nITER = nITER + 1	 
        if(!reOpt) nITER = maxITERS #set as done with optim if not re-optimizing
      } else { #PERFORM OPTIMIZATION
        
        tryCatch( {
          foo <- nlm(f=negloglik_phi, p=phi0, iterlim=1000,steptol=steptol, hessian=TRUE)#,print.level=2)
          Sigma <- solve(foo$hessian)
          
          if(all(diag(Sigma)>=0) && foo$iterations>2) { #} && foo$code%in%c(1,2)) { #REQUIREMENT FOR BEING ACCEPTED
            nITER <- 0 #reset INF if accepted
            likval <- -foo$min #obtain local maximum
            
            #was the maximum (approx) equal the prev: Using decimal numbers as difference 
            isEqual = !is.infinite(maxL) && abs(likval-maxL) < logLik_tolerance #all.equal(likval,maxL, tolerance = 1e-2) # # #was the maximum (approx) equal the prev?
            
            if(isEqual) { 
              if(verbose)  print(paste0("Equal maximum found: loglik=",likval))
              nOK = nOK + 1 #add counter by 1
            } else {  #if values were different
              if(likval>maxL) { #if new value is better
                nOK = 1 #first accepted optimization found
                maxL <- likval #maximized likelihood
                maxPhi <- foo$est #set as topfoo     
                maxSigma <- Sigma 
                if(verbose) print(paste0("New maximum at loglik=",likval))
              } else {
                if(verbose)  print(paste0("Local (non-global) maximum found at logLik=",likval))
              }
            } 
            if(verbose && reOpt) print(paste0(" (",nOK,"/",nDone,") optimizations done"))
            #flush.console()
          } else { #NOT ACCEPTED
            nITER <- nITER + 1 
          }
        },error=function(e) e,finally = {nITER <- nITER + 1} ) #end trycatch (update counter)
      } #end if else 
      
      if(!reOpt || (nOK==0 && nITER>maxITERS)) {
        nOK <- nDone #finish loop (also if  not re-optimizing)
        break #stop loop if too many iterations  
      }
    } #end while loop
  })

  #Convert param back:
  par = convParamBack(phi=maxPhi) #obtain list with parameters (on original scale)
  
  #Remember to subtract with prior for likelihood value (IMPORTANT)
  Amhat = par$Am #estimated Am
  if(priorAsLN) Amhat = log(Amhat)
  log_prior = sum(dnorm(Amhat, prior_MEAN,prior_SD,log = TRUE)) #obtain prior weight (normal assumption)	
  logLik = maxL - log_prior #subtract with prior
  
	#ret <- mget(names(formals()),sys.frame(sys.nframe())) #return all values in argument call
  mlefit$fit = list(par=par,phihat=maxPhi,phiSigma=maxSigma,loglik=logLik,loglik2=maxL) #override param estimates
  mlefit$prepareC$markerEfficiency <- as.numeric(par$Am) #Also update C-object with marker efficiencies
  
  return( mlefit )
}
