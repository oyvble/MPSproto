#' @title inferEvidence
#' @description Function to perform inference of new evidence (MLE based) restricted on calibration inference
#' @details Optimizes only sample specific params (mx,mu,omega). Fixating params from fitted calibration (marker efficiency, noise, stutter expectations). 
#Note that prepareData_prediction should not be applied before this function.

#' @param samples Sample information for evidence profile(s). Must be a list[[sample]][[loc]] = list(adata,hdata)
#' @param popFreq Allele frequencies for a given population. Must be a list[[loc]] = freqs
#' @param refData Sample information for reference profile(s). Must be a list[[sample]][[loc]] = adata
#' @param hypothesis A list which defines the hypothesis to evaluate. Must contain the following elements:
#' NOC Number of contributors in model (Mandatory).
#' condOrder (conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' fst The co-ancestry coefficient. Can be a vector (must contain the marker names). Default is 0.
#' knownRef Specify known non-contributing references from refData (indices). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param calibration A list from calibration results (per marker). Must contain the following elements (in following order):
#' AT Analytical threshold 
#' markerEff Marker efficiency param (possibly also SD if running inferEvidence2)
#' nNoiseParam Parameter for modeling distr of number of noise (Noise model)
#' noiseSizeParam  Parameter for modeling size of noise alleles (Noise model)
#' Reg-coeffs (b0,b1,b2) for different stutter types (each element considers a particular stutter type)
#' @param kit Selected kit for assuming degradation model.
#' @param platform Which platform that is used (MPS or CE)
#' @param nDone Number of optimizations required providing equivalent results (same logLik value obtained)
#' @param delta Scaling of variation of normal distribution when drawing random startpoints. Default is 1.
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @param seed The user can set seed if wanted
#' @param verbose Whether to print out progress
#' @param model Selected model for the signals (read/peak heights): {"GA"=gamma,"NB"=negative binomial}
#' @return Fitted object for evidence (similar as euroformix::contLikMLE)
#' @export
#' @examples
#' \dontrun{ 
#' kit = "ForenSeq"
#' pkg = path.package("MPSproto")
#' calib =  readRDS(paste0(pkg,"/paper_stutterChar/calibrated_MPSproto.RDS"))
#' popFreq =  importMPSfreqs(paste0(pkg,"/paper_stutterChar/freqFile_ForenSeqFWbrack_Norway.csv"))[[1]]
#' gen = genMPSevidence(calib,2,popFreq,mu=1000,omega=0.2,beta=1,kit=kit )
#' plotMPS(gen$samples,gen$refData,AT=10)
#' hyp = list(NOC=2,cond=c(1,0),fst=0.01)
#' fit = inferEvidence(gen$samples,popFreq,gen$refData, hyp,calib)
#' valid = validMLE(fit)
#' }

inferEvidence = function(samples, popFreq, refData=NULL, hypothesis, calibration,  kit=NULL,platform="MPS", nDone=3, delta=1, steptol=1e-4, seed=NULL, verbose=FALSE, model="GA") {
  if(!is.null(seed)) set.seed(seed) #set seed if provided
  
  #pre-step to get Q-assignation:
  AT = sapply(calibration,function(x) as.numeric(x$AT)) #obtain AT
  normalize = FALSE
  if(!is.null(hypothesis$normalize) && hypothesis$normalize) normalize = TRUE #assign true if indicated
  dat = prepareData_prediction(samples,refData, popFreq, AT=AT, normalize=normalize, minF=hypothesis$minF)
  
  #call first:
  c_obj = prepareC_prediction(dat, hypothesis,calibration, kit, platform, model)  #obtain data to be used
  
  #OBTAIN START VALUES FOR PARAMS:
  theta0 = prefit_prediction(c=c_obj) #obtain param start values across the replicates (mu,sigma,beta)
   
  #PREPEARING THE LIKELIHOOD OPTIMZATION (what parameters are provided?):
  NOC = hypothesis$NOC
  condOrder = hypothesis$cond
  nparam = NOC - 1 + length(theta0) #number of params to estimate
  
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
  
  drawPHvars = function(th0) {   #CONSIDER PH prop variables
    th0 = na.omit(th0) #remove possible NAs
    sdPH = delta*0.15*th0 #obtain considered SD of PH props 
    return( log( abs( rnorm(length(th0),th0,sd=sdPH)) ) )  #Obtain random start for mu/sigma/tau, Note using the delta here (should be small)
  }
  
  #Helpfunction to convert transfered param back 
  convParamBack = function(phi) {
    indMxParam = numeric()
    if(NOC>1) indMxParam = 1:(NOC - 1) # number of param
    indMuOmegaBetaParam = (NOC - 1)  + 1:3 #assume common mu,omega param (each replicates)
    
    #obtain params:
    theta_mx = phi[indMxParam]
    theta_mu = phi[indMuOmegaBetaParam[1]]
    theta_omega = phi[indMuOmegaBetaParam[2]]
    theta_beta = phi[indMuOmegaBetaParam[3]]
    theta_beta[is.na(theta_beta)] = 0 #set a number if not given
    
    #Mixture proportions
    mx_vec = c(theta_mx,1) #obtain mixture proportions (default)
    theta_mu = exp(theta_mu)  
    theta_omega = exp(theta_omega)
    theta_beta = exp(theta_beta)
    
    #Transform mix prop:
    if(NOC>1) {
      cs = 0.0; #init cumulative sum
      for(i in 1:(NOC-1)) { #for each mix props
        mx_vec[i] = (1-cs)/(1+exp(-theta_mx[i])) #convert back (ilogit)
        cs = cs +  mx_vec[i]; #add mixture proportion to cumulative sum 
      } 
      mx_vec[NOC] = 1-sum(mx_vec[-NOC]) #last contributor is 1- sum of mix-prop
    }
    return(list(mx=mx_vec,mu=theta_mu,omega=theta_omega,beta=theta_beta))
  }
  
  #function for calling on C-function: Must convert "real domain" values back to model params 
  negloglik_phi <- function(phi) { #assumed order: mixprop(1:C-1),mu,sigma,beta,xi
    parList = convParamBack(phi) #obtain list with parameters (on original scale)
    #write.table(unlist(parList),file="params.txt",append = TRUE) #USED FOR DEBUGGING NB-model ERROR
    logLik = calcLogLikC_prediction(parList,c_obj)
    #print(calcObj$logLik)
    return( -logLik ) #return negative log-likelihood
  }

  #Iterate with several startvalue (mx is random)
  nOK <- 0 #number of times for reaching largest previously seen optimum
  maxL <- -Inf #value of maximum obtained loglik
  maxITERS <- 100 #number of possible times to be INF or not valid optimum before any acceptance
  nITER <- 0 #number of times beeing INF (invalid value)
  
  logLik_tolerance = 0.01 #tolerance of accepting similar liklihood optimization
  if(verbose) print(paste0("Number of parameters to optimize: ",nparam))
  suppressWarnings({
    while(nOK<nDone) {
      
      #FIRST: generate start values (real domain): 
      phi0 = drawMixprop(NOC) #draw mixture proportions
      phi0 = append(phi0, drawPHvars( theta0 ) )#append with the other params
 # convParamBack(phi0)
      timeOneCall = system.time({ #estimate the time for calling the likelihood function one time 
        likval <- -negloglik_phi(phi=phi0)   #check if start value was accepted
      })[3] #obtain time in seconds
      
      if( is.infinite(likval) ) { #if it was infinite (invalid)
        nITER = nITER + 1	 
        
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
            if(verbose) print(paste0(" (",nOK,"/",nDone,") optimizations done"))
            #flush.console()
          } else { #NOT ACCEPTED
            nITER <- nITER + 1 
          }
        },error=function(e) e,finally = {nITER <- nITER + 1} ) #end trycatch (update counter)
      } #end if else 
      
      if(nOK==0 && nITER>maxITERS) {
        nOK <- nDone #finish loop
        maxL <- -Inf #maximized likelihood
        maxPhi <- rep(NA,nparam) #Set as NA
        maxSigma <- matrix(NA,nparam,nparam)#Set as NA
        break #stop loop if too many iterations  
      }
    } #end while loop
  })
  
  #Convert param back:
  par = convParamBack(maxPhi) #obtain list with parameters (on original scale)
  
  ret <- mget(names(formals()),sys.frame(sys.nframe())) #return all values in argument call
  ret$fit=list(par=par,phihat=maxPhi,phiSigma=maxSigma,loglik=maxL)
  ret$samples <- ret$popFreq <- ret$refData <- NULL #remove these
  ret$data = dat #include the processed data instead
  ret$prepareC = c_obj
  ret$model = model #copy model name
    
  #Return results and data:
  return( ret )
}
