#' @title optimizeLik 
#' @description helpfunction to repeat optim steps with random start points
#' @param logLik Objective function to optimize 
#' @param phi0 Start value of parameters (real space)
#' @param maxIter Maximum trials of optimization
#' @param delta SD of randomly drawed start values
#' @param steptol Step tolerance in nlm
#' @param iterlim Iteration limit in nlm
#' @param verbose Whether to print out progress

#' @export 
optimizeLik = function(logLik,phi0,maxIter=NULL,delta=NULL, steptol=NULL, iterlim = NULL, verbose=TRUE) {
  if(is.null(maxIter)) maxIter= 10
  if(is.null(delta)) delta = 1
  if(is.null(steptol)) steptol = 1e-6
  if(is.null(iterlim)) iterlim = 1000

  isOK = FALSE
  iterTry = 0
  bestoptim = list(minimum=Inf)
  while(!isOK) {
    phi0 = phi0 + rnorm(length(phi0),0,delta)
    suppressWarnings({
      tryCatch({
        if(verbose) print(paste0("Trying to optimize function: Iteration ",iterTry + 1,"..."))
        fit = nlm(logLik,phi0,hessian = T,useNegative=T,transform=T,iterlim = iterlim, steptol=steptol) 
        
        #After-check: Ensureing global optimum
        if(fit$min < bestoptim$minimum) {
          bestoptim = fit
          chol(fit$hessian) #check if positive definite
          isOK = TRUE
        } 
        #chol(fit$hessian) #check if positive definite
        #bestoptim = fit
        #isOK = TRUE
        
      },error = function(e) 1)
    })
    iterTry = iterTry + 1
    if(iterTry>=maxIter && bestoptim$minimum < 10^10) {
      isOK = TRUE
      if(verbose) print("Optimization was not successful! Returning from function...")
    }
  }
  bestoptim$iterTry = iterTry
  return(bestoptim)
}