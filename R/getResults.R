#Small helpfunctions to get results after fitting MLE (hp and Hd)

getLR = function(mle1,mle2) exp(mle1$fit$loglik- mle2$fit$loglik)

getlog10LR = function(mle1,mle2) (mle1$fit$loglik- mle2$fit$loglik)/log(10)

getPar = function(mle1,mle2) cbind(Hp=unlist(mle1$fit$par),Hd=unlist(mle2$fit$par))

getLRi = function(mle1,mle2) exp(logLiki(mle1)- logLiki(mle2))

getAIC = function(mle1,mle2=NULL) {
 AIC1 = -2*mle1$fit$loglik + 2*length(mle1$fit$phihat) 
 if(!is.null(mle2)) {
	AIC2 = -2*mle2$fit$loglik + 2*length(mle2$fit$phihat)
	return(AIC2-AIC1)	
 } else {
	return(AIC1)
 } 
}
