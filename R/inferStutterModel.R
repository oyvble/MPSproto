#' @title inferStutterModel 
#' @description A function to perform beta-regression analysis 
#' @details The data table can be stutterdata table obtained from the getStutterData function
#' @param stutterdata Datatable with columns ("Locus","Type","BLMM","Prop.)
#' @param print2pdf Name of pdf file for plots
#' @param signifLevel Significance level for keeping BLMM variable
#' @param addValid Whether to add Validation figures (PP-plots) 
#' @param modPrec Whether the model should include precision for BLMM variable
#' @export 

inferStutterModel = function(stutterdata,print2pdf=NULL,signifLevel=0.05,addValid=FALSE,modPrec=FALSE) {
  blmm_sep = "-" #separator for blmm when multiple

  #Check that all column names are included  
  if(!all(c("Locus","Type","Prop","BLMM")%in%colnames(stutterdata) )) stop("Missing column name!")
  
  regFitList = list()  
  locs = unique(stutterdata$Locus)
  for(loc in locs){
# loc="D13S317" 
# loc=locs[3]    
    sub = subset(stutterdata,stutterdata$Locus==loc)
    regFitList[[loc]] = list()
    
    types = unique(sub$Type)
    for(type in types) {
#  type=types[1]
      sub2 = subset(sub,sub$Type==type)
      blmm = sub2$BLMM #covariate
      fitConstant <- FALSE #whether a constant model should be fitted
      if(any(grepl(blmm_sep,blmm))) {
        blmm = strsplit(blmm,blmm_sep)
        sub2$BLMM1 = as.numeric( sapply(blmm,function(x) x[1]) )
        sub2$BLMM2 = as.numeric( sapply(blmm,function(x) x[2]) )
        
        suppressWarnings({
          tryCatch({
            if(!modPrec) {
              fit = betareg::betareg(as.formula("Prop ~ BLMM1 + BLMM2"), data = sub2)
            } else {
              fit <- betareg::betareg(as.formula("Prop ~ BLMM1 + BLMM2|BLMM1 + BLMM2"), data = sub2)
            }
            pvals = coef(summary(fit))$mean[2:3,4]
            if( any( coef(fit)[2:3]<0) || any(pvals>signifLevel)) fitConstant <- TRUE
          } , error = function(e) {
            fitConstant <<- TRUE
            #stop()
          })
        })
        
      } else {
        sub2$BLMM = as.numeric(blmm)
        
        suppressWarnings({
          tryCatch({
            if(!modPrec) {
              fit <- betareg::betareg(as.formula("Prop ~ BLMM"), data = sub2)
            } else {
              fit <- betareg::betareg(as.formula("Prop ~ BLMM|BLMM"), data = sub2)
            }
            pval = coef(summary(fit))$mean[2,4]
            if( coef(fit)[2]<0 || pval>signifLevel) fitConstant <- TRUE
          } , error = function(e) {
            fitConstant <<- TRUE
            #stop()
          })
        })
      }      
      suppressWarnings({
        if(fitConstant) fit = betareg::betareg(as.formula("Prop ~ 1"), data = sub2)
      })
      regFitList[[loc]][[type]] = fit
      #summary(fit)
    }
  }

  #structure fitted models into table  
  betaTable = NULL  
  np2 = 3 #total number of params (b0,b1,b2)
  for(loc in locs){
    types = names(regFitList[[loc]])
    
    for(type in types) {
      fit = regFitList[[loc]][[type]]
      suppressWarnings({
        linpred = coef(summary(fit))$mean
      })
      np = nrow(linpred) #number of param for mean
      row = rep(NA,2*np2)
      row[1:np] <- linpred[,1]
      row[1:np + np2] <- linpred[,4] #insert pvalues
      df = data.frame(Locus=loc,Type=type,Row=t(row))
      betaTable = rbind(betaTable, df)
    }
  }
  colnames(betaTable)[-(1:2)] = c(paste0("beta_",0:2),paste0("pval_",0:2))
  
  #LAST: VISUALIZE DATA WITH FITTED MODEL
  if(!is.null(print2pdf)) {
    pdf(paste0(print2pdf,".pdf"),height=8,width=8)
    for(loc in names(regFitList)){
 # loc=locs[3]    
      for(type in names(regFitList[[loc]]) )  {
#type= names(regFitList[[loc]])[4]
#loc="D12S391"; type="BW1"
        sub = subset(stutterdata,stutterdata$Locus==loc & stutterdata$Type==type)
        blmm = sub$BLMM #covariate
        if(any(grepl(blmm_sep,blmm))) {
          blmm = strsplit(blmm,blmm_sep)
          sub$BLMM1 = as.numeric( sapply(blmm,function(x) x[1]) )
          sub$BLMM2 = as.numeric( sapply(blmm,function(x) x[2]) )
        } else {
          sub$BLMM = as.numeric(sub$BLMM)
        }
        fit = regFitList[[loc]][[type]] #get fitted model
        counts = nrow(sub)   #number of observations

        pred = predict(fit, type="response") #obtain model expectation
        sub = cbind(sub, fit=pred) #insert for each observations
        
        #Always order x-var 
        if( any(grepl(blmm_sep,sub$BLMM)) ) {
          ord = order(sub$fit,decreasing = FALSE)
        } else {
          ord = order(sub$BLMM,decreasing = FALSE)
        }
        sub$BLMM = factor(sub$BLMM,levels=unique(sub$BLMM[ord]))

        main = paste0(loc,": ",type," (",counts,")")
        p <- ggplot2::ggplot(sub, ggplot2::aes(BLMM, Prop))
        p <- p + ggplot2::ggtitle(main)
        p <- p + ggplot2::labs(y="Stutter Prop.", x = "BLMM")        
        p <- p + ggplot2::ggtitle(paste0(loc,": ",type," (",counts,")"))
        p <- p +  ggplot2::geom_point( ggplot2::aes(y = fit), size = 20,color=2,pch="-")
        p <- p + ggplot2::geom_point() #superimpose points last
        #p <- p + geom_ribbon( aes(ymin = lwr, ymax = upr), alpha = .15) 
        #p <- p +  ggplot2::geom_line( ggplot2::aes(y = fit), size = 1,col=2)
        #p <- p + facet_grid(cols = vars( Type)) + 
        print(p)
        
        #LAST PART: POSSIBLE TO INCLUDE VALIDATION PLOT:
        if(addValid) { #include model validation
          ilogit = function(x) 1/(1+exp(-x))
          X = model.matrix(fit) #get model matrix for each observation

          #Obtain shape parameters for each observations:
          betaHat = fit$coefficients$mean #estimated coefficients in linear predictor
          etahati = X%*%betaHat
          muhat = ilogit(etahati) #obtain estimated expectation
          
          phihat = fit$coefficients$prec #constant?
          if(length(phihat)>1) {
            phihat =  c( exp(X%*%fit$coefficients$prec) ) #estimated dispersion
          }
          shapeParamsHat = phihat*cbind(muhat,1-muhat) #shape and scale (maximum likelihood estimate)
          
          #Uniform quantiles (cumulative increasing)
          N = counts #number of observations
          cumunif = ppoints(N)# ((1:n)-0.5)/n
          
          #PP plot: PLUG IN responses into cumulative beta-distr to obtain probabilities
          pp  = pbeta(fit$y,shapeParamsHat[,1],shapeParamsHat[,2])
          plot(cumunif,sort(pp),xlab="Theoretical",ylab="Observed",main="PP plot");abline(a=0,b=1)
          mtext(main)
          
          xsq <- seq(0,1,l=1000) #Draw envolope lines
          alpha = 0.05
          ysq <- c(alpha,alpha/N) #quantiles to consider
          for(qq in ysq) {
            lines(xsq,qbeta(qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
            lines(xsq,qbeta(1-qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
          }
          legend("bottomright",legend=paste0("1-Envelope-coverage=",c("",paste0(alpha,"/",N,"=")),signif(ysq,2)),col=1:length(ysq),lty=2,cex=1.3)
          #QQ plot: PLUG IN cumunif into inverse cumulative beta-distr to obtain quantiles
          #qq  = qbeta(cumunif,shapeParamsHat[,1],shapeParamsHat[,2])
          #plot(sort(qq),sort(fit$y),xlab="Theoretical",ylab="Observed",main="Q-Q plot");abline(a=0,b=1)
        }
      }
    }
    dev.off()
  }

  return( list(regFitList=regFitList,betaTable=betaTable) )
}

