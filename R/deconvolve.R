#' @title deconvolve
#' @author Oyvind Bleka
#' @description deconvolve ranks the set of the most conditional posterior probability of genotypes the STR DNA mixture given a fitted model under a hypothesis.
#' @details The procedure calculates the likelihood for each single locus. Then it combines the most probable genotypes from each loci to produce a ranked list of deconvolved profiles.
#' 
#' @param mlefit Fitted object using inferEvidence function.
#' @param alpha Required sum of the listed posterior probabilities.
#' @param maxlist The ranked deconvolved profile list will not exceed this number (used to avoid endless search).
#' @return ret A list(table1,table2,table3,table4,rankGi,rankG,pG) where rankG is the ranked genotypes with corresponding probabilities in pG. rankgGi is the same, but per marker. table1 is rankG and pG combined (joint results). table2 uses rankGi to find marginal results for top-genotypes. table3 and table4 shows this marginalized on genotypes and alleles per contributor 
#' @export

# alpha=0.95;maxlist=1000
deconvolve = function(mlefit,alpha=0.95,maxlist=1000){ 
  par = mlefit$fit$par
  c <- mlefit$prepareC #returned from prepareC
  nC = c$NOC #number of contrs
  nM = c$nLocs #number of markers
  locs = c$locNames #loci to evaluate

  #Step 1) Calculate L(E|g,thetahat) for each marker
  deconvlisti = calcLogLikC_prediction(par,c,TRUE) #returnOutcome (all combinations are there)

  #POST PROCESSING:
  kvec <- 1:nC #index of contributors
  colN <- paste0("C",kvec) #column name of contributors
  
  #Step 2) Convert rank-list to list with genotype/allele-names
  for(loc in locs) {
   tmp = deconvlisti[[loc]] #obtain table
   lastCol = ncol(tmp)
   lik = tmp[,lastCol] #obtain likelihood
   deconvlisti[[loc]][,lastCol] = lik/sum(lik) #and convert it to probability
   ord = order(deconvlisti[[loc]][,lastCol],decreasing = TRUE)
   deconvlisti[[loc]] = deconvlisti[[loc]][ord,,drop=FALSE]
   colnames(deconvlisti[[loc]])[-lastCol] = colN
   colnames(deconvlisti[[loc]])[lastCol] = "Probability"
  }
  
  #Step 3) Create table layouts:
  #Helpfunctions to obtain marginal probabilities
  getMarg <- function(x,y) { #get marginal of genotypes
    agg <- aggregate(y,by=list(x),sum) #get probabilities
    ord <- order(agg[,2],decreasing=TRUE)
    agg2 <- agg[ord,,drop=F]
    colnames(agg2) <- c("Genotype","Probability")
    return(agg2)
  }
  getMarg2 <- function(x,y) { #get marginal of alleles
    tmp <- unlist(strsplit(x,"/"))
    unA <- unique(tmp) #unique alleles
    x2 <- t(matrix(tmp,nrow=2))
    prob <- rep(NA,length(unA))  
    for(aa in unA) prob[which(unA==aa)] <- sum(y[rowSums(x2==aa)>0]) #sum probabilities
    ord <- order(prob,decreasing=TRUE)
    agg <- data.frame(Allele=unA[ord],Probability=prob[ord])
    return(agg)
  }
  maxI <- function(p) min(min(which(cumsum(p)>=alpha)),maxlist,length(p))  #helpfunction to obtain a maximum size of a vector (bounded in both length and probability)
  
  #A) Calculate marginal probabilities for all contributors (genotypes and alleles):
  deconvlistic <- list() #genotype list per contributor
  deconvlistica <- list() #allele list per contributor
  cn <-  c("TopGenotype","probability","ratioToNextGenotype") #names for each contributor
  toplist <- list()
  for(loc in locs) {
    deconvlistica[[loc]] <- deconvlistic[[loc]] <- list()
    X <- deconvlisti[[loc]]
    nc <- ncol(X) #number of column
    nr <- nc - 1 #number of contributors
    tab <- matrix(,nrow=3,ncol=nr)
    rownames(tab) <- cn 
    colnames(tab) <- colN
    for(rr in 1:nr) {
      deconvlistic[[loc]][[colN[rr]]] <- tmp <- getMarg(x=X[,rr],y=as.numeric(X[,nc]))
      deconvlistica[[loc]][[colN[rr]]] <- getMarg2(x=X[,rr],y=as.numeric(X[,nc]))
      rat <- ifelse(nrow(tmp)>1, tmp[1,2]/tmp[2,2],NA) #get ratio from first to second genotype
      tab[,rr] <- c(tmp[1,1],signif(tmp[1,2],4),signif(rat,4)) 
    }
    toplist[[loc]] <- tab
  }
  
  #B) Create tables
  table1 <- table2 <- table3 <- table4 <- numeric()
  for(loc in locs) {
    combs <- deconvlisti[[loc]]
    prob <- as.numeric(combs[,ncol(combs)])
    combs[,ncol(combs)] <- signif(prob,4)
    maxind <- maxI(prob)
    combs <- combs[1:maxind,,drop=F]
    table1 <- rbind(table1,cbind(loc,1:nrow(combs),combs))
    table1 <- rbind(table1, rep("",ncol(table1)) )
  }
  colnames(table1)[1:2] <- c("Locus","Rank")
  
  for(loc in locs) table2 <- rbind(table2, c(toplist[[loc]]) )
  colnames(table2) <- paste0(cn,"_",c(t(replicate(length(cn),colN))))
  rownames(table2) <- locs
  
  maxI2 <- function(p) min(max(which(p>(1-alpha))),maxlist,length(p))
  space <- cbind("","","","")
  for(cc in colN) {
    for(loc in locs) {
     tmp <- deconvlistic[[loc]][[cc]]
     maxind <- maxI(p=tmp$Probability)
     newrow <- tmp[1:maxind,,drop=F]
     newrow[,2] <- signif(newrow[,2],4)
     newrows <- as.matrix(cbind(cc,loc,newrow))
     table3 <- rbind(table3,newrows,space)
    
     tmp <- deconvlistica[[loc]][[cc]]
     maxind <- maxI2(p=tmp$Probability)
     newrow <- tmp[1:maxind,,drop=F]
     newrow[,2] <- signif(newrow[,2],4)
     newrows <- as.matrix(cbind(cc,loc,newrow))
     table4 <- rbind(table4,newrows,space)
    }
  }
  colnames(table3)[1:2] <- colnames(table4)[1:2] <- c("Contr.","Locus")
  
  return(list(table1=table1,table2=table2,table3=table3,table4=table4,toprankGi=toplist,rankGi=deconvlisti))
} #end function

