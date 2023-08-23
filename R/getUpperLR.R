#' @title getUpperLR
#' @description Obtaining upper boundary LR for POI 
#' @details The function takes a fitted MLE object under Hd (where hypothesis$knownRef is the POI)
#' 
#' @param mle Fitted object using inferEvidence function.
#' @param scale Whether to scale the RMP when having more than 2 unknowns (recommended)
#' @return Returning upper boundary LR (log10 scale)
#' @export

getUpperLR = function(mle, scale=TRUE) {
  fst0 = mle$hypothesis$fst
  cond = which(mle$hypothesis$cond>0) #obtain conditionals
  POI = mle$hypothesis$knownRef
  if(length(POI)!=1) return(NULL) #print("Couldn't calculate upper boundary LR. Returning...")
    
  locs = names(mle$data)
  rmp = setNames( rep(1,length(locs)) , locs)
  for(loc in locs) {
    datLoc = mle$data[[loc]]
    refA = datLoc$refs[[POI]] #get alleles of ref
    freqRef = datLoc$freq[ match(refA,names(datLoc$freq))] #obtain freq of ref-alleles
    
    #Obtain number of typed alleles
    nTyped = table(c(unlist(datLoc$refs[cond]),refA)) #obtain typed alleles
    totTyped = sum(nTyped) #total typed alleles
    
    gprob = 1 #prod(freqRef) #genotype prob
    for(a in 1:2) {
      aind = which(names(nTyped)==refA[a]) #get correct index
      gprob = gprob * (fst0*nTyped[aind] + (1-fst0)*freqRef[a]) / (1 + (totTyped-1)*fst0); 
      nTyped[aind] = nTyped[aind] + 1; #update allele count 
      totTyped = totTyped + 1; #update total count
    }
    if(refA[1]!=refA[2]) gprob = 2*gprob #get het situation
    rmp[loc] = gprob  
  }
  
  #THIS LAST BLOCK IS TO TAKE INTO ACCOUNT THE POSSIBILITY THAT 
  #THE LR CAN EXCEED 1/RMP WHEN HAVING MORE THAN 2 unknowns
  scale0 = 1 #default is no scaling
  nCond = length(cond) #number of conditionals
  nU = mle$hypothesis$NOC - nCond  #get number of unknowns (Hd)
  if(scale) { #need to scale rmp
    scale0 = (1+(3+2*nCond)*fst0)/(1+(1+2*nCond)*fst0)*(1+(4+2*nCond)*fst0)/(1+(2+2*nCond)*fst0) 
  }
  upperLR = -sum(log10(rmp/scale0))  #get upper limit of LR (log10 scale), fst taken into account
  return(upperLR)
}