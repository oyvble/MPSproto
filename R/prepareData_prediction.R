#' @title prepareData_prediction
#' @author Oyvind Bleka 
#' @description Reorganasing data for calculation preparations per marker (simpler structure)
#' @details Helpfunction to reorganising data for prepareC input 
#' Also takes detection thresholds as argument to remove possible alleles below the thresholds.
#' 
#' @param samples A list with evidence profiles [[sample]][[locus]]$adata/hdata
#' @param refData A list with reference profiles [[reference]][[locus]]...
#' @param popFreq A list with population frequencies [[locus]]
#' @param minF The freq value included for new alleles (new alleles as potential stutters will have 0). Default NULL is using min.observed in popFreq.
#' @param normalize A boolean of whether normalization should be applied or not. Default is FALSE.
#' @param AT A detection threshold value or thresholds per marker. NULL ignores filtering.
#' @param fillHomGen A boolean of whether to fill in homozygote genotypes given with one allele
#' @param verbose Whether to print out if new alleles are imputed
#' @param doQ Whether to do Q-assignation of allele outcome (Recommended!)
#' @return Restructured data in list used as input for functions evaluating the liklihood function
#' @export

#minF=NULL;normalize=FALSE;fillHomGen=TRUE;AT=NULL;verbose=FALSE;doQ=TRUE
prepareData_prediction = function(samples,refData=NULL,popFreq=NULL,minF=NULL,normalize=FALSE,AT=NULL,fillHomGen=TRUE,verbose=FALSE,doQ=TRUE) { #Helpfunction to get data to analyse
  if(is.null(popFreq)) stop("Populatation frequency object must be provided")
  Qallele="99" #Name of allele given if missing in evidence. Default is 99. This is important when considering the degradation model since 99 is closest to maximum allelein a locus. 
  MPSsymbol = ":" #Used for extracting CE for MPS strings. Example is "10:[ATCG]10".
  #NOTE THAT samples can contain the CE (in front). which is necessary for degradation model (should not be included for others)
  
  locs <- names(popFreq) #get loci in popFreq (decides order)
  sampleNames = names(samples) #obtain sample names 
  refNames = names(refData) #obtain reference names
  
  if(is.null(minF)) minF <- min(unlist(popFreq)) #lowest observed frequency if given as NULL
  
  #Update reference data list element order 
  refData2 = NULL
  if(!is.null(refData)) {
    refData2 <- list()
    for(loc in locs)  {
      refData2[[loc]] <- lapply(refData,function(x) { #THIS ENSURES THAT ALL REFERENCE NAMES ARE INCLUDED!
        av = unlist(x[[loc]])#$adata #get alleles
        if(fillHomGen && length(av)==1) av = rep(av,2) #alleles given only ones m
        return(av) #return selected loci
      })
      if(!all(names(refData2[[loc]])%in%refNames)) stop("Some of the markers did not contain reference info!")
    }
  }
  
  #Traverse each locus and insert re-structured data
  datList = list()
  for(loc in locs)  {
#    loc=locs[27]    
    samples2 = list()
    for(sample in sampleNames) { #for each sample
      
      AT0 = 0 #defualt is no threshold
      if(!is.null(AT))  {
        if(length(AT)==1) {
          AT0 = AT #is constant
        } else if(is.vector(AT)) {
          AT0 = AT[loc] #assume a constant per replicate
          if(is.na(AT0)) stop(paste0("AT was not found for marker ",loc))
        } 
      }

      #Extract evidence data above defined threshold:
      samples2[[sample]] = setNames(numeric(),character()) #set as empty by default (IMPORTANT FOR AVOIDING EMPTY MARKER)
      if(is.null(samples[[sample]][[loc]])) next #Skip if no data
      
      av = samples[[sample]][[loc]]$adata #get alleles
      hv = samples[[sample]][[loc]]$hdata #get heights
      if(length(av)!=length(hv)) stop(paste0("At locus ",loc,": Number of alleles in sample was not same as number of heights!"))
      if(length(av)>0) { #if contains alleles empty
        keep = hv>=AT0
        av = av[keep]
        hv = hv[keep]
      } 
      samples2[[sample]] = setNames(hv,av)
    } #end for each sample
    
    alleles = NULL #empty by default
    if(length(samples2)>0) { #skip marker if no evid data were present (THIS IS IMPORTANT!)
      alleles = unique(unlist(lapply(samples2,names))) #obtain observed alleles
      freqs = popFreq[[loc]] #obtain freqs
      
      #INCLUDE THE POSSIBILITY THAT ALLELES CAN HAVE FORMAT "CE:STRING" (KEEP ONLY STRING)
      indChange = grep(MPSsymbol,alleles)
      for(ind in indChange) alleles[ind] = strsplit(alleles[ind],MPSsymbol)[[1]][2] 
  
      newa <- alleles[!alleles%in%names(freqs)]   #get alleles not in popFreq-table
      if(length(newa)>0) {
        tmp <- names(freqs)
        freqs <- c(freqs,rep(as.numeric(minF),length(newa)))
        names(freqs) <-  c(tmp,newa)
        if(verbose) print(paste0("Locus ",loc,": Allele(s) ",paste0(newa,collapse=",")," was inserted with frequency ",minF))
        
        if(as.logical(normalize)) { #Normalization is an option
          freqs <- freqs/sum(freqs) #normalize
          if(verbose) {
            print(paste0("New frequencies for locus ",loc))
            print(freqs) 
          }
        }
      }
    }
    
    #Prepare frequencies:
    freqs2 <- freqs #init
    
    #Perform Q-assigniation for freqs and references (depending on sample observations)
    if(doQ) {
      tmp <- freqs[names(freqs)%in%alleles] #find observed alleles
      if( length(tmp) < length(freqs)) {  
        tmp <- c(tmp,1-sum(tmp))
        names(tmp)[length(tmp)] <- Qallele
      }
      freqs2 = tmp
  
      if(!is.null(refData2)) { #insert 99 as default allele of missing refs
        if(!all(unlist(refData2[[loc]])%in%alleles)) { #there was some missing alleles
            for(k in 1:length(refData2[[loc]])) refData2[[loc]][[k]][! refData2[[loc]][[k]]%in%alleles ] <- Qallele #insert missing     
        }
      }
    }
    
    #insert to retur data list
    datList[[loc]] = list(freq=freqs2,samples=samples2,refs=refData2[[loc]])
    
  } #end for each loci
  return(datList)
}