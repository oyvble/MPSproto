#' @title getStutterIndex 
#' @description Helpfunction for obtaining stutter index between alleles for different types of stutters 
#' @details The nomenlclature for stutters are as follows:
#' BW1 means that the LUS motif is providing a backward stutter (n-1).
#' FW1 means that the LUS motif is providing a forward stutter (n+1).
#' BW2 means that the non-LUS motif is providing a backward stutter (n-1).
#' etc... prefix "D" means double (n-2 or n+2)
#' n0 stutters are named as either FWBW or BWFW depending on which occuring first in the sequence (hence orientation matters!)

#' @param seqs Input sequences
#' @param type Stutter type. Supports (BW1,BW2,FW1,FW2,DBW1,FWBW,BWFW)
#' @param platform Platform type {"MPS","CE"}. NOTE: CE has only LUS based stutter types
#' @param locus Locus names (necessary for some markers of CE to obtain useful BLMM)
#' @examples
#'\dontrun{
#' getStutterIndex(c("[ATCG]10","[ATCG]9"),"BW1")
#' getStutterIndex(c("[ATCG]10","[ATCG]9"),"BW2")
#' getStutterIndex(c("[ATCG]9","[ATCG]5 TTGA [ATCG]10","[ATCG]5 TTGA [ATCG]9"),"BW1")
#' getStutterIndex(c("[ATCG]9","[ATCG]5 TTGA [ATCG]10","[ATCG]5 TTGA [ATCG]9"),"BW2")
#' getStutterIndex(c("[ATCG]9","[ATCG]5 TTGA [ATCG]10","[ATCG]4 TTGA [ATCG]10"),"BW2")
#'}
#' @export 


getStutterIndex = function(seqs,type="BW1",platform="MPS", locus=NULL) {
  sep ="/" #separator if multiple stutter candidates
  sep2 ="-" #separator if multiple motifs of candidates
  seqs = as.character(seqs) #be sure that seqs is a string
  nSequences = length(seqs) #number of sequences
  if(nSequences==1) return(NULL) #return if none to compare with

  #Obtain motif change and rank
  stutterRule = getStutterTypeRule(type,platform)
  if(is.null(stutterRule)) stop("Stutter type not implemented")
  nmotifs = length(stutterRule) - 1 #number of motif changes
  motifchange = stutterRule[1:nmotifs]
  rank = stutterRule[nmotifs+1] #last element is rank
  
  #The idea is to store stutter index between all sequencies (including BLMM)
  SIvec = rep(0,nSequences) #[i] <- j #Give index of parental allele (0=none)
  BLMM <- SIvec #store corresponding parental BLMM 
  #For each sequence obtain number of repeats for each motif:
  if(platform=="CE") {
    suppressWarnings( {   seqs2 = as.numeric(seqs)  }   )
    if(any(is.na(seqs2))) stop("Non-numeric allele found!")
    
    stuttTo1 = as.character(seqs2 + motifchange[1]) 
    ind1 = match( seqs, stuttTo1)
    if( length(motifchange)==2) { #In case of "2BP" (special handling)
      isWhole = round(seqs2)==seqs2 #indicate whole numbers
      ind1[!isWhole[ind1]] = NA #Remove matches which not originate from whole numbers
      stuttTo2 = as.character(seqs2 + motifchange[2]) #Check non-whole numbers
      ind2 = match( seqs, stuttTo2)
      ind2[isWhole[ind2]] = NA #Remove matches which not originate from non-whole numbers
      
      notNA = !is.na(ind2)
      ind1[notNA] = ind2[notNA] #insert 
    }
    SIvec = ind1
    SIvec[is.na(SIvec)] = 0
    
    insBLMM = SIvec>0
    BLMM[insBLMM] = seqs2[SIvec[insBLMM]] #insert parental allele as BLMM
    
    #Special handling of BLMM for some markers:
    if(!is.null(locus)) {
      dec = (BLMM%%1) #obtain decimal part only
      isDec = dec>0 #indiciate alleles with decimal points
      
      if(locus=="D1S1656")  BLMM[isDec] = BLMM[isDec]-(4+dec[isDec])
      if(locus=="SE33") BLMM[isDec] = BLMM[isDec]-(7+dec[isDec])
      if(locus=="D21S11") BLMM[isDec] = BLMM[isDec]-(2+dec[isDec])
      if(locus=="D2S441") { #other rule apply:
        isDec = BLMM>=14 #shift if allele is at least 14 (stutter allele is 13)
        BLMM[isDec] = BLMM[isDec]-3.5
      } 
    }
    
    names(SIvec) <- names(BLMM) <- seqs
    return(list(SI=SIvec,BLMM=BLMM))
  }
  
  #FURTHER IS ONLY FOR MPS (bracket format)
  motifRepList = getMotifRepsSequences(seqs) 
  seqInd = seq_len(nSequences)

  #obtain number of motifs for each sequenes
  nmotifRepList = sapply(motifRepList,length) 
  if( rank>1 && sum(nmotifRepList>1)<=1 ) return(NULL) #return empty if sequence doesn't have more motifs (require at least 2)
    
  for(i in seqInd){ #traverse each sequence (this is row)
#    i=1
    stuttProd = motifRepList[[i]] #potential stutter product to evaluate
    stuttFrom_ind = setdiff(seqInd,i) #sequences to compare with (this is column)
    
    for(j in stuttFrom_ind) {
#j=2 
      stuttFrom = motifRepList[[j]]
      if(length(stuttProd)!=length(stuttFrom)) next #skip if different length

      if(!all(names(stuttProd)==names(stuttFrom)) ) next #skip if different length
      
      blockdiff =  stuttProd - stuttFrom #obtain difference vector (per block)
      blockdiffInds = which(blockdiff!=0)
      nblockDiffs = length(blockdiffInds) #number of block differences
      
      #Require to have either 1 (if ranke=1,2) or 2 (if rank=3) motif differences
      if(nblockDiffs==0) next
      if(rank%in%(1:2) && nblockDiffs>1) next
      if(rank==3 && nblockDiffs!=2) next #require exactly 2 differences
      
      #CHECK THAT MOTIF CHANGE IS CORRECT:
      if(!all(blockdiff[blockdiffInds]==motifchange)) next
      
      #Obtain what is LUS
      LUSval = max(stuttFrom)  #obtain LUS value (can be on multiple positions)
      BLMMs = stuttFrom[blockdiffInds] #obtain missing block motifs
      
      isMatch = TRUE #ok by default
      if(rank==1) { #This looks up LUS
        if(BLMMs!=LUSval) isMatch = FALSE #LUS was not BLMM
        
      } else if(rank==2) {
        LUSind = stuttFrom==LUSval #obtain LUS indices 
        
        #Not a mach if the motif was the only LUS (ok if multiple!)
        if(sum(LUSind)==1 && BLMMs==LUSval ) isMatch = FALSE #LUS was BLMM
      }
      
      if(isMatch) {
        BLMMins = BLMMs #Obtain BLMM of matching Motif(s)
        if( length(BLMMins)>1 ) BLMMins = paste0(BLMMins,collapse=sep2) #keep both motifs
        
        if(SIvec[i]==0) { #nothing insertet already
          SIvec[i] = j
          BLMM[i] = BLMMins #insert number of repeats
        } else {
          #warning("Multiple candidates found!")
          SIvec[i] = paste0(SIvec[i],sep,j)
          BLMM[i] = paste0(BLMM[i],sep,BLMMins)
        }
      }
    }
  }
  names(SIvec) <- names(BLMM) <- seqs

  return(list(SI=SIvec,BLMM=BLMM))
}