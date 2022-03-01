#' @title getStutteredSequence 
#' @description Helpfunction for obtaining different types of stutters 
#' @param seqs Input sequences
#' @param type Stutter type. Supports (BW1,BW2,FW1,FW2,DBW1,TBW1,DFW1,FWBW,BWFW)
#' @param platform Platform type {"MPS","CE"}. NOTE: CE has only consider only some stutter types (BW1,DBW1,FW1)
#' @export 

getStutteredSequence = function(seqs,type="BW1",platform="MPS") {  
  
  #Obtain motif change and rank
  stutterRule = getStutterTypeRule(type,platform)
  if(is.null(stutterRule)) stop("Stutter type not implemented")
  nmotifs = length(stutterRule) - 1 #number of motif changes
  motifchange = stutterRule[1:nmotifs]
  rank = stutterRule[nmotifs+1] #last element is rank
  
  nSeqs = length(seqs)
  motifRepList = getMotifRepsSequences(seqs) #obtain motif counts per sequence
  BLMM = rep(NA,nSeqs) #Oobtain block length of stutter
  for(i in seq_along(motifRepList)){
    #i=1
    ord = order(motifRepList[[i]],decreasing=T) #obtain order
    if(rank%in%(1:2)) { #if rank is defined
      if(rank>length(ord)) { #number of motifs was less than rank
        motifRepList[[i]] = NA 
        next #skip if not possible to obtain corresponding rank index
      } else { #this is default (rank <= number of motifs)
        motifMutateInd = ord[rank] #index to mutate
        BLMM[i] <-  motifRepList[[i]][motifMutateInd] #obtain block length
        motifRepList[[i]][motifMutateInd] = BLMM[i] + motifchange  #update
        
      }
      #If rank not defined (this means complex motifchange is defined)
    } else {
      if(length(motifchange)>length(ord)) {  #number of motifs was less than length of motif change
        motifRepList[[i]] = NA 
        next 
      } else {
        blmm = NULL
        for(j in 1:length(motifchange)) { #for each motif change
          motifMutateInd = ord[j] #index to mutate
          blmm1 <- motifRepList[[i]][motifMutateInd]  #obtain block length
          motifRepList[[i]][motifMutateInd] =  blmm1 + motifchange[j] #update
          blmm = c(blmm,blmm1)
        }
        BLMM[i] = paste0(blmm,collapse="/") #if seveal
      }
    }
  }    
  stutters = convertMotifReps2seq(motifRepList) #obtain stuttered alelels
  return(list(stutter=stutters,BLMM=as.character(BLMM)))
}