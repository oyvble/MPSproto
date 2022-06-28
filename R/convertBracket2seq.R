#' @title convertBracket2seq 
#' @description Helpfunction to convert bracket format into full length sequence
#' @param seqs Input sequence in bracket format
#' @param brksign Bracket format break symbols
#' @export 

convertBracket2seq = function(seqs,brksign=c("\\[","\\]")) {
  # All sequences are split with pattern " ". Output is a list of lists.Every sequence
  # is a list of tetranucleotides!!!
  xList = sapply(seqs,strsplit, split =" ")
  
  # Secondary function that will extract the motifs, repeat them and collapse:
  motivosfun = function(xList) { # Function to go through every row of list xVList
    seq = "" #init full sequence
    for(i in 1:length(xList)) { # Iteration through all the lists (rows) of xList
      tmp = xList[i]
      if( grepl(brksign[1],tmp) ) {
        ret = getMotifReps(tmp,brksign)  # Obtain motif and number of motifs
        motif = ret[1]
        motifReps = ret[2]
        tmp = paste(rep(motif,motifReps),collapse="")
      }
      seq = paste0(seq,tmp) ##
    }
    return(seq)
  }
  
  lapply(xList, motivosfun) ## It passes the function motivosfun to every element of the lists within the list xList
}

