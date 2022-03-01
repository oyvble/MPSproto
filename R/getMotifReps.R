#' @title getMotifReps 
#' @description Helpfunction to Get motifs and repetitions from short bracket format
#' @param block Input block of sequence 
#' @param brksign Bracket format break symbols
#' @export 

getMotifReps = function( block, brksign=c("\\[","\\]") ) {
  if( grepl(paste0(brksign[1]),block) ) {
    tmp2 = strsplit(block,brksign[2])[[1]]
    motif = gsub(brksign[1],"",tmp2[1]) #motif string (without brackets)
    motifReps = tmp2[2] #convert number of repititions to number
  } else {
    motif = block
    motifReps = 1 #only 1 motif repetition
  }
  return(c(motif,motifReps)) #return motif string and number of motif
}