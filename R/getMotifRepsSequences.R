#' @title getMotifRepsSequences 
#' @description Helpfunction for obtaining number of motif repeats for each sequence and blocks in each sequence
#' @param seqs Input sequences in string format
#' @export 

getMotifRepsSequences = function(seqs) {
  splittedSeqs = strsplit(seqs," ") #split up sequences into 
  MotifRepsSeqs = list()
  for(i in seq_along(splittedSeqs)){
    nMotifReps = sapply(splittedSeqs[[i]],getMotifReps)
    nreps = as.integer(nMotifReps[2,])
    names(nreps) = nMotifReps[1,]
    MotifRepsSeqs[[i]] =   nreps
  }
  return(MotifRepsSeqs)
}