#' @title convertMotifReps2seq 
#' @description Helpfunction to convert sequence in motifReps format back to sequence
#' @param motifRepList List with motif blocks
#' @param brksign Bracket format break symbols
#' @export 

convertMotifReps2seq = function(motifRepList,brksign=c("[","]")) {
  seqs = rep(NA,length(motifRepList))
  for(i in 1:length(motifRepList)){
    if(is.na(motifRepList[[i]][1])) next
    motifcount = motifRepList[[i]] #obtain motif counts (each motif)
    inclBrk = motifcount>1 #indices including bracket
    rmMotif = motifcount<1 #motif should be removed
    motifs = names(motifcount)
    motifs[inclBrk] = paste0(brksign[1],motifs[inclBrk],brksign[2],motifcount[inclBrk])
    seqs[i] = paste0(motifs[!rmMotif],collapse=" ")
  }  
  return(seqs)
}