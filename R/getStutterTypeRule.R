#' @title getStutterTypeRule 
#' @description Helpfunction to get motifchange for a given stutter type
#' @details The nomenlclature for stutters are as follows:
#' BW1 means that the LUS motif is providing a backward stutter (n-1).
#' FW1 means that the LUS motif is providing a forward stutter (n+1).
#' BW2 means that the non-LUS motif is providing a backward stutter (n-1).
#' etc... prefix "D" means double (n-2 or n+2)
#' n0 stutters are named as either FWBW or BWFW depending on which occuring first in the sequence (hence orientation matters!)
#' @param type Stutter type. Supports (BW1,BW2,FW1,FW2,DBW1,FWBW,BWFW)
#' @param platform Platform type {"MPS","CE"}. NOTE: CE has only LUS based stutter types
#' @export 

getStutterTypeRule = function(type,platform="MPS") {  
  #Recognize last letter (is it a number?)
  rank = substr(type,nchar(type),nchar(type)) #obtain rank of repeated block (should be LUS=1 or not=2)
  suppressWarnings( {   rank = as.integer(rank)  }   )
  if(is.na(rank)) rank = 3 #set as rank 3 (this is either FWBW or BWFW)

  type0 = type #copy type
  if(rank < 3) type0 = substr(type,1,nchar(type)-1) #extract type  
  if(platform=="CE" && rank==3) type0 = type #CE variant did not contain  number in suffix

  motifchange = NULL
  if(type0=="BW") motifchange = -1 #stutter on longest (or second longest) motif
  if(type0=="FW") motifchange = +1  #stutter on longest (or second longest) motif
  if(type0=="DBW") motifchange = -2 #stutter on longest motif
  if(type0=="DFW") motifchange = +2 #stutter on longest motif
  if(type0=="TBW") motifchange = -3 #stutter on longest motif

  #Handle unique stutter types for each platform
  if(platform=="MPS") {
	if(type0=="FWBW") motifchange = c(+1,-1) #two motifs are stuttering at the same time
	if(type0=="BWFW") motifchange = c(-1,+1) #two motifs are stuttering at the same time
  } else if(platform=="CE") {
	if(type0=="2BP") motifchange = c(-0.8,-0.2) #Can happen in two ways (whole number and not):
	#example: 17 -> 16.2 OR 25.2 -> 25
  }
  
  if(is.null(motifchange)) {
    return(NULL)
  } else {
    return(c(motifchange,rank))
  }
}