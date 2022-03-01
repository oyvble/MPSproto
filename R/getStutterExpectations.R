
#' @title getStutterExpectations
#' @author Oyvind Bleka
#' @description Obtain expected stutter proportions
#' @details Gives an overview over stutter proportion between each of the sequences
#' 
#' @param c Output from prepareC_prediction
#' @return A list
#' @export

getStutterExpectations = function(c) {
  stuttList = list()
  for(m in seq_len(c$nLocs) ) {
    loc = c$locNames[m]
    av = c$alleleNames[[loc]] #allele vector
    
    nStutt = c$nStutters[m] #number of stutters
    stuttTab = NULL
    if(nStutt>0) {
      startInd = c$startIndMarker_nStutters[m] + seq_len(nStutt)
      stuttType = c$stuttType[startInd ]
      stuttFrom = c$stuttFromInd[startInd ] + 1
      stuttTo = c$stuttToInd[startInd ] + 1 
      stuttExp = c$stuttExp[startInd]
      
      stuttTab = cbind(stuttType,av[stuttFrom],av[stuttTo],stuttExp)
      colnames(stuttTab) = c("Type","From","To","Stutterprop.")
    }
    stuttList[[loc]] = stuttTab
  }
  return(stuttList)
}
