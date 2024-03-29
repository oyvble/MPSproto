#' @title importMPSsample
#' @description Converting table to list format (helpfunction)
#' @param fn A file name to import data from (Evid or refs)
#' @return outL A data list
#' @export

importMPSsample = function(fn) {
  table = tableReader(fn)
  colname = colnames(table) #colnames in file
  lind = grep("marker",tolower(colname),fixed=TRUE) #locus col-ind
  if(length(lind)==0) lind = grep("loc",tolower(colname),fixed=TRUE) #try another name
  sind = grep("sample",tolower(colname),fixed=TRUE) #sample col-ind
  if(length(sind)>1)  sind = sind[grep("name",tolower(colname[sind]),fixed=TRUE)] #use only sample name
  
  A_ind = grep("allel",tolower(colname),fixed=TRUE) #allele col-ind
  if(length(A_ind)==0) A_ind  = grep("seq",tolower(colname),fixed=TRUE) #try another name (sequence)
  H_ind = grep("heig",tolower(colname),fixed=TRUE) #height col-ind
  if(length(H_ind)==0) H_ind  = grep("cov",tolower(colname),fixed=TRUE) #try another name (coverage)
  if(length(H_ind)==0) H_ind  = grep("int",tolower(colname),fixed=TRUE) #try another name  (intensity)
  if(length(H_ind)==0) H_ind  = grep("read",tolower(colname),fixed=TRUE) #try another name  (read)
  
  locs = unique(toupper(table[,lind])) #locus names: Use uniques and Convert to upper case
  samplenames = unique(as.character(table[,sind])) #sample names
  outL = list() #Init outList (insert non-empty characters):
  for(samplename in samplenames) { #for each sample in matrix
    outL[[samplename]] = list() #one list for each sample
    for(loc in locs) { #for each locus
      rowinds = which(table[,sind]==samplename & toupper(table[,lind])==loc) #get row index in table for given sample and locus
      if(sum(rowinds)==0) next #no data found (skip marker)
      
      alleles <- heights <- numeric() #init vector with data
      for(rowind in rowinds) { #for each row indices
        keep <- !is.na(table[rowind,A_ind]) & table[rowind,A_ind]!="" #get boolean vector of alleles to keep
        alleles <- c(alleles, table[rowind,A_ind[keep]] ) #extract allele(s), be be levels
        if(length(H_ind)>0) {
          heights <- c(heights, table[rowind,H_ind[keep]]) #extract peak heights (intensities)
        }
      } #end for each rows
      
      insAlleles <- character() #default if no data
      insHeights <- numeric() #default  if no data
      if(length(alleles)>0) { #if any data to store:
        if(length(A_ind)>0) insAlleles = as.character(alleles)
        if(length(H_ind)>0) insHeights = as.numeric(as.character(heights))
      }
      if(length(A_ind)>0) outL[[samplename]][[loc]]$adata = insAlleles
      if(length(H_ind)>0) outL[[samplename]][[loc]]$hdata = insHeights
    } #end for each loci
  } #end for each samples
  return(outL)
}