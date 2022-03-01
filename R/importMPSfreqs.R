#' @title importMPSfreqs
#' @description importMPSfreqs is a function to import allele frequencies directly from file(s) supporting text format
#' @details
#' The function reads local text file as input, and return a list of populations.
#'
#' @param fn File name(s).
#' @return A list giving the population frequencies for every population
#' @export

importMPSfreqs = function(fn) {
 if(is.null(fn) || length(fn)==0) {
  stop("No frequency file were given!")
 }
 dbList <- list() #list of populations
 for(ff in fn) { #for each file
   tab=tableReader(ff)
   Anames = tab[,1] #first column is allele frequeneies
   tab = tab[,-1,drop=FALSE] 
   freqlist = list()
   for(j in 1:ncol(tab)) { #for each locus
     tmp = tab[,j]
     tmp2 = tmp[!is.na(tmp) & as.numeric(tmp)>0] #require that allele is not NA and is>0
     names(tmp2) = Anames[!is.na(tmp)]
     freqlist[[j]] = tmp2
   }
   names(freqlist) = toupper(colnames(tab)) #LOCUS-names are assigned as Upper-case! This is important to do!
   pop = unlist(strsplit(basename(ff),"\\."))[1] #get population name (remove '.' symbols)
   dbList[[pop]] <- freqlist
 } #end for each file ff in fn
 return(dbList)
}
