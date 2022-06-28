#' @title sample_listToTable
#' @author Oyvind Bleka
#' @description Converting profiles from list to table format (helpfunction)
#' @param outL A list on form outL[[samplename]][[locusname]]$adata,outL[[samplename]][[locusname]]$hdata (Evid or refs)
#' @return table 
#' @export

sample_listToTable = function(outL) {
  sn = names(outL) #obtain sample names
  aM = 0   #count number of max allele data:
  hM = 0   #count number of max allele heights:
  checkElem = outL[[1]][[1]] #check first element
  isEvid = "adata"%in%names(checkElem) #Require adata
  if(isEvid && !"hdata"%in%names(checkElem)) stop("Missing hdata name in list structure!")
  hM <- 0 #default is no height data
  for(ss in sn) { #for each sample
    if(isEvid) {
      aM = max(unlist( lapply(outL[[ss]],function(x) length(x$adata)) ),aM)
      hM = max(unlist( lapply(outL[[ss]],function(x) length(x$hdata)) ),hM)
    } else {
      aM = max(unlist( lapply(outL[[ss]],function(x) length(unlist(x)) ),aM))
    }
  }
  #create tables:
  table=numeric()
  for(ss in sn) { #for each sample
    newsample=numeric() #for allele
    ln = names(outL[[ss]])
    for(loc in ln) {
      if(isEvid) {
        newrow = outL[[ss]][[loc]]$adata
      } else {
        newrow = unlist(outL[[ss]][[loc]])
      }
      newsample = rbind(newsample, c(newrow,rep("",aM-length(newrow))))
    }
    newsample2=numeric() #for heights
    if(hM>0) {
      for(loc in ln) {
        newrow = outL[[ss]][[loc]]$hdata
        newsample2 = rbind(newsample2, c(newrow,rep("",hM-length(newrow))))
      }      
    }
    table = rbind(table,cbind(ss,ln,newsample,newsample2))
  }
  cn = c("SampleName","Marker", paste("Allele",1:aM,sep=""))
  if(hM>0) cn = c(cn,paste("Height",1:hM,sep=""))
  colnames(table)  = cn
  return(table)
} #end of functions

