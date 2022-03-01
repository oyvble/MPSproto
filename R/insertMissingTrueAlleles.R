#' @title insertMissingTrueAlleles 
#' @description Indicate missing alleles of true contributor and insert into table if missing
#' @param dat Datatable with columns (RefName,SampleName,Locus,Allele,Coverage,Dose)
#' @export 

insertMissingTrueAlleles = function(dat) {

  #Prepare format (as dataframe):
  dat = data.frame(RefName=dat[,1],SampleName=dat[,2],Locus=dat[,3],Allele=dat[,4],Coverage=as.integer(dat[,5]),Dose=as.integer(dat[,6]))

  refs = unique(dat$RefName)
  samples = unique(dat$SampleName)
  locs = unique(dat$Locus)
  
  for(ref in refs) { #for each unique reference
    samplesRef = unique(dat$SampleName[dat$RefName==ref]) #obtain samples
    for(loc in locs) {
      sub = dat[dat$RefName==ref & dat$Locus==loc,,drop=F]
      trueAlleles = unique(sub$Allele[sub$Dose>0]) #obtain unique true alleles for contributor
      
      #CHECK IF INFO ABOUT TrueAllele (i.e. dose>0) is missing for any of samples. Insert utilizing other info!
      tab = table( factor(sub$SampleName,levels=samplesRef),sub$Allele)
      tab = tab[,colnames(tab)%in%trueAlleles,drop=F] #obtain table with doses
      missInfo = which(tab==0,arr.ind = T)  #obtain index of missing info
      dose = unique(tab[tab>0])
      if(length(dose)>1) stop("Two different doses observed. SampleName possibly not unique!")
      nMisses = nrow(missInfo) #insert missing alleles
      if(nMisses>0) {
        for(missind in 1:nMisses) { #for each missing allele
          insMissing = c(ref, rownames(tab)[missInfo[missind,1]],loc,colnames(tab)[missInfo[missind,2]],0,dose)
          dat = rbind(dat,insMissing) #Insert missing row
        }
      } #  
    }
  }
  dat$Coverage = as.integer(dat$Coverage)
  dat$Dose = as.integer(dat$Dose)
  return(dat)
}
