#' @title getFilteredData 
#' @description Obtain filtered version of data where parental alleleles are not too close
#' @details The purpose of this filter is to avoid that two parental alleles interfers with each others in the stutter analysis
#' We check the distance between the two alleles:
#' The marker is not used if one is in the back stutter of the other (dist=1)
#' The CE is tailored to utilize as much data as possible
#
#' @param dat Datatable with columns ("SampleName","Locus","Dose")
#' @param platform MPS or CE
#' @export 


#dat = readRDS("C:\\Users\\oyvbl\\Dropbox\\Forensic\\MPSproject\\MPSproto\\ForenSeq_RebeccaData_Sensitivitydata\\ForenSeq_RebeccaSensitivityData.RDS")

getFilteredData = function(dat, platform="MPS") {
	#Filter situations where two alleles are too close (this is tailored!)
	keepInd = rep(TRUE,nrow(dat))
	samples = unique(dat$SampleName)
	locs = unique(dat$Locus)
	for(sample in samples) {
#sample=samples[2]
	  for(loc in locs) {
#loc=locs[1]
		indrows = dat$SampleName==sample & dat$Locus==loc
		sub = subset(dat,indrows) 
		
		donor = sub$Allele[sub$Dose>0]
		if(length(donor)<2 || length(donor)==nrow(sub)) next 
		
		if(length(donor)>2) stop("Not possible with more than 2 true alleles")

		if(platform=="CE") {
	
		  dist = abs(diff(as.numeric(donor))) #use absolute value
		  if(dist==1) { #remove all if within 1 allele
			keepInd[indrows] = FALSE #dont keep 
			
			#Special handling if 2-3 in difference        
		  } else if(dist%in%c(2,3)) {
			firstAllele = min(as.numeric(donor))
			av = as.numeric(sub$Allele)
			
			if(dist==2) rmAlleles = av>firstAllele #remove alleles after first allele
			if(dist==3) rmAlleles = av==(firstAllele+1) #remove allele which is 'first allele'+1

			#keep first allele only (but remove all alleles after)
			keepInd[which(indrows)[rmAlleles]] = FALSE #dont keep alleles after first allele
		  } else if(dist<2) { #special handling (other possible stutters)
			dec = round( dist%%1*10 )
			if( dec%in%c(2,8) ) { #check if in 2bp position of eachother
			  firstAllele = min(as.numeric(donor))
			  secondAllele = max(as.numeric(donor))
			  firstIsDec = firstAllele%%1>0
			  
			  alleleRM = NULL #default is not alleles removed
			  if(dist<1) {
				  alleleRM = secondAllele-1 #this is allele to remove if found (masked with 2bp)
			  } else { #if distance is more than 1
				if(firstIsDec) {
				  alleleRM = c(firstAllele + 1,secondAllele-2) #these are allele to remove 
				}
			  }
			  rmAlleles = sub$Allele==alleleRM #remove alleles after first allele
			  keepInd[which(indrows)[rmAlleles]] = FALSE #dont keep alleles after first allele
			} 
		 }
		}else if(platform=="MPS") {
		  indBW1 = getStutterIndex(donor,"BW1")$SI #get stutter index of type BW1
		  indBW2 = getStutterIndex(donor,"BW2")$SI #get stutter index of type BW1
		  if(!all(indBW1==0)) { #if one is in stutter position of the other
		    keepInd[indrows] = FALSE 
		  } else if(!is.null(indBW2)) {
		    if(!all(indBW2==0)) keepInd[indrows] = FALSE 
		  }
		}
		
	 } #end marker
	} #end sample
	return(dat[keepInd,,drop=FALSE])
}