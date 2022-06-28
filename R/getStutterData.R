#' @title getStutterData 
#' @description helpfunction to obtain stutter data table used for making inference
#' @details Does not require getFilteredData first anymore
#' The nomenlclature for stutters are as follows:
#' BW1 means that the LUS motif is providing a backward stutter (n-1).
#' FW1 means that the LUS motif is providing a forward stutter (n+1).
#' BW2 means that the non-LUS motif is providing a backward stutter (n-1).
#' etc... prefix "D" means double (n-2)
#' n0 stutters are named as either FWBW or BWFW depending on which occuring first in the sequence (hence orientation matters!)
#' 
#' @param dat Datatable with columns (RefName,SampleName,Locus,Allele,Coverage,Dose)
#' @param minStuttOccurence The number of required observations per stutter type (remaining put to noise)
#' @param platform Platform type {"MPS","CE"}. NOTE: CE has only consider only some stutter types (but also include -2bp)
#' @param stutterTypes Stutter types. See getStutterTypeRule for supported types
#' @param verbose Whether to print out progress
#' @export 

getStutterData = function(dat, minStuttOccurence=5, platform="MPS",stutterTypes=NULL,verbose=FALSE) {
  sep = "/"
  if(!platform%in%c("CE","MPS")) stop("Platform not supported. Please use either CE or MPS!")
  isCE = platform=="CE"
  isMPS = platform=="MPS"
  
  #ASSUME COMMON STUTTER TYPES FOR EACH PLATFORM
  if(is.null(stutterTypes)) {
    stutterTypes=c("BW1","FW1","DBW1","DFW1","TBW1") #also for CE
    if(isMPS) {
      stutterTypes=c(stutterTypes,"BW2","FWBW","BWFW") #exclusive for MPS
    } else if(isCE) {
      stutterTypes = substr(stutterTypes, 1,nchar(stutterTypes)-1) #remove rank letters
      stutterTypes = c(stutterTypes,"2BP") #-2BP stutters (important for SE33)
    }
  }
  
  #Extracting unique samples and loci
  samples = unique(dat$SampleName) #obtain unique samples
  locs = unique(dat$Locus) #obtain unique loci
  nLocs = length(locs) #number of loci
  nSamples = length(samples) #number of samples
  
  #Define column name for output tables
  cn = c("Locus","Sample","Type","Allele","Coverage","Parental","Ratio","Prop","BLMM")  
  
  if(verbose) print("Performing indexing for different stutter types...")
  stuttertable <- noisetable <-  NULL #output tables
  for(loc in locs) { #traversing each locus
#    loc=locs[3]
    if(verbose) print(loc)
    for(sample in samples) {
# sample=samples[2]
      sub = dat[dat$Locus==loc & dat$SampleName==sample,,drop=F]
      
      donor = sub$Allele[sub$Dose>0] #donor allele
      isArtefact = sub$Dose==0 #get index of artefact
      artefacts = sub$Allele[isArtefact] #artefacts
      
      if(length(donor)<2 || length(donor)==nrow(sub)) next  #skip if homozygous OR all alleles are donors (no artefacts)
      if(length(donor)>2) stop( paste0("Sample=",sample," Locus=",loc,": Not possible with more than 2 true alleles") )
      
      #STEP 1: Obtain canidate stutter types for each artefacts
      candStuttType = rep(NA,length(artefacts)) #get candidate stutter types for each artefact
      candBLMM <- candParental <- candStuttType
      for(type in stutterTypes) {
        for(donor1 in donor) { #traverse one donor allele at the time
    #      donor1=donor[1]
          stuttObj = MPSproto::getStutteredSequence(donor1,type) #get expected stutter
          stuttAllele = stuttObj$stutter #get stutter
          if(is.na(stuttAllele)) next
          indStutt = stuttAllele==artefacts #get which could be parental allele
          if(!any(indStutt)) next
          
          isNAins = is.na(candStuttType[indStutt])
          
          indStutt1 = indStutt[isNAins] #non inserted before
          indStutt2 = indStutt[!isNAins] # already included before
          
          candStuttType[indStutt1] = type
          candStuttType[indStutt2] = paste0(candStuttType[indStutt2],sep,type)

          candBLMM[indStutt1] = stuttObj$BLMM
          candParental[indStutt1] = donor1
        }
      }		  

      
      #Alleles not having any parental allele are indicated as Noise
      isNoise = is.na(candStuttType)
      if(any(isNoise)) {
        isNoise2 = which(isArtefact)[isNoise]
        noiseAlleles = sub$Allele[isNoise2]

        newrows = data.frame(loc,sample,"NOISE",noiseAlleles,sub$Coverage[isNoise2],NA,NA,NA,NA)
        colnames(newrows) = cn
        
        #Check if the noise allele is a 1bp error
        donor_seqs = MPSproto::convertBracket2seq(donor)
        noise_seqs = MPSproto::convertBracket2seq(noiseAlleles) #obtain sequence of noise
        nchanges = adist(noise_seqs,donor_seqs)  #obtain sequence of donors
        is1errorTab = which(nchanges==1,arr.ind = TRUE) #check which is 1 bp error
        for(i in seq_len(nrow(is1errorTab)) ) {
          inds = is1errorTab[i,]
          parental = donor[inds[2]]
          ratio = newrows$Coverage[inds[1]]/sub$Coverage[sub$Allele==parental] #calc ratio
          
          newrows$Type = "ERR1" #one base error
          newrows$Parental[inds[1]] = parental
          newrows$Ratio[inds[1]] = ratio
          newrows$Prop[inds[1]] = 1/(1+1/ratio) #convert to proportion
          
          #get broken part (BLMM)
          allele1 = MPSproto::getMotifRepsSequences(newrows$Allele)[[1]] #errored
          allele2 = MPSproto::getMotifRepsSequences(parental)[[1]] #parental
          for(j in seq_along(allele2)) {
            blmm = allele2[j] #obtain motif length
            if(blmm!=allele1[j]) break #found if repat length was not the same
          }
          newrows$BLMM = blmm
        }
        
        noisetable = rbind(noisetable, newrows)
      }
      
      
      isOKstutter = !grepl(sep,candStuttType) & !isNoise
      isOKstutterInd = which(isOKstutter)
      
      #Extract stutter allele info
      for(ind in isOKstutterInd) { #traverse each OK stutter
        ToInd = which(isArtefact)[ind] #get index of stutter
        
        if(sub$Allele[ToInd]!=artefacts[ind]) stop() # must be true
        type = candStuttType[ind]
        FromBLMM = candBLMM[ind] 
        Parental = candParental[ind] #obtain perantal allele
        
        FromInd = Parental==sub$Allele
        coverage_stutter = sub$Coverage[ToInd] #this is coverage of stutter
        stutter_rate = sub$Coverage[ToInd]/sub$Coverage[FromInd]
        stutter_prop = 1/(1/stutter_rate + 1) #convert to stutter prop
        newrow = data.frame(loc,sample,type,artefacts[ind],coverage_stutter,Parental,stutter_rate,stutter_prop,FromBLMM)
        colnames(newrow) = cn
        stuttertable  = rbind(stuttertable,newrow)
      } #end for each stutter type
    } #end for each sample
  } #end for each locus
  
  #LAST: IF TOO FEW STUTTER TYPES OBTAINED FOR A MARKER
  #THESE ARE MOVED TO NOISE TABLE
  if(minStuttOccurence>0) {
    if(verbose) print("Attaching too few stutter observations to noise data...")
    #create a overview table:
    tab = table(stuttertable$Locus,stuttertable$Type)
    mustMove = which(tab <= minStuttOccurence & tab>0,arr.ind = TRUE)
    for(i in seq_len( nrow(mustMove)) ) {
      stuttType = colnames(tab)[mustMove[i,2]]
      stuttLoc = rownames(tab)[mustMove[i,1]]
      
      #INDICATE
      moveAsNoise = stuttertable$Locus==stuttLoc & stuttertable$Type==stuttType
      
      #ADD TO NOISE
      noisetable = rbind(noisetable,stuttertable[moveAsNoise,,drop=FALSE])
      
      #ADD TO NOISE
      stuttertable = stuttertable[!moveAsNoise,,drop=FALSE] #update by removing
    }
  }
  
  return(list(stutterdata=stuttertable,noisedata=noisetable))
}
