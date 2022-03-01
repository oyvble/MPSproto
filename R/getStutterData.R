#' @title getStutterData 
#' @description helpfunction to obtain stutter data table used for making inference
#' @details The nomenlclature for stutters are as follows:
#' BW1 means that the LUS motif is providing a backward stutter (n-1).
#' FW1 means that the LUS motif is providing a forward stutter (n+1).
#' BW2 means that the non-LUS motif is providing a backward stutter (n-1).
#' etc... prefix "D" means double (n-2)
#' n0 stutters are named as either FWBW or BWFW depending on which occuring first in the sequence (hence orientation matters!)
#' 
#' @param dat Datatable with columns (RefName,SampleName,Locus,Allele,Coverage,Dose)
#' @param stutterTypes Stutter types. Supports (BW1,BW2,FW1,FW2,DBW1,TBW1,FWBW,BWFW). Also special -2bp stutters for CE
#' @param platform Platform type {"MPS","CE"}. NOTE: CE has only consider only some stutter types (but also include -2bp)
#' @param verbose Whether to print out progress
#' @export 

#platform="MPS";verbose=TRUE;stutterTypes=NULL
#platform="CE";verbose=TRUE;stutterTypes=NULL
getStutterData = function(dat,stutterTypes=NULL,platform="MPS",verbose=FALSE) {

  if(!platform%in%c("CE","MPS")) stop("Platform not supported. Please use either CE or MPS!")
  isCE = platform=="CE"
  isMPS = platform=="MPS"
  if(is.null(stutterTypes)) {
  	stutterTypes=c("BW1","FW1","DBW1")#,"DFW1","TBW1") #also for CE
      if(isMPS) {
        stutterTypes=c(stutterTypes,"BW2","FW2","FWBW","BWFW") #exclusive for MPS
      } else if(isCE) {
      stutterTypes = substr(stutterTypes, 1,nchar(stutterTypes)-1)
  	  stutterTypes = c(stutterTypes,"2BP") #-2BP stutters (important for SE33)
  	}
  }
  
  #Extracting unique samples and loci
  samples = unique(dat$SampleName) #obtain unique samples
  locs = unique(dat$Locus) #obtain unique loci
  nLocs = length(locs) #number of loci
  nSamples = length(samples) #number of samples
  
  if(verbose) print("Performing indexing for different stutter types...")
  datatable = NULL #create datatable for Observations and expectations, and stutter info
  noisetable = NULL
  for(loc in locs) { #traversing each locus
#    loc=locs[5]
    if(verbose) print(loc)
    for(sample in samples) {
# sample=samples[2]
      sub = dat[dat$Locus==loc & dat$SampleName==sample,,drop=F]
      if( nrow(sub)==0 || !any(sub$Dose==0) ) next #SKIP IF NO OBSERVED STUTTERS

#      stop()
      seqs = sub$Allele #extract all alleles (should be string)
      #sub$Coverage
      SIlist <- BLMMlist <- list() #store stutter index info and BLMM info in list
      for(type in stutterTypes) {
#        type = stutterTypes[3]
#        type="2BP" 
        stuttInfo = getStutterIndex(seqs=seqs,type=type,platform=platform, locus=loc) #obtain stutter indices
        if(is.null(stuttInfo)) next #skip stuttertype 
        
        SIlist[[type]] = stuttInfo$SI
        BLMMlist[[type]] = stuttInfo$BLMM
      }    
      #data.frame(SIlist)
      #data.frame(BLMMlist)

      #NB: Assume that missing true alleles are already inserted in Dose information:
      trueAllele_ind = which(sub$Dose>0) #obtain index of true alleles
      isStuttAllele = setNames( sub$Dose==0,seqs) #obtain stutter alleles

      #Extract info where a stutter allele links to a true allele (for each stutter type):
      subtable = NULL
      for(type in names(SIlist)) {
#   type=names(SIlist)[4]
        #OBTAIN CANDIDATES
        ToInd = which(isStuttAllele) #stutter to - index
        FromInd = SIlist[[type]][isStuttAllele] #stutter from -index
        BLMMstutt = BLMMlist[[type]][isStuttAllele] #obtain BLMM of stuttered
        #seqs
        #sub2 = subset(sub,select=c(Allele,Coverage))
        #View( cbind(seq_len(nrow(sub2)),sub2, FromInd[match(sub2$Allele,names(FromInd))] ) )
        #if(any(grepl("/",FromInd))) stop("Multiple From index")
  
        #STORE ONLY WHICH HAS FROM ALLELE AS PARENtAL ALLELE
        isOK = rep(FALSE,length(FromInd))
        for(i in seq_along(isOK) ) {
#i=4          
          stuttFromVec = unlist(strsplit(as.character(FromInd[i]),"/")) #obtain index which allele stutter from
          blmmVec = unlist(strsplit(as.character(BLMMstutt[i]),"/")) #obtain index which allele stutter from
          isFromTrueAllele = stuttFromVec%in%trueAllele_ind  #obtain which are from true allele
          if(sum(isFromTrueAllele)==1) { #possible if one solution parental (Store this)
            FromInd[i] = stuttFromVec[isFromTrueAllele]
            BLMMstutt[i] = blmmVec[isFromTrueAllele] #store correct BLMM
            isOK[i] = TRUE
          }
        }      
        if(sum(isOK)==0) next #skip if none is OK
        ToInd = as.integer(ToInd[isOK])
        FromInd = as.integer(FromInd[isOK])
        FromBLMM = BLMMstutt[isOK] #store BLMM of parental allele
        
        coverage_stutter = sub$Coverage[ToInd] #this is coverage of stutter
        stutter_rate = sub$Coverage[ToInd]/sub$Coverage[FromInd]
        stutter_prop = 1/(1/stutter_rate + 1) #convert to stutter prop
        newrow = data.frame(loc,sample,type,seqs[ToInd],seqs[FromInd],FromBLMM,coverage_stutter,stutter_rate,stutter_prop)
        colnames(newrow) = c("Locus","Sample","Type","Stutter","Parental","BLMM","Coverage","Rate","Prop")
        subtable  = rbind(subtable,newrow)
      } #end for each stutter type
      datatable  = rbind(datatable,subtable)
      
      #Alleles not being parental allele or stutter are indicated as Noise
      isNoise = !seqs%in%unique(c(subtable$Stutter,subtable$Parental, seqs[trueAllele_ind]))
      if(any(isNoise)) {
        newrow = data.frame(loc,sample,seqs[isNoise],sub$Coverage[isNoise])
        colnames(newrow) = c("Locus","Sample","Sequence","Coverage")
        noisetable = rbind(noisetable, newrow)
      }
    } #end for each sample
  } #end for each locus
  
  return(list(stutterdata=datatable,noisedata=noisetable))
}
