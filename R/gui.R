#' @title gui
#' @author Oyvind Bleka
#' @description gui is a GUI wrapper (gWidgets) for the MPSproto package
#' @param envirfile A Rdata file including a saved environment of a project
#' @param envir A environment object
#' @export


#library(MPSproto);envirfile=NULL;envir=NULL
gui = function(envirfile=NULL,envir=NULL) {
  #size of main window
  mwH <- 500
  mwW <- 1000
  
  #type of gwidgets-kit
  #library(gWidgetstcltk)
  options(guiToolkit="tcltk")
  
  pack = "MPSproto" #name of package and program
  #version:
  version = packageVersion(pack) #follows same version as package number
  softname <- paste0(pack," v",version)
  
  #GUI-Restriction on maximum number of contributors
  maxKsetup <- 4 
  CEsep = ":" #CE separator
  
  #Spacing between widgets
  spc <- 10
  .sep <- .Platform$file.sep # Platform dependent path separator. 
  pgkPath <- path.package(pack, quiet = FALSE) # Get package path.
  
  #Option files are stored in system (opt-settings): These are default values
  optList = list(
    optFreq =list(fst=0,freqsize=0,minF=0,normalize=1), #option when new frequencies are found (size of imported database,minFreq), and missmatch options
    optMLE = list(nDone=3,delta=1,maxIter=100,maxThreads=0,seed=0,steptol=1e-3,equaltol=0.01), #options when optimizing,validation (nDone,delta)
    optDC = list(alphaprob=0.99,maxlist=20), #options when doing deconvolution
    optKit = list(platform="MPS",kit="NONE"),
    optFreqFile =  NULL, #default frequency file 
    optCalibFile = NULL #default calibration (NEW)
  )  
  
  configList = list() #list of file name for each element

  #TRAVERSE EACH config FILE (configPrefix) AND CHECK IF CONFIG FILE XISTS
  for(optName in names(optList) ) {
    tmp <- paste0("config",gsub("opt","",optName))
    configFile <- configList[[optName]] <- paste0(pgkPath,.sep,tmp)
    if( file.exists(configFile) ) {  #use default values if not existing
      #optF <- scan(file=configFile,what=character(),quiet=TRUE)
      optF <- readLines(configFile,warn = FALSE) #,what=character(),quiet=TRUE)
      for(i in 1:length(optF)) {
        suppressWarnings({
          tmp = as.numeric(optF[i])
        })
        if(is.na(tmp)) tmp = optF[i] #keep original
        if(length(tmp)==0) tmp = NULL #no value
        optList[[optName]][[i]] = tmp #convert direclty from file
      } 
      if(length(optF)==1) optList[[optName]] = optList[[optName]][[1]]
    }
      
  }

  
  #####################
  #create environment #
  #####################
  
  #Helpfunction to obtain the calibration object 
  .getCalibrated = function(fn) {
    calib = readRDS(fn)  #load object
    return(calib)
  }
  
  if(!is.null(envir)) {
    mmTK = envir #use environment direclty
  } else if(!is.null(envirfile)) {
    load(envirfile) #loading environment
  } else { #NEW OBJECT
    mmTK = new.env( parent = globalenv() ) #create new envornment object (must be empty)
  
    #Toolbar options: can be changed any time by using toolbar
    for(optName in names(optList) ) assign(optName,optList[[optName]],envir=mmTK)  #store to envir

    #initializing environment variables
    assign("workdir",NULL,envir=mmTK) #store work dir
    assign("version",NULL,envir=mmTK) #Store version
    
    #Check frequency file and try import it.:    
    popFreq = NULL #try load frequency data from selected file
    if(!is.null( optList$optFreqFile ) ) {
      tryCatch({
        popFreq = MPSproto::importMPSfreqs( optList$optFreqFile)[[1]] 
      }, error=function(e) print("System frequency file not valid."))
    }
    assign("popFreq",popFreq,envir=mmTK)  #store frequency in 

    #Check calibration file and try import it.:    
    calibrated = NULL #try load calibration data from selected file
    if(!is.null( optList$optCalibFile ) ) {
      tryCatch({
        calibrated = .getCalibrated( optList$optCalibFile )
      }, error=function(e) print("System calibration file not valid."))
    }
    assign("calibrated",calibrated,envir=mmTK) 
    
        
    #imported data:
    assign("mixData",NULL,envir=mmTK)  #DATA FOR MIXTURES
    assign("refData",NULL,envir=mmTK)  #DATA FOR REFERENCES
    assign("relSettings",NULL,envir=mmTK)  #this is settings for relatedness (NOT IMPLEMENTED)
        
    #models: stored setups for model specification
    assign("calcList",NULL,envir=mmTK) #store model of calculated results  (stored in a table)
    #STRUCTURE: calcList[[index]] = list(meta,hd,hp) #hp can be NULL if only deconvolution (POI leaved empty)
    
    #results: stored results after calculations
    #assign("resEVID",NULL,envir=mmTK) #assign evidence weighting results (i.e. calculated LR with MLE estimates)
    assign("resDC",NULL,envir=mmTK) #assign deconvolved results (i.e. ranked tables of results)
    assign("resCompare",NULL,envir=mmTK) #assign information when doing model selection searchs (NOT USED)
  }

  ####################################
  #auxiliary functions and variables:#
  ####################################
  

  #helpfunction to get small number from log-value
  .getSmallNumber = function(logval,sig0=2,scientific="e") {
    log10base = logval/log(10) #convert to 10 base
    power = floor(log10base) #get power number
    remainder = log10base - power
    return( paste0( round(10^remainder,sig0),scientific,power)) #representation of very small numbers (avoid underflow)
  }
  
   .NAtoSign <- function(x) {
    x[is.na(x)] <- "-" #NB: New version of gtable does not accept NA values
    return(x)
   }
  # helptext = function(obj,txt) { gWidgets2::addHandlerRightclick(obj,handler=function(h,...) { gWidgets2::gmessage(txt,title="Detailed information") }) }
  .helptext = function(obj,txt) { gWidgets2::tooltip(obj) = txt } 
  
  #Function to get data from environment
  #sel used to select a specific datasubset
  .getData = function(type,sel=NULL) {
   Data <- NULL
   if(type=="mix") Data <- get("mixData",envir=mmTK) #assign kit to mmTK-environment
   if(type=="ref") Data <- get("refData",envir=mmTK) #assign kit to mmTK-environment 
   if(!is.null(sel)) return(Data[sel]) #returns only selected datasubset
   return(Data)
  }
  
  #function for inserting sample/ref/db-names into existing gWidgets2::gcheckboxgroup
  .getDataNames_type = function(type) {
    subD <- .getData(type)
    if(!is.null(subD)) { return( names(subD))
    } else { return("") }
  }
  
  #Function which takes rownames and adds to first column
  .addRownameTable = function(tab,colname=NULL) {
    tmp <- colnames(tab)
    tab <- cbind(rownames(tab),tab)
    if(is.null(colname)) colname = "X"
    colnames(tab) <- c(colname,tmp)
    return(tab)
  }
  
  #save result table to file:
  .saveTable = function(tab,sep="txt") {
    tabfile  = .mygfile(text="Save table",type="save") #csv is correct format!
    if(length(tabfile)==0) return()
     if(length(unlist(strsplit(tabfile,"\\.")))==1) tabfile = paste0(tabfile,".",sep)
     if(sep=="txt" | sep=="tab") write.table(tab,file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) 
     if(sep=="csv") write.table(tab,file=tabfile,quote=FALSE,sep=",",row.names=FALSE) 
     print(paste("Table saved in ",tabfile,sep=""))
  } #end file
  
  #Helpfunction to print tables (R v3.5.0) has problem showing tables in the GUI.
  .printTable = function(x) {
    print(cbind(rownames(x),x),row.names=FALSE)
  }

  
###################################################################
#############GUI HELPFUNCTIONS#####################################
###################################################################

 #TAKEN FROM CASESOLVER: This function is written since the encoding in  gWidgets2::gfile is fixed to UTF-8 which doesn't handle special letters
 .mygfile <- function(text,type,filter=list(),initf=NULL) { #Possible bug: text ignored when type="selectdir"
   file <- gWidgets2::gfile(text=text,type=type,filter=filter,initial.filename=initf)
   Encoding(file) <- options()$encoding #Set to local encoder: Handle special cases.
   return(file)
 }
 
 #Helpfunction to get focus 
 .getFocus = function() {
   gWidgets2::visible(mainwin) <- TRUE
   gWidgets2::focus(mainwin) <- TRUE
 }
 
 #Menu bar file-lists:
 .f_setwd = function(h,...) {
  dirfile = .mygfile(text="Select folder",type="selectdir")
  if(length(dirfile)==0) return()
  setwd(dirfile)
  assign("workdir",dirfile,envir=mmTK) #assign working directory
 }
  
 .f_openproj = function(h,...) {
  projfile = .mygfile(text="Open project",type="open", filter=list("Project"=list(patterns=list("*.Rdata")) , "All files"=list(patterns=list("*"))))
  if(length(projfile)==0) return()
  gWidgets2::dispose(mainwin)
  MPSproto::gui(projfile) #send environment into program
 }
 
 .f_saveproj = function(h,...) {
  projfile = .mygfile(text="Save project",type="save")
  if(length(projfile)==0) return()
  
   if(length(unlist(strsplit(projfile,"\\.")))==1) projfile = paste0(projfile,".Rdata") #add extension if missing
  # tmp = sort(sapply(mmTK,object.size)/1e6,decreasing=TRUE)
   #.selectDataToModel(h=list(action="SAVE")) #store data before saving project
   #save(mmTK,file=projfile,compress="xz") #save environment, #Dont compress!
   save(mmTK,file=projfile,compress="xz",eval.promises=FALSE,precheck=FALSE,compression_level=2)
   print(paste("Project saved in ",projfile,sep=""))
  
 }

 .f_quitproj = function(h,...) {
  ubool <- gWidgets2::gconfirm("Do you want to save project?",title="Quit Program",icon="info")
  if(ubool) {
    .f_saveproj(h)
  } else { 
   print("Program terminated without saving")
  }
  gWidgets2::dispose(mainwin) #remove window!
 }
 
 #Throw error message
 .NAerror <- function(what) {
   gWidgets2::gmessage(paste0(what," must be specified as a valid value"),title="Wrong input",icon="error")
   #stop("Wrong user-input")
 }
 
 #helpfunction to get value from user
 .getValueUser <- function(txt="",val=0) {
   val2 <- gWidgets2::ginput(txt, text=val, title="User input",icon="question")
   return(val2)   
 }

 #helpfunction to get value in from user and store
 #what1="optFreq";what2="fst"
 .setValueUser <- function(what1,what2,txt,allowNULL=FALSE,allowText=FALSE) {
   listopt <- get(what1,envir=mmTK) #get object what 1.
   val <- listopt[[what2]]
   if(is.null(val)) val ="" #gwidgets2 does not handle NULL, must use empty string
   sw <- gWidgets2::gwindow(title="User input",visible=FALSE, width=300,height=50)
   grid <- gWidgets2::glayout(spacing=0,container=sw )
   grid[1,1] <- gWidgets2::glabel(txt, container=grid)
   grid[1,2] <- gWidgets2::gedit(text=val,container=grid,width=30)
   grid[2,1] <- gWidgets2::gbutton("OK", container=grid,handler = function(h, ...) { 
    GUIval = gWidgets2::svalue(grid[1,2]) #obtain GUI value
    if(allowNULL && GUIval=="") { #if accepting empty string
      tmp = NULL #Insert NULL
    } else {
      tmp <- GUIval
      if(!allowText) {
        tmp <- as.numeric(GUIval) #insert new value
        if(is.na(tmp)) {
          .NAerror(what2)
          return()
        }
      }
    }
    listopt[[what2]] <- tmp
    assign(what1,listopt,envir=mmTK) #assign user-value to opt-list
    write(unlist(listopt),file=configList[[what1]],ncolumns = 1) #store selected values to file   
    gWidgets2::dispose(sw)
   } )
   grid[2,2] <- gWidgets2::gbutton("Cancel", container=grid,handler = function(h, ...) { gWidgets2::dispose(sw) } )
   gWidgets2::visible(sw) <- TRUE
 }

 #helpfunction to check value of x
 #replaces checkProb,checkPositive,checkPosInteger
 .checkValue = function(x,type,what,strict=FALSE) {
   if(is.na(x)) .NAerror(what)
   isOK = TRUE
   if(type=="prob" && (x < 0 || x>1)) isOK = FALSE
   if(type=="posInt" && (x < 1 || round(x)!=x)) isOK = FALSE
   if(type=="pos") {
     if(x < 0) isOK = FALSE
     if(strict && x==0) isOK = FALSE
   } 
   
   if(!isOK) {
     switch(type,
            prob={msg="must be specified in interval [0,1]"},
            pos={msg="must be a positive number"},
            posInt={msg="must be a positive integer"}
     )
     gWidgets2::gmessage(paste0(what," ",msg,""),title="Wrong input",icon="error")
     stop("Wrong user-input")
   }
 }

 #helpfunction for printing evidence sample to terminal
 .printEvid = function(subD) {
  locs <- names(subD) #get unique loci
  mixtab <- matrix(ncol=2,nrow=length(locs))
  for(loc in  locs) { #for each locus
        mixA <- subD[[loc]]$adata
        mixH <- subD[[loc]]$hdata
        if(!is.null(mixA)) mixtab[which(loc==locs),1] <- paste0(mixA ,collapse="/")
        if(!is.null(mixH)) mixtab[which(loc==locs),2] <- paste0(mixH ,collapse="/")
  }
  rownames(mixtab) <- locs
  colnames(mixtab) <- c("Allele","Height")
  .printTable(mixtab)
 }  
 
 #helpfunction for printing reference sample to terminal
 .printRefs = function(refD,refSel=NULL) {
   locs <- unique(unlist(lapply(refD,names))) #get unique loci
   reftab <- matrix(ncol=length(refSel),nrow=length(locs)) #last row is RMP
   for(rsel in refSel) {
    for(loc in  locs) { #for each locus
      refA <-refD[[rsel]][[loc]]$adata
      if(!is.null(refA)) {
       reftab[which(loc==locs),which(rsel==refSel)] <- paste0(refA ,collapse="/")
      }
    }
   }
   rownames(reftab) <- locs
   colnames(reftab) <- refSel 
   .printTable(reftab)
  }
 
 
########################################################################################################################################## 
  
##################################### 
###########GUI WINDOW STARTS#########
##################################### 
 
 ##########
 #Menu bar#
 ##########
 mblst = list( #project saving and so on
  File=list(  
    gWidgets2::gaction('Set directory',handler=.f_setwd),
    gWidgets2::gaction('Open project',handler=.f_openproj),
    gWidgets2::gaction('Save project',handler=.f_saveproj),
    #gWidgets2::gaction('Settings',handler=.f_settings, action="GLOBAL"),
    gWidgets2::gaction('Quit',handler=.f_quitproj,icon="close")
  ),
  Frequencies=list(
    gWidgets2::gaction('Set Fst of calculation',handler=function(h,...) {   #moved from settings
      .setValueUser(what1="optFreq",what2="fst",txt="Set Fst-coefficient") 
    }),
    gWidgets2::gaction('Set size of frequency database',handler=function(h,...) {  
      .setValueUser(what1="optFreq",what2="freqsize",txt="Set size of imported freq database \n(Min observed used if not spesified):") 
    }),
    gWidgets2::gaction('Set minimum frequency',handler=function(h,...) {  
      .setValueUser(what1="optFreq",what2="minF",txt="Set minimum freq for new alleles\n(Min observed used if not spesified):") 
    }),
    gWidgets2::gaction('Set whether to normalize frequencies',handler=function(h,...) {  
      .setValueUser(what1="optFreq",what2="normalize",txt="Should frequencies always add up to one\nafter including rare alleles? (1=YES,0=NO)") 
    })
  ),
  Optimization=list(
    gWidgets2::gaction('Set number of optimizations',handler=function(h,...) {  
      .setValueUser(what1="optMLE",what2="nDone",txt="Set required number of (identical) optimizations:") 
    }),
    gWidgets2::gaction('Set variation of randomizer',handler=function(h,...) {  
      .setValueUser(what1="optMLE",what2="delta",txt="Set variance of start point randomizer:") 
    }),
    gWidgets2::gaction('Set max number of iterations',handler=function(h,...) {  
      .setValueUser(what1="optMLE",what2="maxIter",txt="Set max number of iterations:") 
    }),
    gWidgets2::gaction('Set maximum threads for computation',handler=function(h,...) {  
      .setValueUser(what1="optMLE",what2="maxThreads",txt="Set max number of threads to be used in parallelisation\n(all used if zero):") 
    }),
    gWidgets2::gaction('Set seed of randomizer',handler=function(h,...) { 
      .setValueUser(what1="optMLE",what2="seed",txt="Set seed of randomizer\n(Not used if zero):",allowNULL=TRUE) 
   }),
   gWidgets2::gaction('Set accuracy of optimization',handler=function(h,...) { 
     .setValueUser(what1="optMLE",what2="steptol",txt="Set accuracy of optimization (steptol, see ?nlm):") 
   })
  ),
  Deconvolution=list(
    gWidgets2::gaction('Set required summed probability',handler=function(h,...) {  
      .setValueUser(what1="optDC",what2="alphaprob",txt="Set required summed posterior genotype-probability of list:") 
    }),
    gWidgets2::gaction('Set max listsize',handler=function(h,...) {  
      .setValueUser(what1="optDC",what2="maxlist",txt="Set size of maximum elements in deconvoluted list:") 
    })
  )
 )

##################################################################################################
########### Program starts #######################################################################
##################################################################################################

 #change working directory to the one stored in mmTK-environment
 wd=get("workdir",envir=mmTK) #assign working directory to mmTK-environment
 if(!is.null(wd) && dir.exists(wd)) setwd(wd)
 
 #Main window:
 mainwin <- gWidgets2::gwindow(softname, visible=TRUE, width=mwW,height=mwH)
 gWidgets2::gmenu(mblst,container=mainwin)
 nb = gWidgets2::gnotebook(container=mainwin)
 tabimport0 = gWidgets2::ggroup(horizontal=FALSE,spacing=10,container=nb,label="Data") #tab2: (imports all files)
 tabmodel = gWidgets2::glayout(horizontal=FALSE, spacing=spc,container=nb,label="Model") #tab3: specify model used in weight-of-evidence (INT/MLE) or in a Database search 
 tabMLE = gWidgets2::ggroup(horizontal=FALSE, spacing=spc,container=nb,label="Results")#,expand=T,fill=T) #results from MLE
 tabDC = gWidgets2::ggroup(horizontal=FALSE,spacing=spc,container=nb,label="Deconvolution") #results from a deconvolution


####################################################
###############Tab 1: Import Data:##################
####################################################

 #When program starts, import assumed model for EVID.

  editboxsize = 3 #length of boxes 
 #b) load/save profiles/database: Supports any filetype
 
  .f_importprof = function(h,...) {
    type=h$action #get type of profile
# type = "mix"
    #  proffile = .mygfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab"))))
    proffile = .mygfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*"))))
    if(length(proffile)==0) return() 
    
    Data = MPSproto::importMPSsample(proffile) #load profile (converts from table to list format)

    #MUST UTILIZE LUSstrR if sequences are given
    #converting to "CE:bracket" format
    
    #Special handling
    if(type=="ref") { #check that homozygote alleles are given twice
      txt2 <- "The number of alleles in a genotype must be 2."
      txt1 <- paste0("Only one allele was given for a genotype in a reference profile. ",txt2, " Hence the second allele was automatically set as the first allele.") 
      txt3 <- paste0("Too many alleles where given for a genotype in a reference profile. ",txt2) 
      miss <- FALSE #indicator whether hom. are missing
      
      for(kn in names(Data)) { #for each profile
       for(loc in names(Data[[kn]])) { #for each profile
         if( length(Data[[kn]][[loc]]$adata)>2 )  {
           gWidgets2::gmessage(txt3,"Wrong file-input",icon="error")
           break #breaking loop if wrong input
         }
         if( length(Data[[kn]][[loc]]$adata)==1) {
          Data[[kn]][[loc]]$adata <- rep(Data[[kn]][[loc]]$adata,2) #duplicated
          miss <- TRUE
         }
          Data[[kn]][[loc]]$hdata = NULL #Updated v2.1.0: Remove peak heights if these have been added for references.
       }
     }
     if(miss) gWidgets2::gmessage(txt1,"Warning",icon="info")
    } 
  
    #get already stored data:
    if(type=="mix") Data2 <- .getData("mix") #get data from mmTK-environment
    if(type=="ref") Data2 <- .getData("ref") #get data from mmTK-environment
  
    oldNames = names(Data2) #old samples
    if(is.null(Data2)) { #if no previous already there
      Data2 <- Data
    } else {
      for(kn in names(Data)) Data2[[kn]] <- Data[[kn]] #insert dataframe for each profile
    }
    if(type=="mix")  assign("mixData",Data2,envir=mmTK) #assign data to mmTK-environment
    if(type=="ref")  assign("refData",Data2,envir=mmTK) #assign data to mmTK-environment
    sampleNames = names(Data2) #Updated samplenames
    
    #Update table
    if(type=="mix")  tabimportC[1,1][] <- cbind(sampleNames)
    if(type=="ref")  tabimportC[1,2][] <- cbind(sampleNames)
  }


  #Compare references against each evidence
 .f_compare = function(h,...) {
   evidD = .getData("mix") #get selected references
   refD <- .getData("ref")
   
   nEvids = length(evidD)
   nRefs = length(refD)
   if(nEvids==0 || nRefs==0) {
     print("Missing data to compare...")
     return()
   }
   
   nLocs <- MAC <- matrix(0,nrow=nEvids,nRefs,dimnames = list(names(evidD),names(refD)))
   for(evid in names(evidD) ) { #for each selected evidence 
     subD <- evidD[[evid]] #selected samples
     locs <- names(subD) #obtain locs

     for(loc in locs) { #for each locus
        evidA = subD[[loc]]$adata #allele vector of evid
        if( length(evidA)==0) next #don't count missing markers
        
        if( any( grepl(":",evidA) )) {
          tmp = strsplit(evidA,":") #remove 
          evidA = sapply(tmp,function(x) x[2])
        } 
        for(ref in names(refD) ) { #for each selected evidence 
          refA <- refD[[ref]][[loc]]$adata
          if(is.null(refA) || length(refA)==0) next #skip if no data
          
          nLocs[evid,ref] = nLocs[evid,ref] + 1 #count locus
          MAC[evid,ref] <- MAC[evid,ref] + sum(refA%in%evidA) #count number of alleles
       }
     } 
   }
   matchrate <- MAC/(2*nLocs)
   nMissMatches = 2*nLocs - MAC
   print("---NUMBER OF MISSMATCHES FOR EACH REFERENCES/SAMPLES:")     
   .printTable(nMissMatches)
 }
     
 #prints evidence, references, EPG, databases and population frequencies
 .f_viewdata = function(h,...) {
  sep = "--------------------------------------------"
  help_gmessage = function(what) gWidgets2::gmessage(paste0("Please import and select ",what),icon="info")

  #View evidence:
  mixSel <- gWidgets2::svalue(tabimportC[1,1])  #get selected references
  refSel <- gWidgets2::svalue(tabimportC[1,2])  #get selected references
  
  evidD = .getData("mix",mixSel) #get selected references
  refD <- .getData("ref",refSel)
  nrefs = length(refD) #number of refs
  
  #PRINT TO CONSOLE:
  for(evidName in names(evidD) ) {
    subD <- evidD[[evidName]] #selected samples
    print("------------------------------------")
    print(paste("Samplename: ",evidName,sep=""))
    .printEvid(subD)
  }
  if(nrefs>0) .printRefs(refD,refSel)
  
  #FINALLY: SHOW EPGs (but first, modify allles of referenes)
  if(length(mixSel)>0) {
    refD2=refD #modify alleles of refernces
    if(nrefs>0) { #if any referenes
      CEsep = ":"
      locs = unique(unlist(lapply(refD,names))) #get loci to consider
      for(loc in locs) {
        #loc=locs[1]
        alleles = unique( unlist(lapply(evidD,function(x) x[[loc]]$adata) ))
        if( !any(grepl(CEsep,alleles)) ) next
        convtab = matrix(unlist(strsplit(alleles,CEsep)),nrow=2) #convert table
        
        for(ref in refSel) {
          #ref  = refSel[1]
          av <- unlist(refD2[[ref]][[loc]]) 
          found = av%in%convtab[2,]
          ind = match(av[found],convtab[2,]) #get col inds
          av[found] = paste0(convtab[1,ind],CEsep,av[found])
          refD2[[ref]][[loc]]$adata = av #insert again
        }
      }
    } 
    MPSproto::plotMPS(evidD,refD2)
  } 
  
 }  #end viewdata
 

 .f_viewcalib = function(h,...) {
   calibrated = get("calibrated",envir=mmTK) #assign calibrated object
   locs = names(calibrated)
   cns = unique(unlist(sapply(calibrated,names)))
   tab = NULL
   for(loc in locs) {
     row = rep("",length(cns))
     for(c in seq_along(cns)) {
       elem = calibrated[[loc]][[cns[c]]] #extract
       if(!is.null(elem)) row[c] = paste0(round(elem,3),collapse="/")
     }
     tab = rbind(tab,row)
   }
   rownames(tab) = locs
   colnames(tab) = cns
   print(tab)
   
   dbwin <- gWidgets2::gwindow("Calibrated model", visible=FALSE)#,height=mwH)
   gWidgets2::gtable( .addRownameTable(tab) ,container=dbwin) #create table
   gWidgets2::visible(dbwin) <- TRUE
   
 }
 
 #helpfunction used to extract selected importdata-elements to further model-setup
 .selectDataToModel <- function(h,....) {
   #All: popFreq must be imported!
   #EVID: must have both mixture and reference profiles
   #DC: Deconvolution requires only mixtures. Reference profiles is optional
   popFreq <- get("popFreq",envir=mmTK)
   mixSel <- refSel <-  numeric()
   if(length(tabimportC[1,1][])>0) mixSel <- gWidgets2::svalue(tabimportC[1,1])  #get selected mixtures
   if(length(tabimportC[1,2][])>0) refSel <- gWidgets2::svalue(tabimportC[1,2])  #get selected references
   
   if(is.null(popFreq)) {
     gWidgets2::gmessage("No frequencies was specified!\n Please import table with population frequencies.")
   } else if(length(mixSel)==0) {
     gWidgets2::gmessage("Please import and select evidence-profile!")
   } else {
     refreshTabModel(mixSel,refSel) #refresh table with selected data
     gWidgets2::svalue(nb) <- 2 #change tab of notebook
   }
 } #end selectDataToModel
 
 #Choose box and import button
 .setAsImported = function(fn,widget) { #set as imported
   gWidgets2::svalue(widget) <- "selected"
   gWidgets2::font(widget) <- list(weight="bold",size=11,color="green")
   .helptext(widget,fn) #set file name in .helptext
 }
 
 #helpfunction to store changed kit or platform type to environment and to file
 .f_selectedKit = function(h,...) {
   optKit = get("optKit", envir=mmTK) #obtain stored object
   if(h$action=="platform") optKit$platform = gWidgets2::svalue(tabimportA[3,2])
   if(h$action=="kit") optKit$kit = gWidgets2::svalue(tabimportA[4,2])
   assign("optKit",optKit,envir=mmTK) #assign to environment
   write(unlist(optKit),file=configList$optKit)   #STORE TO SYSTEM FILE
 }
 
 ###############
 #start layout:#
 ###############
 tabimport = gWidgets2::ggroup(container=tabimport0,horizontal = FALSE)
 #tabimport2 = gWidgets2::ggroup(container=tabimport0)
 
 tabimportA = gWidgets2::glayout(spacing=3,container=gWidgets2::gframe("Step 1) Import model data",container=tabimport)) #kit and population selecter
 tabimportB = gWidgets2::glayout(spacing=3,container=gWidgets2::gframe("Step 2) Import Profiles",container=tabimport)) #kit and population selecter
 tabimportC = gWidgets2::glayout(spacing=3,container=gWidgets2::gframe("Step 3) Select Evidence(s) and Reference(s)",container=tabimport,expand=T,fill=T),expand=T,fill=T) #evidence,ref dataframe
 tabimportD = gWidgets2::glayout(spacing=3,container=gWidgets2::gframe("Step 4) Select Interpretation",container=tabimport)) #Tasks button
 
 #SELECTION OF POPULATION FREQUENCIES 
 tabimportA[1,1] = gWidgets2::gbutton(text="Import frequencies",container=tabimportA,handler=
  function(h,...) {
    fn = .mygfile(text="Select file",type="open",filter = list(`All files` = list(patterns = c("*"))))
    if(length(fn)==0) return()
    popFreq = MPSproto::importMPSfreqs(fn)[[1]]
    
    assign("popFreq",popFreq,envir=mmTK) #assign popFreq
    assign("Freqfile",fn,envir=mmTK) #store filename
    .setAsImported(fn,tabimportA[1,2] ) #set label as imported
    write(fn,file= configList$optFreqFile )   #STORE TO SYSTEM FILE
  })
 .helptext(tabimportA[1,1],paste0("Choose a frequency file (LRmix/EFM format)."))
 tabimportA[1,2] <- gWidgets2::glabel("none", container=tabimportA) #Create label of whether freq data is imported:
 if(!is.null(get("popFreq",envir=mmTK))) .setAsImported( optList$optFreqFile, tabimportA[1,2]  )

 #SELECTION OF CALIBRATED OBJECT
 tabimportA[2,1] = gWidgets2::gbutton(text="Import calibrated",container=tabimportA,handler=
  function(h,...) {
    fn = .mygfile(text="Select file",type="open",filter = list(`All files` = list(patterns = c("*"))))
    if(length(fn)==0) return()
    calibrated = .getCalibrated(fn)
    assign("calibrated",calibrated,envir=mmTK) #assign calibrated object
    assign("calibFile",fn,envir=mmTK) #store filename
    .setAsImported(fn,tabimportA[2,2]) #set label as imported
    write(fn,file= configList$optCalibFile )   #STORE TO SYSTEM FILE
  })
 .helptext(tabimportA[2,1],paste0("Choose a calibration file (Rds format)."))
 tabimportA[2,2] <- gWidgets2::glabel("none", container=tabimportA) #Create label of whether calibration data is imported:
 if(!is.null(get("calibrated",envir=mmTK))) .setAsImported( optList$optCalibFile, tabimportA[2,2]  )

 #ADDING VIEW BUTTONS
 tabimportA[1,3] = gWidgets2::gbutton(text="View",container=tabimportA,handler=function(h,...) { print(get("popFreq",envir=mmTK)) })
 tabimportA[2,3] = gWidgets2::gbutton(text="View",container=tabimportA,handler=.f_viewcalib) 
 
 #ADDING platform and kit selections
 tabimportA[3,1] = gWidgets2::glabel("Platform",container=tabimportA)
 selPlatform = get("optKit", envir=mmTK)$platform
 tabimportA[3,2] = gWidgets2::gradio(c("MPS","CE"),select=selPlatform, container=tabimportA,horizontal = TRUE, handler=.f_selectedKit, action="platform")

 tabimportA[4,1] = gWidgets2::glabel("Kit:",container=tabimportA)
 listKits <- c("NONE",MPSproto::getMPSkit()) #list kit names (shortname)
 selKit = which(listKits==get("optKit", envir=mmTK)$kit)
 if(length(selKit)==0) selKit = 1 #first by default
 
 tabimportA[4,2] = gWidgets2::gcombobox(items=listKits, width=100, selected = selKit , editable = FALSE , container = tabimportA, handler=.f_selectedKit, action="kit")
   
 #Set other import buttons
 tabimportB[1,1] = gWidgets2::gbutton(text="Import evidence",container=tabimportB,handler=.f_importprof,action="mix")
 .helptext(tabimportB[1,1],"Imports evidence  profile(s) from a selected file. \n\nThe column names must contain 'sample..', 'marker', 'allele..', 'height..'.")
 
 tabimportB[1,2] = gWidgets2::gbutton(text="Import reference",container=tabimportB,handler=.f_importprof,action="ref")
 .helptext(tabimportB[1,2],"Imports reference profile(s) from a selected file. \n\nThe column names must contain 'sample..', 'marker', 'allele..'.")
 
 
 tabimportC[1,1] = gWidgets2::gcheckboxgroup(items= .getDataNames_type("mix"), container = tabimportC, use.table=TRUE)
 tabimportC[1,2] = gWidgets2::gcheckboxgroup(items= .getDataNames_type("ref"), container = tabimportC, use.table=TRUE)
 
 gWidgets2::size(tabimportC[1,1]) <- gWidgets2::size(tabimportC[1,2]) <- c(450,150) #c(600,150)
 
 #Button-choices further (INTERPRETATION:
 tabimportD[1,1] = gWidgets2::gbutton(text="View data",container=tabimportD,handler=.f_viewdata)
 .helptext(tabimportD[1,1],"Print all selected data to console. Selected Evidence profiles are shown as EPGs (replicates for those in same group).")

 tabimportD[1,2] = gWidgets2::gbutton(text="Compare data",container=tabimportD,handler=.f_compare)
 .helptext(tabimportD[1,2],"Obtain number of missmatches between all references and evidences.")
  
 tabimportD[1,3] = gWidgets2::gbutton(text="Interpretation",container=tabimportD,handler=.selectDataToModel)
 .helptext(tabimportD[1,3],"A module for calculating the Likelihood Ratio or Deconvolution for selected sample(s) (treated as replicates).")
 
 tabimportD[1,5] = gWidgets2::gbutton(text="RESTART",container=tabimportD,handler=function(h,...) {
   gWidgets2::dispose(mainwin) #remove window!
   MPSproto::gui() # open gui again
 })
 
 
####################################################################################################################
#######################################Tab 2: Hypotheses:##########################################################
#####################################################################################################################

#  mixSel= names(.getData("mix"));refSel=names(.getData("ref"));
  refreshTabModel = function(mixSel,refSel) { 
    booltxt = c(TRUE,FALSE)
    names(booltxt) = c("Yes","No")
    
    nRefs = length(refSel)
    nEvids = length(mixSel)
    mixData2 = get("mixData",envir=mmTK)[mixSel]
    refData2 = NULL
    if(nRefs>0) refData2 = get("refData",envir=mmTK)[refSel]
    
    #Obtain all alleles:
    allAlleles <- unlist(lapply(mixData2,function(x) unlist(lapply(x,function(y) y$adata)))) #get allele names of sample
    
    #helpfunction which takes GUI settings and stores them in "set'type'"
    storeSettings = function() {
      
      #EXTRACT HYPOTHESIS SPECIFICATION (tabmodelA1, tabmodelA2)
      #get specified preposition 
      NOC = as.integer(gWidgets2::svalue(tabmodelA1[1,2])) #number of contributors is fixed! 
      .checkValue(NOC,"posInt","Number of contributors under Hp/Hd")
      nC_hp <- nC_hd <- NOC #number of contributors in model:
      
      condOrder_hp <- condOrder_hd <- rep(0,nRefs)
      knownref_hp <- knownref_hd <- NULL #Typed profiles (known non-contributors)
      for(row in seq_along(refSel)) { #traverse each reference
        valhd <- as.integer(gWidgets2::svalue(tabmodelA2$hd[row,2])) #Hd
        condOrder_hd[row] <- valhd + valhd*max(condOrder_hd)
        valhp <- as.integer(gWidgets2::svalue(tabmodelA2$hp[row,2]))  #Hp
        condOrder_hp[row] <- valhp +  valhp*max(condOrder_hp)
      }
      knownref_hp <- which(condOrder_hp==0) #those not conditioned on under Hp
      if(length(knownref_hp)==0) knownref_hp <- NULL
    
      #THIS IS DONE FOR BOTH DC AND EVID
      knownref_hd <- which(condOrder_hd==0) #those not conditioned on under Hd
      if(length(knownref_hd)==0) knownref_hd <- NULL

      #get specified preposition 
      nHp =  sum(condOrder_hp>0) #number of conditionals
      nHd =  sum(condOrder_hd>0) #number of conditionals
      if( nHp > NOC || nHd > NOC) {
        gWidgets2::gmessage("The number of contributors was speified too low!",title="Wrong setting",icon="error")
        return(0) #return with error
      }
      if(nHp==0) nC_hp = NULL #SET TO NULL TO INDICATE THAT Hp is not to be calculated
      
      #get input to list: note: "fit_hp" and "fit_hd" are list-object from fitted model
      optFreq = get("optFreq",envir=mmTK)
      fst = optFreq$fst #obtain selected fst
      normalize = as.logical(optFreq$normalize)
      minF = optFreq$minF
      hyp = list(hp = list(NOC=nC_hp,condOrder=condOrder_hp,knownRef=knownref_hp,fst=fst,normalize=normalize,minF=minF),
                 hd = list(NOC=nC_hd,condOrder=condOrder_hd,knownRef=knownref_hd,fst=fst,normalize=normalize,minF=minF))

      #Extract selected model options
      model = list(MOD=gWidgets2::svalue(tabmodelB[2,2]), #distribution
                   DEG=booltxt[[gWidgets2::svalue(tabmodelB[1,2])]], #degradation
                   EXT=booltxt[[2]]) #svalue(tabmodelB[3,2])]]) #Extra model fit?

      #Store to a listed object
      set <- list(samples=mixData2,refData=refData2,model=model,hyp=hyp)     
      calcList = get("calcList",envir=mmTK) #obtain already calculated objects
      if(is.null(calcList)) calcList = list()
      calcList[[length(calcList)+1]] = set #insert object
      assign("calcList",calcList,envir=mmTK) #store data to envir
      
      return(1) #success
    } #end store settings from GUI to environment
    
  
    #type={EVID",DC"}
    gWidgets2::visible(mainwin) <- FALSE
    # dispose(tabmodel) 
    tabmodeltmp <- gWidgets2::glayout(spacing=spc,container= tabmodel[1,1] <- gWidgets2::ggroup(container=tabmodel) ) 
    #tabmodelCC = gWidgets2::glayout(spacing=10,container=(tabmodeltmp[1,1] <-gWidgets2::gframe(spacing=10,container=tabmodeltmp)))  
    tabmodelA = gWidgets2::glayout(spacing=5,container=(tabmodeltmp[1,1] <-gWidgets2::gframe("Hypothesis specification",container=tabmodeltmp))) 
    tabmodelB = gWidgets2::glayout(spacing=1,container=(tabmodeltmp[2,1] <-gWidgets2::gframe("Model specification",container=tabmodeltmp))) 
    tabmodelC = gWidgets2::glayout(spacing=1,container=(tabmodeltmp[3,1] <-gWidgets2::gframe("Further:",container=tabmodeltmp)))  
    
    #Obtain evid settings
    edwith = 6 #edit width
    NOCsel = nRefs + 1
	
    #Hypothesis selection: subframe of A
    Krange <- 1:maxKsetup #default Contr range
    txt <- "Contributor(s) under H"
    tabmodelA1 = gWidgets2::glayout(spacing=0,container=(tabmodelA[1,1] <-gWidgets2::gframe("Number of contributors",container=tabmodelA))) 
    tabmodelA1[1,1] <- gWidgets2::glabel("NOC:",container=tabmodelA1)
    tabmodelA1[1,2] <- gWidgets2::gcombobox(items=Krange,selected=NOCsel,editable=TRUE,container=tabmodelA1)
    gWidgets2::size(tabmodelA1[1,2]) = 3 #set with
    

#    tabmodelTmp= tabmodelA2$hp
    .createHypSetup = function(hypsel,checked=TRUE) { #Helpfunction
      tabmodelTmp = tabmodelA2[[hypsel]]
      for(rsel in refSel) { #indicate all selected references
        rowind = which(rsel==refSel)
        tabmodelTmp[rowind,1]  <- gWidgets2::glabel(paste0("C",rowind),container=tabmodelTmp)
        tabmodelTmp[rowind,2]  <- gWidgets2::gcheckbox(rsel,container=tabmodelTmp,checked=checked)
      }
    }
    
    tabmodelA2 = list()
    tabmodelA2$hp = gWidgets2::glayout(spacing=0,container=(tabmodelA[2,1] <-gWidgets2::gframe(paste0(txt,"p:"),container=tabmodelA))) 
    .createHypSetup("hp",TRUE)

    
    tabmodelA2$hd = gWidgets2::glayout(spacing=0,container=(tabmodelA[3,1] <-gWidgets2::gframe( paste0(txt,"d:"),container=tabmodelA)))
    .createHypSetup("hd",TRUE)
    
    #Model specification (replicates): 
    tabmodelB[1,1] <- gWidgets2::glabel(text="Degradation",container=tabmodelB)
    tabmodelB[1,2] <- gWidgets2::gradio(names(booltxt),selected = 2,horizontal = TRUE,container=tabmodelB)
    tabmodelB[2,1] <- gWidgets2::glabel(text="Model",container=tabmodelB)
    tabmodelB[2,2] <- gWidgets2::gradio(c("GA","NB"),selected = 1,horizontal = TRUE,container=tabmodelB)
    #tabmodelB[3,1] <- gWidgets2::glabel(text="Extended",container=tabmodelB)
    #tabmodelB[3,2] <- gWidgets2::gradio(names(booltxt),selected = 2,horizontal = TRUE,container=tabmodelB)

    #BLOCK FOR DEACTIVATING DEGRADATION OPTION
    allowDEG = TRUE
    optKit = get("optKit",envir=mmTK) #obtain kit info
    if(is.null(optKit$kit) || optKit$kit=="NONE") allowDEG = FALSE #Turn off DEG option if kit not defined
    if(allowDEG && any(!grepl(CEsep,allAlleles))) allowDEG = FALSE #Turn off DEG option if any alleles is missing CE
    if(!allowDEG) gWidgets2::enabled(tabmodelB[1,2]) = FALSE #Turn off DEG option
    
    
   #Calculation button:  
    tabmodelC[1,1] = gWidgets2::gbutton(text="CALCULATE",container=tabmodelC,handler=
    function(h,...) {
      ret = storeSettings() #store settings
      if(ret==0) return() #don't continue if value is 0
      doCalculate() #refresh MLE fit tab (i.e. it fits the specified model)
      gWidgets2::svalue(nb) <- 3 #go to mle-fit window (for all cases) when finished
    }) #end cont. calculation button

  
   gWidgets2::visible(mainwin) <- TRUE
   gWidgets2::focus(mainwin) <- TRUE
} #end refresh setup tab-frame

 
 ########################################################################################
 ############################  CALCULATION STEP  ########################################
 ########################################################################################

 #BEGIN MAIN FUNCTION (wrapper:  
 doCalculate = function(resID=NULL) {  #resID is index to show (elemnent in calcList)
   dec <- 4 #number of significant numbers to have in MLE print
   gWidgets2::visible(mainwin) <- FALSE
   
   calcList <- get("calcList",envir=mmTK)  #obtain list with calculations:    
   nCalcs = length(calcList)  #number of calculations done
   if(is.null(resID)) resID = nCalcs #uselast element if not selected
   
   set = calcList[[resID]] #use selected element
   if(is.null(set)) return() #NO DATA ARE GIVEN (return...)

   #PREPARE CALCULATIONS: DONE IF mlefit_hd is not obtained#
   #take out relevant parameters from stored list
   mod <- set$model #obtain param settings of model
   hyp <- set$hyp      #obtain hypotheses
   par <- set$param #obtain parameters
   calcHp = !is.null( hyp$hp$NOC) #whether Hp is considered (NOT NULL)
   
   calcMLE = function(hyp) { #helpfunction for optimizing likelihood function (hypothesis setup is changed)
     optMLE = get("optMLE",envir=mmTK) #Obtain optimizing options
     optFreq = get("optFreq",envir=mmTK)
     optKit = get("optKit",envir=mmTK) 
     popFreq = get("popFreq",envir=mmTK) 
     calibration = get("calibrated",envir=mmTK) 
     
     seed0 = optMLE$seed
     if(seed0==0) seed0=NULL #convert seed to NULL (0 means none)
     kit0 = optKit$kit
     if(kit0=="NONE") kit0 = NULL
     #print(paste0(optMLE$nDone," random startpoints with variation ",optMLE$delta," are applied in the optimizer.")) 
     
     MPSproto::inferEvidence(set$samples, popFreq, set$refData, hyp, calibration,
                             kit = kit0, optKit$platform, nDone = optMLE$nDone, delta = optMLE$delta,
                             steptol = optMLE$steptol, seed = seed0, verbose = FALSE, model =  mod$MOD)
   } 
   
   #fit under hp: (only for evidence)
   mlefit_hp <- NULL #not used otherwise
   if( calcHp  ) { #considering HP
     #nUhp <- mod$nC_hp-sum(mod$condOrder_hp>0) #number of unknowns
     print("Calculating under Hp...")
     time <- system.time({ mlefit_hp <- calcMLE(hy=hyp$hp) })[3]
     print(paste0("Optimizing under Hp took ",format(time,digits=5),"s"))
     if(!is.null(set$mlefit_hp) && set$mlefit_hp$fit$loglik>mlefit_hp$fit$loglik )  mlefit_hp <- set$mlefit_hp #the old model was better
   }
   
   #fit under hd: (does it for all methods)
   print("Calculating under Hd...")
   time <- system.time({    mlefit_hd <- calcMLE(hy=hyp$hd) })[3]
   print(paste0("Optimizing under Hd took ",format(time,digits=5),"s"))
   if(!is.null(set$mlefit_hd) && set$mlefit_hd$fit$loglik>mlefit_hd$fit$loglik )  mlefit_hd <- set$mlefit_hd #the old model was better
   
   fixmsg <- "The specified model could not explain the data.\nPlease re-specify the model."
   if(is.infinite(mlefit_hd$fit$loglik)) {
     gWidgets2::gmessage(fixmsg,title="Wrong model specification (Hd)",icon="error")
   } else if(calcHp && !is.infinite(set$mlefit_hd$fit$loglik) && is.infinite(mlefit_hp$fit$loglik)) {
     gWidgets2::gmessage(fixmsg,title="Wrong model specification (Hp)",icon="error")
   }
   
   #store MLE result: also store best mle-values once again (possible with re-running optim)
   set$mlefit_hp=mlefit_hp #store fitted mle-fit
   set$mlefit_hd=mlefit_hd #store fitted mle-fit
   
   #Store by-products of calculation:
   set = .getMeta(set)
   
   #Update calculated object to environment:
   calcList[[resID]]  = set #insert
   assign("calcList",calcList,envir=mmTK)  #obtain list with calculations:    
   
   refreshResults(resID)
 } #ENDING CALCULATIONS
 
 #helpfunction to get different "meta values" after calculations
 .getMeta = function(set) {
   #set = get("calcList",envir=mmTK)[[1]]
   fithp = set$mlefit_hp
   fithd = set$mlefit_hd
   calcHp = !is.null(fithp) #check if calculated under hp
   
   #Get log-lik values
   loglikHp = fithp$fit$loglik
   loglikHd = fithd$fit$loglik
   loglikiHp = NULL
   if(calcHp) loglikiHp = MPSproto::logLiki(fithp)
   loglikiHd = MPSproto::logLiki(fithd)
   
   #Obtain LR-values:
   set$log10LR <- (loglikHp - loglikHd)/log(10) #obtain LR (log10)
   
   set$log10LRtxt = ""
   if(length(set$log10LR)>0) set$log10LRtxt <- round(set$log10LR,2)
   set$upperLR <- MPSproto::getUpperLR(fithd,scale=TRUE) #obtain upper LR
   set$upperLRtxt = ""
   if(length(set$upperLR)>0 && set$upperLR!="") set$upperLRtxt = round(set$upperLR,2) #needs to be a numeric!
   set$log10LRi = (loglikiHp - loglikiHd)/log(10)
   
   #Get params
   thetaHp = fithp$fit$par
   thetaHd = fithd$fit$par
   nparHp = length(fithp$fit$phihat) #number of unknown params
   nparHd = length(fithd$fit$phihat) #number of unknown params
   
   #get adjusted loglik ("aic" alike)
   set$aic_hp <- loglikHp + nparHp
   set$aic_hd <- loglikHd + nparHd
   
   #Obtain text
   evidNames = names(set$samples) #extract reference names
   refNames = names(set$refData) #extract reference names
   condHp <- set$hyp$hp$condOrder
   condHd <- set$hyp$hd$condOrder
   iscondHp <- which(condHp>0)
   iscondHd <- which(condHd>0)
   POIind = setdiff(iscondHp,iscondHd) #get contribution indices under Hp but not hd
   CONDind = intersect(iscondHp,iscondHd) #get contribution indices under both Hp and Hd
   
   #INSERT
   set$evidNames = paste0(evidNames,collapse="/")
   set$POI = paste0(refNames[POIind],collapse="/")
   set$COND = paste0(refNames[CONDind],collapse="/")
   set$NOC = length(thetaHd$mx) #number of contributors
   set$POI_mx = thetaHp$mx[POIind] #get Mx of POI (under hp)
   set$POImxtxt = ""
   if(length(set$POI_mx)>0)  set$POImxtxt = paste0(signif(set$POI_mx*100,2),"%")
   set$COND_mx = thetaHp$mx[CONDind] #get Mx of COND (under hp)
   set$MOD = set$model$MOD
   set$DEG = ifelse(set$model$DEG,"Yes","No")
   set$EXT = ifelse(set$model$EXT,"Yes","No")
   #set$model #other info already inside
   
   #obtain hyp text:
   nUhp = set$NOC - length(iscondHp)
   nUhd = set$NOC - length(iscondHd)
   
   #Naming Components of Mx estimates:
   set$MxRefsHp = thetaHp$mx
   set$MxRefsHd = thetaHd$mx
   names(set$MxRefsHp)[condHp] = refNames[iscondHp]
   names(set$MxRefsHd)[condHd] = refNames[iscondHd]
   if(nUhp>0) names(set$MxRefsHp)[which(is.na(names(set$MxRefsHp)))] = paste0("Unknown ",1:nUhp)
   if(nUhd>0) names(set$MxRefsHd)[which(is.na(names(set$MxRefsHd)))] = paste0("Unknown ",1:nUhd)
   
   getHypTxt = function(ref,nU) {
     txt=""
     if(length(ref)>0) txt = paste0(txt, paste0(ref,collapse="/"))
     if(nU>0) {
       tmp = paste0(nU," unknown")
       if(nU>1) tmp = paste0(tmp,"s") #plural
       if(length(ref)>0) tmp = paste0(" + ",tmp) #additional to refs
       txt = paste0(txt, tmp)
     }
     return(txt)
   }
   set$hypTxt =  c(NA,getHypTxt(refNames[iscondHd],nUhd)) #create hyp-text
   if( calcHp  ) set$hypTxt[1] = getHypTxt(refNames[iscondHp],nUhp)
   names(set$hypTxt) = paste0("H",c("p","d"),": ",set$hypTxt) #insert
   
   return(set) #return updated object
 } #end get meta
 
########################################################################################
############################Tab 3: RESULT TABLE:########################################
########################################################################################
#This section involved calculations and obtaining different results
 
   #Helpfunction to show table in GDF (editable cells)
   .showGDFtable = function(title,table) {
     setwin2 <- gWidgets2::gwindow( title ,visible=FALSE) 
     guitab <- gWidgets2::gdf(items=table,container = setwin2) 
     gWidgets2::visible(setwin2) <- TRUE
     gWidgets2::focus(setwin2) <- TRUE
   }
 
  #helpfunction ran when call deconvolution
  .doDC <- function(mleobj) {
     dcopt <- get("optDC",envir=mmTK) #options when Deconvolution
     dcobj <- MPSproto::deconvolve(mlefit=mleobj,alpha=dcopt$alphaprob,maxlist=dcopt$maxlist) 
     DCtable1 <- .addRownameTable(dcobj$table2)
     colnames(DCtable1)[1] <- "Locus"
     DCtable2<-dcobj$table1
     DCtable3<-dcobj$table3
     DCtable4<-dcobj$table4
     assign("resDC",list(DCtable1,DCtable2,DCtable3,DCtable4),envir=mmTK) #assign deconvolved result to environment
     .refreshTabDC() #update table with deconvolved results
     gWidgets2::svalue(nb) <- 4 #go to deconvolution results window (for all cases) when finished     
     gWidgets2::focus(mainwin) <- TRUE
   }

  
  .f_showPerMarkerLR = function(h,...) {
    LRi = h$action
    print("------------LR-perMarker-----------------")
    tab = cbind(LRi,log10(LRi))
    colnames(tab) = c("LR","log10LR")
    print(tab)
  }
  
  #helpfunction to get ID-element of tabMLEresGUI/calcList 
  .getID = function() { 
    sel = gWidgets2::svalue(tabMLEresGUI) #get selected row
    if(is.null(sel)) return(0) #ID=0 means none
    return(as.integer(gsub("#","",sel)))
  }
  
  #Helpf 
  .f_exporttable = function(h,...) {
    resTable = tabMLEresGUI[]
    if(is.null(resTable)) return()
    .saveTable(resTable) #save table
  }

  .f_createreport = function(h,...) {
    print("Creating report")

    
  }

  #helpfunction to calculate extended
  .f_calcext = function(h,...) {
    
  }
  
  #helpfunction to get different options (a new window)
  .f_getoptions = function(h,...) {
    id = .getID()
    if(id==0) return()
    calcRes <- get("calcList",envir=mmTK)[[id]] #obtain table
    
    optwin <- gWidgets2::gwindow(paste0("More options for Hyp. set #",id) ,visible=FALSE)
    optgrp = gWidgets2::ggroup(horizontal = FALSE,spacing=5,container=optwin) 
    opthead = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=optgrp) #header with text
    opthead[1,1] = gWidgets2::glabel("Sample(s):",container=opthead)
    opthead[1,2] = gWidgets2::glabel(calcRes$evidNames,container=opthead)
    opthead[2,1] = gWidgets2::glabel("Hp:",container=opthead)
    opthead[2,2] = gWidgets2::glabel(calcRes$hypTxt[1],container=opthead)
    opthead[3,1] = gWidgets2::glabel("Hd:",container=opthead)
    opthead[3,2] = gWidgets2::glabel(calcRes$hypTxt[2],container=opthead)
    opthead[4,1] = gWidgets2::glabel("Upper log10LR:",container=opthead)
    opthead[4,2] = gWidgets2::glabel(calcRes$upperLRtxt,container=opthead)
    
    optgrid = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=optgrp) 
    optgrid[1,1] <- gWidgets2::gbutton("Show LR-per marker",container=optgrid,handler=function(h,...) {
      s0 = 2
      LRi = calcRes$log10LRi
      tab <- cbind(round(LRi,s0), signif(10^(LRi),s0+1))  #obtain LR per marker
      colnames(tab) = c("log10LR","LR")
      print(tab) #print to console
      .showGDFtable("LR per marker", .addRownameTable(tab,"Marker"))
      barplot(LRi,xlab="",ylab="log10 LR",main="",las=2) #also show barplot
    })
    
    optgrid[2,1] <- gWidgets2::gbutton("Show estimated params.",container=optgrid,handler=function(h,...) {
      hypsettxt = paste0(names(calcRes$hypTxt),collapse="/")
      #Create param table
      parHp = calcRes$mlefit_hp$fit$par
      parHd = calcRes$mlefit_hd$fit$par
      MxHp = parHp$mx
      MxHd = parHd$mx
      
      tab = round(cbind(unlist(parHp),unlist(parHd)),3)
      colnames(tab) = c("Hp","Hd")
      .showGDFtable(hypsettxt,.addRownameTable(tab,colname = "Param"))
      print(tab) #to console
    
      #Show mixture proportions in a separe plot:
      if(require(plotrix)) {
        MxHp2 = sort(calcRes$MxRefsHp,decreasing = TRUE) #sort to make same colors
        MxHd2 = sort(calcRes$MxRefsHd,decreasing = TRUE) #sort to make same colors
        s0 = 2
        labelHp = paste0(names(MxHp2)," ",round(MxHp2*100),"%") #Obtain labels
        labelHd = paste0(names(MxHd2)," ",round(MxHd2*100),"%") #Obtain labels
        plotfn = paste0("WOEpie",id) #file name
        graphics.off() #remove other plots first
        png(plotfn,height=400,width=1000)
        par(mfrow=c(1,2))
        marg = c(0,5,0,6)
        plotrix::pie3D(MxHp,radius=0.9,labels=labelHp,explode=0.1,main="Hp",mar=marg)
        plotrix::pie3D(MxHd,radius=0.9,labels=labelHd,explode=0.1,main="Hd",mar=marg)
        mtext(calcRes$evidNames,cex=1.5,outer=T,line=-2)
        dev.off()
        
        TopWin <- gWidgets2::gwindow(paste0("Mixture proportions for ",hypsettxt))#,width=800,height=400)
        gWidgets2::gimage(plotfn,container=TopWin)
        gWidgets2::focus(TopWin) = TRUE
        file.remove(plotfn) #remove file after
      }
      
    })
    
    optgrid[3,1] <- gWidgets2::gbutton("Show model validation",container=optgrid,handler=function(h,...) {
      alpha0 = 0.01
      
      hypsettxt = paste0("Hyp. set #",id)
      #hyp = gWidgets2::gconfirm("Select hypothesis (Yes=Hp, No=Hd)","Choose hypothesis",icon = "question")
      plotfn1 =  paste0("WOEvalidHp",id) #file name
      plotfn2 =  paste0("WOEvalidHd",id) #file name
      graphics.off() #remove other plots first
      size0 = 600
      
      png(plotfn1,height=size0,width=size0)
      MPSproto::validMLE(calcRes$mlefit_hp,"Hp",alpha=alpha0,createplot = TRUE, verbose = TRUE)
      dev.off() 
      
      png(plotfn2,height=size0,width=size0)
      MPSproto::validMLE(calcRes$mlefit_hd,"Hd",alpha=alpha0,createplot = TRUE, verbose = TRUE)
      dev.off() 

      TopWin <- gWidgets2::gwindow(paste0("Model validation for ",hypsettxt))#,width=800,height=400)
      ggrp <- gWidgets2::ggroup(container=TopWin)
      gWidgets2::gimage(plotfn1,container=ggrp)
      gWidgets2::gimage(plotfn2,container=ggrp)
      gWidgets2::focus(TopWin) = TRUE
      file.remove(plotfn1) #remove file after
      file.remove(plotfn2) #remove file after
    })
    
    optgrid[4,1] <- gWidgets2::gbutton("Show model fit",container=optgrid,handler=function(h,...) {
      hyp = gWidgets2::gconfirm("Select hypothesis (Yes=Hp, No=Hd)","Choose hypothesis",icon = "question")
      if(hyp) {
        MPSproto::plotTopMPS( calcRes$mlefit_hp )
      } else {
        MPSproto::plotTopMPS( calcRes$mlefit_hd )
      }
    })
    optgrid[5,1] <- gWidgets2::gbutton("Deconvolve",container=optgrid,handler=function(h,...) {
      hyp = gWidgets2::gconfirm("Select hypothesis (Yes=Hp, No=Hd)","Choose hypothesis",icon = "question")
      if(hyp) {
        resDC = MPSproto::deconvolve( calcRes$mlefit_hp )
      } else {
        resDC = MPSproto::deconvolve( calcRes$mlefit_hd )
      }
      print("Deconvolve finished")
      assign("resDC",resDC,envir=mmTK) #get deconvolved results
      .refreshTabDC()
      gWidgets2::svalue(nb) <- 4 #go to DC result tab when done
      .getFocus()
    })
    gWidgets2::visible(optwin) <- TRUE
  } #END OPTION WINDOW
  
  #BEGIN GUI HERE:
  tabMLEgrid = gWidgets2::glayout(horizontal = FALSE,spacing=5,container=tabMLE) #kit and population selecter
  tabMLEgrid[1,1] <- gWidgets2::gbutton(text="Export table",container=tabMLEgrid,handler=.f_exporttable)#
  tabMLEgrid[1,2] <- gWidgets2::gbutton(text="Create report",container=tabMLEgrid,handler=.f_createreport)#  
  tabMLEgrid[1,3] <- gWidgets2::gbutton(text="Infer LSAE (EXT)",container=tabMLEgrid,handler=.f_calcext)#  
  tabMLEgrid[1,4] <- gWidgets2::gbutton(text="More options",container=tabMLEgrid,handler=.f_getoptions)#  
  
  #This panel is for selecting different calculated results and showing comparison of calculated models
  tabMLEres <- gWidgets2::ggroup(container=tabMLE,expand=TRUE,fill=TRUE)
  tabMLEresGUI <- gWidgets2::gtable(items="",multiple = FALSE,container=tabMLEres, handler=.f_getoptions)
  gWidgets2::add(tabMLEres,tabMLEresGUI,expand=TRUE,fill=TRUE)#add to frame

  
  refreshResults = function(selrow=NULL) {  #update result table
    calcList <- get("calcList",envir=mmTK) #obtain table
    nres = length(calcList) #number of results
    if(nres==0) return() #no results
    
    #Curate table:
    cns = c("Sample(s)","PoI","Cond","NOC","MOD","DEG","EXT","LR","Mx","AIC")
    #x = calcList[[1]]
    resTable = t(sapply(calcList,function(x) c(x$evidNames,x$POI,x$COND,x$NOC,x$MOD,x$DEG,x$EXT,x$log10LRtxt, x$POImxtxt ,round(x$aic_hd,2))))
    colnames(resTable) = cns
    resTable = cbind(ID=paste0("#",1:nres),resTable)
    
    tabMLEresGUI[] <-  resTable  
    gWidgets2::size(tabMLEresGUI) <- list(column.widths=c(30,250,100,220,40,40,40,40,60,60,80)) 
    
    if(!is.null(selrow)) gWidgets2::svalue(tabMLEresGUI) <- selrow #highlight row
    .getFocus()
  }
  refreshResults() #Show already calculted evidence-results when program starts

##############################################################
###############Tab 5: Deconvolution results:##################
##############################################################

 .f_savetableDC = function(h,...) {
   if(is.null(DCtable)) {
    gWidgets2::gmessage("There is no deconvolution results available!")
   } else {
    .saveTable(DCtable[], "txt") #save deconvolution results
   }
 }
 .refreshTabDC = function(dctype=1) { #1=table1 (top marginal results),2=table2 (joint results), 3=table3 (all marginal results per genotype), 4=table4 (all marginal results per alleles)
   DCtables <- get("resDC",envir=mmTK) #get deconvolved results
   if(!is.null(DCtables)) {
     DCtable[] = .NAtoSign(DCtables[[dctype]])  #update Table
    }
 }

 #CREATE DECONV-GUI 
 tabDCa = gWidgets2::glayout(spacing=1,container=tabDC) #table layout
 tabDCb = gWidgets2::ggroup(spacing=1,container=tabDC,expand=T,fill=T)
 itemvecDC = c("Top Marginal","All Joint","All Marginal (G)","All Marginal (A)")
 tabDCa[1,1] <- gWidgets2::glabel("Select layout:",container=tabDCa)
 tabDCa[1,2] <-  gWidgets2::gradio(items=itemvecDC,selected=1,horizontal=TRUE,container=tabDCa,handler=function(x) {
   .refreshTabDC( which(itemvecDC==gWidgets2::svalue(tabDCa[1,2])) )
 })
 tabDCa[2,1] <- gWidgets2::gbutton(text="Save table",container=tabDCa,handler=.f_savetableDC)  
 
 #ADD DC TABLE
 DCtable = gWidgets2::gtable(items="",multiple = TRUE,container = tabDCb,expand=T,fill=T)
 gWidgets2::add(tabDCb,DCtable,expand=T,fill=T)
 .refreshTabDC() #open results when program starts

 for(nbvisit in 4:1) gWidgets2::svalue(nb) <- nbvisit #visit tables 
 .getFocus()

} #end funcions
