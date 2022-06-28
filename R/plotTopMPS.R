#' @title plotTopMPS
#' @author Oyvind Bleka
#' @description MPS data visualizer (interactive)
#' @details Plots the expected peak heights of the top genotypes. The peak heights for corresponding alleles (one sample) are superimposed.
#' 
#' @param MLEobj Fitted object using inferEvidence
#' @param DCobj An object returned from devonvolve: Must be run with same object as mlefit
#' @param locYmax A boolean of whether Y-axis should be same for all markers (FALSE) or not (TRUE this is default)
#' @param options A list of possible plot configurations. See comments below
#' @param withStutterModel Whether taking the stutter model into account
#' @param returnOnly Whether only returning fig (not plotted)
#' @return sub A plotly widget
#' @export

#DCobj=NULL;grpsymbol=":";locYmax=TRUE;options=NULL;withStutterModel=TRUE
plotTopMPS = function(MLEobj,DCobj=NULL,locYmax=TRUE,options=NULL,withStutterModel=TRUE, returnOnly=FALSE) {
 MPSsymbol = ":" #Used for extracting CE for MPS strings. Example is "10:[ATCG]10".
 Qallele = "99"
 if(is.null(options$h0)) { h0 = 300 } else { h0 = options$h0 } # 5500/nrows #standard height for each dye (depends on number of rows? No)
 if(is.null(options$w0)) { w0 = 1800 } else { w0 = options$w0 } # standard witdh when printing plot
 if(is.null(options$marg0)) { marg0 = 0.015 } else { marg0 = options$marg0 } #Margin between subplots
 if(is.null(options$txtsize0)) { txtsize0 = 12 } else { txtsize0 = options$txtsize0 } #text size for alleles
 if(is.null(options$locsize0)) { locsize0 = 20 } else { locsize0 = options$locsize0 } #text size for loci
 if(is.null(options$minY)) { minY = 100 } else { minY = options$minY } #default minimum Y-axis length
 if(is.null(options$ymaxscale)) { ymaxscale = 1.06 } else { ymaxscale = options$ymaxscale } #y-axis scaling to the locus name positions
 if(is.null(options$grptype)) { grptype="group" } else { grptype = options$grptype }#,"stack" "group" is default 

 dat = MLEobj$data #obtain data
 
 sn = names(dat[[1]]$samples) #get samples names
 nS = length(sn) #number of replicates
 locs = names(dat) #get locus names
 nL = length(locs)

 refn = names(dat[[1]]$refs)
 refn = refn[MLEobj$hypothesis$cond>0] #ref names conditioned on
 nrefs = length(refn)

 if(is.null(DCobj)) DCobj <- deconvolve(MLEobj,maxlist=1) #get top candidate profiles
 
 #extract info from DC (deconvolution) object
 topG <- sapply(DCobj$rankGi,function(x) x[1,-ncol(x),drop=F])
 if(is.null(nrow(topG))) topG <- t(topG) #consider as matrix
 pG <- as.numeric(sapply(DCobj$rankGi,function(x) x[1,ncol(x)])) #probabilities
 names(pG) = toupper(names(DCobj$rankGi)) #assign loci names

 #estimates:
 nC <- MLEobj$hypothesis$NOC #number of contributors
 par <- MLEobj$fit$par #get estimates parameters
 mx <- par$mx
 mu <- par$mu
 #omega <- par$omega

 theta_Am <- MLEobj$prepareC$markerEfficiency #Object already stored in mlefit. returned from prepareC function
 names(theta_Am) = MLEobj$prepareC$locNames

 #Degradation model 
 beta <- 1
 kitinfo = NULL
 hasKit = !is.null(MLEobj$kit) #is kit defined? THen grpsymbol is included!
 if(hasKit) {
  beta <- par$beta #extract degrad param
  kitinfo = getMPSkit(MLEobj$kit) #names(kitinfo)
 }
 
 #Stutter model: Need to structure data
 if(withStutterModel) stuttList = getStutterExpectations(MLEobj$prepareC)

#  Ccols <- c("cyan","green","coral","gold","hotpink","darkorange","lightgoldenrod") #c("black","gray","brown","darkorange") #contributor cols
 Ccols <- c("blue","green","coral","gold","hotpink","darkorange","lightgoldenrod") #c("black","gray","brown","darkorange") #contributor cols


 df = NULL#store data: (sample,marker,allele,height)
 for(ss in sn) { #create a seperate EPG plot for each samples
# ss=sn[1]
  for(loc in locs) {
# loc=locs[1]
    muLoc = mu*theta_Am[loc] #use estimated marker efficency here
	    
    edat = dat[[loc]]$samples[[ss]] #get evid data   
    if(nrefs==0) {
      rdat = NULL
    } else {
      rdat = dat[[loc]]$refs[refn] #get 
    }

    if(is.null(edat) && is.null(rdat)  ) next #skip if no data (evid or ref)

    hv = as.numeric(edat) #coverage
    av <- names(edat) #alleles
    #Check MPSsymbol is used in front
    if(length(av)>0 && all(grepl(MPSsymbol,av))) {
      tmp = strsplit(av,MPSsymbol) 
      avCE = sapply(tmp,function(x) x[1]) #extract RU allele
      avSEQ = sapply(tmp,function(x) paste0(x[-1],collapse=MPSsymbol)) #collapse other levels if several
    } else {
      avSEQ <- avCE <- av #sequence as before
    }
    av = avSEQ
    av2 = unique(unlist(rdat)) #alleles of referene
    G = unlist(topG[,colnames(topG)==toupper(loc)]) #get genotypes
    av2 = unique( c(av2,unlist(strsplit(G,"/"))) ) #get all alleles 
    adda = av2[!av2%in%av] #obtain missing allele

    #add missing:
    if(length(adda)>0) {
     av = c(av,adda)
     hv  =  c(hv, rep(0,length(adda)) ) #add zero height to
     avCE = c(avCE,adda)
     avSEQ = c(avSEQ,adda)
    }
 
    if(length(av)==0) { #add dummy variables if no alleles
      av <- ""
      hv <- 0 
      av1 <- av2 <- rep("",length(av))
    } else { #otherwise if observed:

     #sort alleles increasingly (handle strings)
      suppressWarnings({ av1n = as.numeric(avCE)}) #all numbers?
      if(any(is.na(av1n))) av1n = avCE #get back to string (order wrt sequence lengt)
      ord = order(av1n) 
      av = av[ord]
      hv = hv[ord]
      av1 = avCE[ord]
      av2 = avSEQ[ord]
    }

    #ref text under each allele (follow original av)
    reftxt <- rep("",length(av))
    if(nrefs>0) {
     for(rr in 1:nrefs) { #for each ref
      indadd = which(av%in%unlist(rdat[[rr]])) #index of alleles to add to text
      hasprevval = indadd[nchar(reftxt[indadd])>0] #indice to add backslash (sharing alleles)
      reftxt[ hasprevval ] = paste0(reftxt[ hasprevval ],"/")      
      reftxt[indadd] = paste0( reftxt[indadd], rr)
     }
    }

    if(hasKit) { #if kit info provided (for degradation)
     tmp = kitinfo[toupper(kitinfo$Marker)==loc,]
     
     ind = match(av1,tmp$Allele) #get index to extract
     bv = tmp$Size[ind]  #get sizes directly from lookup
     isna = which(is.na(ind)) #which alleles are missing?
     if(length(isna)>0) avuse = as.numeric(tmp$Allele) #alleles available (called only once)

     for(missind in isna) {#impute missing bp:
      if( av[missind]==Qallele) { #if it was the Qallele
       bv[missind] = max(tmp$Size) #Use maximum bp
      } else {
       newa = as.numeric(av1[missind ]) #get allele name
       bv[missind] = tmp$Size[which.min(abs(newa - avuse ))]  #use base pair of closest allele
      }
     }
    }

    #GET EXPECTation (Cumulative over all contributors)
    EYmat <- matrix(0,nrow=length(av),ncol=nC) #before and after stutter (cumulative for each contributor)
    for (aa in av) { # Loop over all alleles in locus
       ind = which(av==aa)
       contr <- sapply(strsplit(G,"/"),function(x) sum(x%in%aa)) #get contribution
       EY = cumsum(contr*mx)*muLoc #cumulative sum
       if(hasKit && beta!=1) {
         EY <- EY*beta^((bv[ind]-125)/100) #expected peak heights for each contributors
       }
       EYmat[ind,] <- EY 
    }

    if(withStutterModel) {
      stuttTab = stuttList[[loc]]
      if(!is.null(stuttTab)) {
        #scale expectations with stutter:
        nStutt = nrow(stuttTab)
        
        EYmat2 <- EYmat
        for(s in seq_len(nStutt) ) {
          fromInd = match(stuttTab[s,2],av)
          toInd = match(stuttTab[s,3],av)
          stuttProp = as.numeric(stuttTab[s,4])
          EYmat2[toInd,] = EYmat2[toInd,] + stuttProp*EYmat[fromInd,] #OBTAINED stutters
          EYmat2[fromInd,] = EYmat2[fromInd,] - stuttProp*EYmat[fromInd,] #SUBTRACTED stutters
        }
        EYmat <- EYmat2 #override with stutter products
      }
    }
    EYv <- apply(EYmat,1,function(x) paste(x,collapse="/"))

    df_new = data.frame(Sample=ss,Marker=loc,Allele=av,Height=hv,reftxt=reftxt,Allele1=av1, Allele2=av2,EXP=EYv,stringsAsFactors=FALSE)
    
    df = rbind(df,   df_new)
  } #end for each loci
 } #end for each samples
 #df[,-c(3,7)]
 #colnames(df)
 EYmax <- max(as.numeric(unlist( strsplit(df$EXP,"/") )))
 ymax1 <- ymaxscale*max(minY,EYmax,df$Height) #global max y

#GRAPHICAL SETUP BASED ON SELECTED KIT:
ncols = 5 #number of locs per row (depend on amax)?
maxA = max(aggregate(df$Height,by=list(df$Sample,df$Marker),length)$x)
if(maxA<=2 && nL>40) ncols=10 #Number of cols should depend on allele outcome (SNP vs STR)
if(!is.null(options$ncols)) ncols=options$ncols #use number of column specified in options 
h1 = h0*ncols #standard height for each graph (DEPEND ON NUMBER COLS)	
nrows0=ceiling((nL+1)/ncols) #number of rows to use (use 10 per column)
if(nrefs>0) nrows0=ceiling((nL+1)/ncols) #number of rows to use (use 10 per column)

hline <- function(y = 0, color = "black",xr=0:1) {
  list(
    type = "line", 
    x0 = xr[1], 
    x1 = xr[2], 
    y0 = y, 
    y1 = y, 
    line = list(color = color,dash = 'dot',width=2)
  )
}
  
transdeg = .8 #transparancy degree

#SEPARATE PLOTS

for(ss in sn) {
 #ss = sn[1]
 #locs = locs[1:30]
 plist = list() #create plot object for each marker
 for(loc in locs) {
#loc  = locs[21]
   AT0 <- MLEobj$calibration[[loc]]$AT #temporary on analytical threshold

   dfs = df[df$Sample==ss & df$Marker%in%loc,] #extract subset 
   dfs$Allele = as.character(dfs$Allele)
   dfs$Allele1 = as.character(dfs$Allele1)
   dfs$Allele2 = as.character(dfs$Allele2)
  
   nA = length(dfs$Allele)
   av1 = unique(dfs$Allele1) #get unique alleles
   nA1 = length(av1) #number of unique
   xpos = 0:(nA-1) #position of alleles (standard)
   reptab = table(dfs$Allele1)
   nR = max(reptab) #number of layers
   repcols = rep("black",nR) #gray.colors(nR, start = 0.3, end = 0.9) #color level will be adjust regarding #layers

   repcol = rep(NA,nA)#get layer index
   for(aa in av1) {
    ind = which(dfs$Allele1==aa)
    repcol[ind] = 1:length(ind)
   }
   xpos1 = rep(NA,nA)
   for(i in 1:nA) xpos1[i] = which(dfs$Allele1[i]==av1)-1
   xpos0 = 0:(nA1-1)
   atxtL = nchar(av1) #get allele length

   p = plotly::plot_ly(dfs,height=h1,showlegend = TRUE,colors=repcols[1:nR] )
   p = plotly::add_trace(p,type = "bar", x = xpos1,y=~Height,name=repcol,hoverinfo="y+text",hoverlabel=list(font=list(size=12),namelength=1000),text =~Allele,color=as.factor(repcol))

   #ADD EXPECTATIONS
   shifts = seq(-0.25,0.25,l=nR)
   if(nR==1) shifts = 0

   Ev = strsplit(dfs$EXP,"/") #extract
   shapeList = list()  #add opacity shapes
   cc = 1 #counter
   for(i in 1:length(Ev)) {
    x0 = xpos1[i] #get allele position (centered)
    x1 = x0+shifts[repcol[i]] #layer decides layer pos
    Ev2 = c(0,Ev[[i]])
    for(j in 1:length(Ev[[i]])) { #for each contributor 
     if(Ev2[j+1]==Ev2[j]) next #skip if equal
      shapeList[[cc]] = list(type = "rect",fillcolor = Ccols[j], line = list(color = Ccols[j],width=0.1), opacity = transdeg,x0 =x1-1/(2*(nR+1)), x1 = x1+1/(2*(nR+1)),y0 = Ev2[j], y1 = Ev2[j+1])#,xref="x", yref = "y")
      cc = cc + 1
    }
   }

   #Add threshold lines to shapes
   if(!is.null(AT0))  {
      shapeList[[cc]] = hline(AT0,xr=c(-0.5,nA1-0.5))
      cc = cc + 1
   }
   if(locYmax)  ymax1 = ymaxscale*max( na.omit( c(minY,AT0,as.numeric(unlist(Ev)),dfs$Height) ))  #get max 

   dye2 = "black" 
    prob = pG[names(pG)==toupper(loc)]
    if( length(prob)>0 ) {
     if( prob>=.95 ) {
      dye2 = "forestgreen"
     } else if(prob>=.9) {
      dye2 = "orange"
     } else {
      dye2 = "red"
     }
     p = plotly::add_annotations(p, x=(nA1-1)/2 ,y=ymax1,text=loc,showarrow=FALSE,font = list(color = dye2,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAME
    }
   if(max(atxtL)<=5) p = plotly::add_annotations(p, x=xpos0 ,y=rep(0,nA1),text=av1 ,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-10)  #ADD ALLELE NAMES
   #Rotate alleles if many alleles? Using tickangle: layout(xaxis=list(tickangle=-45))
   #p = plotly::add_annotations(p, x=(nA1-1)/2 ,y=ymax1,text=loc,showarrow=FALSE,font = list(color = 1,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAMES 

   if(nrefs>0) { #need to take care of different layers (shift on x-axis) + missing PHs
     refuse = which(dfs$reftxt!="")  #!duplicated(dfs$Allele1) #only put reference under relevant alleles
     for(rr in refuse) { #for each refs (some are missing PH)
      x0 = xpos1[rr] #get allele position (centered)
      x1 = x0+shifts[repcol[rr]] #layer decides layer pos
      p = plotly::add_annotations(p, x=x1,y=0,text=dfs$reftxt[rr],showarrow=FALSE,font = list(color = repcols[repcol[rr]],family = 'sans serif',size = txtsize0),yshift=-25)  #ADD ALLELE NAMES
      if(dfs$Height[rr]==0)  p = plotly::add_trace(p,type = "scatter",mode="markers", x = x1,y=0,name=repcol[rr],hoverinfo="y+text",hoverlabel=list(font=list(size=12),namelength=1000),text=dfs$Allele[rr],color=as.factor(repcol[rr]))  #add a point to missing alleles (to get hovering)
     } #end for each refuse
   } #end if references
   p = plotly::layout(p,xaxis = list(showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline = TRUE,title = ""),shapes=shapeList)
   plist[[loc]] <- p
 } #end for each loc
# subplot(plist, nrows = nrows0, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE)%>%plotly::layout(title=ss,barmode = grptype)%>%hide_legend()

 #In last plot we will show the Mx and the contributors 
 p = plotly::plot_ly(x = xpos,y=rep(0,length(xpos)),height=h1,showlegend = FALSE,type="scatter",mode="markers")
 p = plotly::add_annotations(p, x=tail(xpos,1),y= c(ymax1-ymax1/6*(1:length(mx))),text= paste0("Contr.C",1:length(mx),"(",Ccols[1:length(mx)],")=",signif(mx,3)),showarrow=FALSE,font = list(family = 'sans serif',size = 15),xshift=0,xanchor = 'right') 

 if(nrefs>0) p = plotly::add_annotations(p, x=0,y=c(ymax1-ymax1/6*(1:nrefs)-ymax1/12),text= paste0("Label ",1:nrefs,": ",refn),showarrow=FALSE,font = list(family = 'sans serif',size = 15),xshift=0,xanchor = 'left')  #ADD ALLELE NAMES
 p = plotly::layout(p,xaxis = list(showline=FALSE, showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline=FALSE,showticklabels = FALSE,title = ""))#,colorway =dye2) 
 plist[[nL+1]] = p

 sub = plotly::subplot(plist, nrows = nrows0, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE)
 sub = plotly::layout(sub, title=ss,barmode = grptype)
 sub = plotly::hide_legend(sub)
 sub = plotly::config(sub, scrollZoom=TRUE, displaylogo=FALSE,modeBarButtonsToRemove=c("lasso2d","select2d","hoverClosestCartesian","hoverCompareCartesian","toggleSpikelines"),toImageButtonOptions=list(width=w0))

 if(!returnOnly) print(sub)
 } #end for each samples
 if(returnOnly) return(sub) #return last created

} #end function
