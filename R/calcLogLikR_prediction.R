#' @title calcLogLikR_prediction
#' @description Likelihood calculation of (evidence) prediction model (uses calibration info) in R
#' @param par parameters A list(mx,mu,omega,beta), The ordinary scale (ttheta)
#' @param c data object (returned from prepareC_calibration)
#' @export 


calcLogLikR_prediction = function(par,c) {
  #keep_threshold A threshold used to filter which joint genotype iterations
  
  #MARKER EFFICIENCY (Am):
  theta_Am = c$markerEfficiency #fixed and known
  nJointCombs = c$nGenos^c$NOU #obtain number of combinations per markers
  startIndMarker_nJointCombs = cumsum(c(0,nJointCombs)) 
  nCombs = sum(nJointCombs) #total number of combinations

  #START FUNCTION (MIMICING C++ script)
  NOC = c$NOC #copy
  nLocs = c$nLocs #copy

  ######################
  #Structure parameters#
  ######################
  
  theta_mx = par$mx
  theta_mu = par$mu
  theta_omega = par$omega
  theta_beta = par$beta
  
  model = c$model #obtain model
  if(model=="GA") {
    scale0 = theta_mu*theta_omega^2 #obtain scale param (same for all obs)
    shape0 = 1/theta_omega^2 #shape for 'full het allele'
  } 
  if(model=="NB") {
    size0 = theta_mu/(theta_mu*theta_omega^2-1) #obtain size param (same for all obs)
    shape0 = theta_mu #this is now expectation
  } 
  
  
  #print(stutterTypes)
  lik_outcomeList = list() #storing all outcomes in the list (to be returned)
  for(m in seq_len(nLocs)) { #run through each locus
#  m=1
    m0 = m - 1 #starts from index 0
    #if(verbose) print(paste0(round(m/nLocs*100),"%"))
    nAlleles = c$nAlleles[m] #number of alleles
    numGenos1p = c$nGenos[m] #as.integer(nAlleles*(nAlleles+1)/2) #number of genocombs
    startIndMarker_nAlleles0 = c$startIndMarker_nAlleles[m]; 
    startIndMarker_nAllelesReps0 = c$startIndMarker_nAllelesReps[m]; 
    startIndMarker_nStutters0 = c$startIndMarker_nStutters[m];
    startIndMarker_outG1allele0 = c$startIndMarker_outG1allele[m];
    startIndMarker_outG1contr0 = c$startIndMarker_outG1contr[m];
    
    
    #Obtain locus specif params (AT,Am,dropin model):
    nRep0 = c$nSamples[m]
    NOK0 = c$NOK[m] #number of knowns
    NOU = NOC - NOK0 #number of unknowns
    AT0 = c$AT[m] #obtain analytical threshold
    nNoiseParam0 = c$nNoiseParam[m] #Dropin count param (geometric)
    fst0 = c$fst[m]
    
    #Create scale parameters
    shape1 = theta_Am[m]*shape0 #shape for 'full het allele'
    
    #prepare vectors for knowns vs unknowns    
    GindKnown <- kindKnown <- rep(0,NOK0)
    kindUnknown <- rep(0,NOU)
    cc <- jj <- 1
    for(kk in seq_len(NOC)) { #for each contributors:
        knownGind0 = c$knownGind[ NOC*m0 + kk ] #copy genotype index  (KNOWN)			
        if(knownGind0>=0) { #If contributor is known (genotype given)
          GindKnown[cc] = knownGind0 #copy genotype index
          kindKnown[cc] = kk ; #insert contributor index
          cc = cc + 1 #update counter for knowns
        } else { #if contributor is unknown (genotype not given)
          kindUnknown[jj] = kk; #insert contributor index
          jj = jj + 1 #update counter for unknowns
        }
    }
    
    nTyped = c$nTyped[m];			
    maTypedvec <- shapevK <- shapev0 <- rep(0,nAlleles)
    for (aa in seq_len(nAlleles)) { #traverse each observed alleles (indices)
      maTypedvec[aa] = c$maTyped[ startIndMarker_nAlleles0 + aa ]; #copy previously typed alleles
      shapev0[aa] = exp( log(shape1) + c$basepair[ startIndMarker_nAlleles0 + aa ]*log(theta_beta) ); #scaling with degradation model (assumed already scaled)
    }
    
    
    #Sum up contribution for each alleles (Taking into account mix proportions): ONLY CALCULATED FOR KNOWN CONTRIBUTORS INITIALLY
    for (kk in seq_len(NOK0)) { #for each known contributors 
      for (aa in seq_len(nAlleles)) { #traverse each alleles (indices)			
        contr = c$outG1contr[ startIndMarker_outG1contr0 + nAlleles*GindKnown[kk] + aa] * theta_mx[kindKnown[kk]]; #contr from contr k to allele aa
        shapevK[aa] = shapevK[aa] + contr
      }
    }		
    
    
    nGjoint = c$nJointCombs[m] #number of combinations to traverse
    lik_outcomeList[[m]] = rep(NA,nGjoint)
    for(iter1 in seq_len(nGjoint)) { #iterate through all samples (for each locus)
#      iter1=1
      iter = iter1 - 1       
      shapev = shapevK #copy
      maTypedvec2 = maTypedvec; #Creating copy of counter (per allele)
      nTyped2 = nTyped; #creating copy of total counter
      
      jointGind = rep(0,NOC); #index of contributors
      genoProd = 1 #init geno prob product
      inserted = FALSE; #boolean of whether all digits are inserted (Rest is zero padded)
      modrest = iter; #used to keep remained after modulo (init as iter number)
      for (k in seq_len(NOU)) { #for each unknown contributors (summing up wrt both contr (outG1contr) and mx (mixprop)): Need each contr to derive shapev
        if (!inserted) { #if not all digits inserted
          if ( k>1 ) {
            modrest = as.integer((modrest - jointGind[k - 1]) / numGenos1p); #extract remaining, divide to get to next digit (necessary with int converion?)
          }
          jointGind[k] = modrest %% numGenos1p; #INSERT NUMBER: convert number to "numGenos1p" basis 	
          if (modrest < numGenos1p) { #check if rest is smaller than base (only run if not inserted)
            inserted = TRUE;  #then all digits are inserted
          }	
        }
      
        #Calculate shape param
        for(a in seq_len(nAlleles)) {
          #shapev[a] = shapev[a] + outG1contr[jointGind[k]+1,a] * theta_mx[k]; #contr from contr k to allele a. NOTICE THE k+*NOK shift!
          contr = c$outG1contr[ c$startIndMarker_outG1contr[m] + nAlleles*jointGind[k] + a] * theta_mx[ kindUnknown[k] ]; #contr from contr k to allele a. NOTICE THE k+*NOK shift!
          shapev[a] = shapev[a] + contr
        }
        
        #CALCULATE GENOTYPE PROBS OF UNKNOWNS (MAY BE RELATED != -1)
        startInd_alleles = startIndMarker_outG1allele0 + 2*jointGind[k];
        for(a in seq_len(2)) { #loop through both alleles
          aind = c$outG1allele[ startInd_alleles + a] + 1  #get allele index of genotype g_a (add index of 1)
          genoProd =  genoProd*(fst0*maTypedvec2[aind] + (1-fst0)*c$freq[startIndMarker_nAlleles0 + aind]) / (1 + (nTyped2-1)*fst0); 
          maTypedvec2[aind] = maTypedvec2[aind] + 1; #update allele count for particular genotype
          nTyped2 = nTyped2 + 1; #update total count
        }
        if( c$outG1allele[ startInd_alleles + 1 ] != c$outG1allele[ startInd_alleles + 2 ] ) { #heterozygous variant
          genoProd = 2*genoProd; #multiply by 2  if het. variant
        } 
      } #end for each unknown
      #print(jointGind)
 
      #Scaling shape parameter with degrad model
      for (a in seq_len(nAlleles)) { #traverse each observed alleles (indices), having PH>0
        shapev[a] = shapev[a]*shapev0[a]; #scaling with degradation model (assumed already scaled)
      }
      
      #Create shape parameters
      shapev2 = shapev #make copy
      
      #loop through stutter relations to obtain modified shape vector
      for(stuttind in seq_len( c$nStutters[m])) { 
#        stuttind=1
        stuttind1 = c$startIndMarker_nStutters[m] + stuttind
        
        stuttFromInd = c$stuttFromInd[stuttind1] + 1 #NB: CAREFUL WITH INDEX
        
        if(shapev[stuttFromInd]>0) {
          stuttToInd = c$stuttToInd[stuttind1] + 1 #NB: CAREFUL WITH INDEX
          stuttExp = c$stuttExp[stuttind1]  
          
          shapev2[stuttToInd] = shapev2[stuttToInd] + stuttExp*shapev[stuttFromInd] #OBTAINED stutters
          shapev2[stuttFromInd] = shapev2[stuttFromInd] - stuttExp*shapev[stuttFromInd] #SUBTRACTED stutters
        }
      }

      #calculate evidence likelihood (traverse through each alllee)
      val = 0 #calculate likelihood of observed coverages
      nDropin = rep(0,nRep0) #count number of dropin
      for (a in seq_len(nAlleles)) { #traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
        for(r in seq_len(nRep0)) { #traverse each replicates (observed alleles indicated by PH)
          cind = startIndMarker_nAllelesReps0 + (a-1)*nRep0 + r; #get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
          peak = c$peakHeights[cind]; #obtain coverage/peak
          
          #IF NOT DROPOUT
          if( peak >= AT0 ){ #if PH>0 (this is same as PH>=AT since threshold has already been applied)
            if(shapev2[a] > 0) { #If contribution and PH>0  					//CONTRIBUTION SET (A)
              
              if(model=="GA") val = val + dgamma(peak,shapev2[a],scale=scale0,log=T) 
              if(model=="NB") val = val + dnbinom(as.integer(peak),size=size0,mu=shapev2[a],log=T) 
              
            } else { #If contribution and PH>0 	#DROPIN SET (C)	
              
              val = val + c$noiseSizeWeight[ cind ]; #likelihood for noise coverage
              nDropin[r] = nDropin[r] + 1; #count dropin for particular replicate															   
            }							
            
            #OTHERWISE IT IS DROPOUT
          } else if(shapev2[a] > 0) { #IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)
            
            if(model=="GA") val = val + pgamma(AT0,shapev2[a],scale=scale0,log.p = T); 
            if(model=="NB") val = val + pnbinom(as.integer(AT0-1), size=size0,mu=shapev2[a],log.p = T); 
            
          }
        }#end for each reps
      } #end for each alleles
      
      #Liklihood OF NUMBER OF NOISE DROP-INS: Traverse each replicates (observed alleles indicated by PH)
      for(r in seq_len(nRep0)) { #traverse each replicates
        val = val + log(nNoiseParam0) + nDropin[r]*log(1-nNoiseParam0); #likelihood for number of noise alleles (geometrical distribution)															   													   
      }

      #INCLUDE LIKELIHOOD
      lik_outcomeList[[m]][iter1] = exp(val)*genoProd #calculate P(E|g)P(g)
    } #end for each joint genotype iter
  }   #end for each locus
  names(lik_outcomeList) = c$locNames
  logLik_marker = log(sapply(lik_outcomeList,sum)) #summing across loci and samples
  logLik = sum(logLik_marker)
  
  return( list(logLik=logLik, logLik_marker=logLik_marker, lik_outcomeList=lik_outcomeList) )
} #end function

