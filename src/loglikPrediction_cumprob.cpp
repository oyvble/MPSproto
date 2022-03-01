//This script contains functions for calculating the likelihoods function.
//AUTHOR: Oyvind Bleka, November 2020
/*ABOUT:

- No structuring of data is needed (already done in prepareC)
- All markers are looped outside the inner large-sum loop.
- Possible to traverse certain genotype combinations (tremendous speedup for 3 and more contributors)

*/


#include <vector> //vector storage
#include <cmath> //includes lgamma
#include <thread> //used to obtain number of logical processes
#include <omp.h> //parallelization
#include <Rmath.h> //includes pgamma

using namespace std;

/*MAIN IMPORTS DATA AN STRUCTURES THE DATA FOR MARKER BASED FUNCTIONS*/
extern "C" {
		

//The P(E|theta)= sum_g P(E|g,theta)*P(g) expression
void loglikPrediction_cumprob(double *pvalVEC, double *maxY, int *nJointCombs, int *NOC, int *NOK, 
	double *mixprop, double *mu, double *sigma, double *beta, double *theta_Am, //eps is stutter proportion vector (1 element for each kind)
	double *AT, double *fst, double *nNoiseParam, double *noiseSizeWeightLong, double *noiseSizeWeightLong2, //note the additional weight 
   	int *nMarkers, int *nRep, int *nAllelesVEC, int *startIndMarker_nAlleles, int *startIndMarker_nAllelesReps,	
	double *peaksLong, double *freqsLong, double *nTypedLong, double *maTypedLong, double *basepairLong,
	int *nGenos, int *outG1allele, int *outG1contr, int *startIndMarker_outG1allele, int *startIndMarker_outG1contr,
	int *nStutters, int *stuttFromIndVEC, int *stuttToIndVEC, double *stuttExpVEC, int *startIndMarker_nStutters, 	
	int *knownGind) {
	
	
	int numThreads = thread::hardware_concurrency();
	int useThreads = numThreads; //number of threads to use
	omp_set_num_threads(useThreads);  //set number of threads to use 
	
	//parameter values (per marker):
	double mu0 = mu[0]; //assume constant
	double sigma0 = sigma[0]; //assume constant
	double beta0 = beta[0]; //assume constant

	//Prepare transformation of parameters:
	double theta_sigmasq = sigma0*sigma0;
	double scale0 = mu0*theta_sigmasq; //obtain scale param (same for all obs)
	double const1 = 1/scale0; //constant 1 (used in pgamma,dgamma)
	double const2 = log(scale0); //constant 2 (used in dgamma)
	double shape0 = 1/theta_sigmasq; // //theta_Am[locind]/theta_omegasq; //shape for 'full het allele'    

	//Calculating the loglik (over all markers), which is returned
	
	const double smalltol = 1.0e-30; //a tiny number > 0 (avoiding zero roundoff errors)
	for(int locind=0; locind< *nMarkers; locind++) {	//for each marker:
	
		//OBTAIN CONSTANTS FOR SPECIFIC MARKER:
	
		//default settings
		double fst0 = fst[locind]; //theta-correction param
		double AT0 = AT[locind]; //analytical threshold
	    double nNoiseParam0 = nNoiseParam[locind]; //Number of noise count param (geometric)
		double shape1 = shape0*theta_Am[locind]; //scale shape param with marker efficiency (theta_Am)
		
		//Prepare dimensions and data vectors
		int NOK0 = NOK[locind]; //obtain number of contributors (may be different for different markers)
		int NOU = *NOC - NOK0; //number of unknowns (may be different for markers)
		int nRep0 = nRep[locind]; //number of replicates for specific marker (may be different for different markers)
		
		//Prepare dimen
		int numGenos1p = nGenos[locind]; //Number of genotypes (1 contributor)
		int nAlleles = nAllelesVEC[locind]; //number of alleles
		//int nAlleles2 = nAlleles2VEC[locind]; //number of alleles (including potential stutters)
		//int nPS = nAlleles2 - nAlleles; //number of potential stutters
													
		int startIndMarker_nAlleles0 = startIndMarker_nAlleles[locind]; //index start number of alleles
		int startIndMarker_nAllelesReps0 = startIndMarker_nAllelesReps[locind]; //index start number of alleles (taking into account reps)
		
		int startIndMarker_nStutters0 = startIndMarker_nStutters[locind];
		int startIndMarker_outG1allele0 = startIndMarker_outG1allele[locind];
		int startIndMarker_outG1contr0 = startIndMarker_outG1contr[locind];
				
		//Prepare vector for known contributors (and also unknown): Need to know positions!
		vector<int> GindKnown(NOK0,0); //genotype index in vector 
		vector<int> kindKnown(NOK0,0); //contributor index in vector
		vector<int> kindUnknown(NOU,0); //contributor index in vector
		
		int aa, kk, rr, aaind; //indices for alleles, contributors and replicates
		int cc = 0; //counters for known
		int jj = 0; //counter for unknowns
		int knownGind0; //obtain genotype index of known 
		for(kk=0; kk< *NOC; kk++) { //for each contributors:
			knownGind0 = knownGind[ (*NOC)*locind + kk ]; //copy genotype index  (KNOWN)			
			if(knownGind0>=0) { //If contributor is known (genotype given)
				GindKnown[cc] = knownGind0; //copy genotype index
				kindKnown[cc] = kk; //insert contributor index
				cc++; //update counter for knowns
			} else { //if contributor is unknown (genotype not given)
				kindUnknown[jj] = kk; //insert contributor index
				jj++; //update counter for unknowns
			}
		}
				
		//Prepare variables (before sum-iterations):
		double nTyped = nTypedLong[locind];			
		vector<double> maTypedvec(nAlleles,0.0);
		vector<double> shapevK(nAlleles, 0.0); //Precalculation for known contributors:
		vector<double> shapev0(nAlleles, shape1); //degrad scaling of shape (init vector)
		bool useDeg = beta0 != 1.0; //whether to include degrad
		for (aa = 0; aa < nAlleles; aa++) { //traverse each observed alleles (indices), having PH>0		
			maTypedvec[aa] = maTypedLong[ startIndMarker_nAlleles0 + aa ]; //copy previously typed alleles
			if(useDeg) { //only apply degradation slope if necessary
				shapev0[aa] *= exp(basepairLong[ startIndMarker_nAlleles0 + aa ]*log(beta0) ); //scaling with degradation model (assumed already scaled)
			}
		}

		//Sum up contribution for each alleles (Taking into account mix proportions): ONLY CALCULATED FOR KNOWN CONTRIBUTORS INITIALLY
		for (kk = 0; kk < NOK0; kk++) { //for each known contributors 
			for (aa = 0; aa < nAlleles; aa++) { //traverse each alleles (indices)			
				shapevK[aa] += outG1contr[ startIndMarker_outG1contr0 + nAlleles*GindKnown[kk] + aa] * mixprop[kindKnown[kk]]; //contr from contr k to allele aa
			}
		}		

		//For each observations we traverse the modified "large sum" 				
		for (rr = 0; rr < nRep0; rr++) { //traverse each replicate		
			for (aa = 0; aa < nAlleles; aa++) { //traverse each observed alleles (indices), also the Q-allele. Potential not necessary!
				aaind = startIndMarker_nAllelesReps0 + aa*nRep0 + rr; //get index vectorized vector which includes (allele,rep) indices as [(1,1), (1,2), (1,3), (2,1), (2,2) etc]
				if( peaksLong[aaind] < AT0 ) { //no PH observation to consider
					continue; //skip if not observed
				}


				//CALCULATING LARGE SUM:
				double bigsum = 0.0;  //total sum over all genotypes

				#pragma omp parallel for reduction(+:bigsum) //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
				for (int iter = 0; iter < nJointCombs[locind]; iter++) { //for each combined/joint genotype outcome		
					//int iter = combUseVEC[ startIndMarker_nJointCombs0 + iter2]; //obtain correct iteration index to consider

					//Following variables must be declared for each iteration (non-shared):
					int a,k,r; //used to traverse alleles(a), contributors(k) and replicates(r)
					int aind; //used as index for alleles in genotypes
					vector<int> jointGind(NOU, 0); //This will be the permuation contribution index for the unknown inds (directly corresponds to indices of Gmarg)
					vector<double> shapev = shapevK; //make a copy of existing shapevector
					vector<double> maTypedvec2 = maTypedvec; //Creating copy of counter (per allele)
					double nTyped2 = nTyped; //creating copy of total counter

					//jointGind = digits(iter, base = numGenos1p, pad = NOU); #Equivalent operaion
					int startInd_alleles; //init for allele index
					double genoProd = 1.0; //calculating the genotype probability of the unknowns
					int modrest = iter; //used to keep remained after modulo (init as iter number)
					bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
					for (k = 0; k < NOU; k++) { //for each unknown contributors (summing up wrt both contr (outG1contr) and mx (mixprop)): Need each contr to derive shapev
						if (!inserted) { //if not all digits inserted
							if ( k>0 ) {
								modrest = int((modrest - jointGind[k - 1]) / numGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
							}
							jointGind[k] = modrest % numGenos1p; //INSERT NUMBER: convert number to "numGenos1p" basis 	
							if (modrest < numGenos1p) { //check if rest is smaller than base (only run if not inserted)
								inserted = true;  //then all digits are inserted
							}	
						} //else { jointGind[k] = 0; //zero pad	}
				
						//Sum up contribution for each alleles (Taking into account mix proportions)
						for (a = 0; a < nAlleles; a++) { //traverse each alleles (indices)
							shapev[a] += outG1contr[  startIndMarker_outG1contr0 + nAlleles*jointGind[k] + a] * mixprop[kindUnknown[k]]; //contr from contr k to allele a. NOTICE THE k+*NOK0 shift!
						}
										
						//CALCULATE GENOTYPE PROBS OF UNKNOWNS (MAY BE RELATED != -1)
						startInd_alleles = startIndMarker_outG1allele0 + 2*jointGind[k];
						for(a = 0;a<2; a++) { //loop through both alleles
							aind = outG1allele[ startInd_alleles + a];  //get allele index of genotype g_a
							genoProd *= (fst0*maTypedvec2[aind] + (1-fst0)*freqsLong[startIndMarker_nAlleles0 + aind]) / (1 + (nTyped2-1)*fst0); 
							maTypedvec2[aind] += 1; //update allele count for particular genotype
							nTyped2 += 1; //update total count
						}
						if( outG1allele[ startInd_alleles ] != outG1allele[ startInd_alleles+1 ] ) { // heterozygous variant
							genoProd *=2; //multiply by 2  if het. variant
						} 
					} //end for each unknown contributor

					//////////////////////////////////////////////
					//Calculating the inner sum-part begins here//
					//////////////////////////////////////////////

					//Scaling shape parameter with degrad model
					for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), having PH>0
						shapev[a] *= shapev0[a]; //scaling with degradation model (assumed already scaled)
					}
					
					vector<double> shapev2 = shapev;  //make a copy of existing shapevector
					
					//Scaling shape parameter with stutter model:			
					int stuttind,stuttind2,stuttFromInd,stuttToInd; //init vars
					double stuttExp; //this is the stutter proportion param to use
					for(stuttind=0; stuttind < nStutters[locind]; stuttind++) { 
						stuttind2 = startIndMarker_nStutters0 + stuttind; //obtain index for correct marker
						stuttFromInd = stuttFromIndVEC[stuttind2];  //NB: CAREFUL WITH INDEX
							
						if( shapev[stuttFromInd]>smalltol) { //ONLY NECESSARY TO PROVIDE MODIFICATION IF shapeval>0
							stuttToInd = stuttToIndVEC[stuttind2];  //NB: CAREFUL WITH INDEX
							
							//Obtaining expected stutter proportion (previously modelled) here:
							stuttExp = stuttExpVEC[stuttind2]; //NB: CAREFUL WITH INDEX					
							shapev2[stuttToInd] += stuttExp*shapev[stuttFromInd]; // #OBTAINED stutters
							shapev2[stuttFromInd] -= stuttExp*shapev[stuttFromInd]; // #SUBTRACTED stutters
						}
					}
						
					//Summing up contribution of each alleles:
					vector<int> nDropin(nRep0,0); //count number of drop-in for each replicate
					double logevidProb = 0.0; //evidence probability (weight)
					int cind; //cumulative allele indexing for PHs (for traversing over all replicates)
					double peak; //obtaining observeed peak height
					for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
						for(r = 0; r < nRep0; r++) { //traverse each replicates (observed alleles indicated by PH)
							cind = startIndMarker_nAllelesReps0 + a*nRep0 + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
							peak = peaksLong[cind]; //obtain coverage/peak
							
							//CHECK IF THE LOOPED INDEX (aaind) IS SAME AS LOOPED ALLELE (cind):
							if( aaind!=cind ) { //proceed as before if different index							
							
								//IF NOT DROPOUT
								if( peak > smalltol ){ //if PH>0 (this is same as PH>=AT since threshold has already been applied)
									if(shapev2[a] > smalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
										logevidProb += -lgamma(shapev2[a]) - shapev2[a]*const2 + (shapev2[a]-1)*log(peak) - peak*const1; //gamma/lgamma is found in cmath							
									} else { //If contribution and PH>0 	//DROPIN SET (C)	
										logevidProb += noiseSizeWeightLong[ cind ]; //likelihood for noise coverage
										nDropin[r] += 1; //count dropin for particular replicate															   
									}							
									
								//OTHERWISE IT IS DROPOUT
								} else if(shapev2[a] > smalltol) { //IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)
									logevidProb +=  pgamma(AT0,shapev2[a],scale0,1,1); //Add log(dropout-probability)
								}
								
							} else { //modify code if allele index is same (assures that Yvec[cind]>=(*AT) (see early in loop), hence no dropout possible
							
								if(shapev2[a]>smalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
									if( *maxY > smalltol ) {
										peak = *maxY; //replace PH with the one  given as argument
									}
									logevidProb += log( pgamma(peak, shapev2[a],scale0, 1, 0) - pgamma(AT0-1, shapev2[a],scale0, 1, 0));
										
								} else { //If contribution and PH>0 //DROPIN SET (C)
									logevidProb += noiseSizeWeightLong2[ cind ]; //likelihood for observed dropin (this will be the cumulative expression)						
									nDropin[r] += 1; //count dropin for particular replicate											
									//logevidProb += log( 1 - exp(- (*lambda)*(Yvec[cind]-*AT))); //calc logarithm of cumulative expression
								}		
							}
						} //end for each replicates (r)
					} //end for each observed alleles
					
								
					//Liklihood OF NUMBER OF NOISE/DROP-IN: Traverse each replicates (observed alleles indicated by PH)
					for(r = 0; r<nRep0; r++) {
						logevidProb += log(nNoiseParam0) + nDropin[r]*log(1-nNoiseParam0); //likelihood for number of noise alleles (geometrical distribution)															   													   
					}

					//FINAL INSERTION OF VALUES:
				   double val = exp(logevidProb)*genoProd; //calculate P(E|gj)P(gj)	
				   
				   //Calculate inner sum:
				   bigsum += val;
				   //bigsumVEC[startIndMarker_nJointCombs0 + iter] = val; //insert element
				} //end for each combination iterations (bigsum)
				pvalVEC[aaind] = bigsum; //insert p-value value	
			} //end for each allele
		} //end for eachreplicate
	} //end for each marker
} //end main function

} //end external