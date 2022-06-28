//This script contains functions for calculating likelihood values.
/*ABOUT:
- No structuring of data is needed (already done in prepareC)
- All markers are looped outside the inner large-sum loop.
- Possible to traverse certain genotype combinations

-> Stores likelihood values for each genotype combination in parallelization
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
void loglikPrediction_allcomb2_GA(double *bigsumVEC, int *nJointCombs, int *NOC, int *NOK, 
	double *mixprop, double *mu, double *omega, double *beta, double *theta_Am, //eps is stutter proportion vector (1 element for each kind)
	double *AT, double *fst, double *nNoiseParam, double *noiseSizeWeightLong, //Noise parameters are given
   	int *nMarkers, int *nRep, int *nAllelesVEC, int *startIndMarker_nAlleles, int *startIndMarker_nAllelesReps,
	double *peaksLong, double *peaksLong2, double *freqsLong, double *nTypedLong, double *maTypedLong, double *basepairLong,
	int *nGenos, int *outG1allele, int *outG1contr, int *startIndMarker_outG1allele, int *startIndMarker_outG1contr,
	int *nStutters, int *stuttFromIndVEC, int *stuttToIndVEC, double *stuttExpVEC, int *startIndMarker_nStutters, 	
	int *knownGind, int *startIndMarker_nJointCombs) {
	
	
	int numThreads = thread::hardware_concurrency();
	int useThreads = numThreads; //number of threads to use
	omp_set_num_threads(useThreads);  //set number of threads to use 
	
	//parameter values (per marker):
	double mu0 = mu[0]; //assume constant
	double omega0 = omega[0]; //assume constant
	double beta0 = beta[0]; //assume constant

	//Prepare transformation of parameters:
	double omegasq = omega0*omega0;
	double scale0 = mu0*omegasq; //obtain scale param (same for all obs)
	//double const1 = 1/scale0; //constant 1 (used in pgamma,dgamma)
	double const0 = log(scale0); //constant 2 (used in dgamma)
	double shape0 = 1/omegasq; // //theta_Am[locind]/omegasq; //shape for 'full het allele'    

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
		int startIndMarker_nJointCombs0 = startIndMarker_nJointCombs[locind]; //obtain start index
		int startIndMarker_nStutters0 = startIndMarker_nStutters[locind];
		int startIndMarker_outG1allele0 = startIndMarker_outG1allele[locind];
		int startIndMarker_outG1contr0 = startIndMarker_outG1contr[locind];
				
		//Prepare vector for known contributors (and also unknown): Need to know positions!
		vector<int> GindKnown(NOK0,0); //genotype index in vector 
		vector<int> kindKnown(NOK0,0); //contributor index in vector
		vector<int> kindUnknown(NOU,0); //contributor index in vector
		
		int aa, kk; //indices for alleles, contributors and replicates
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

		//PREPARE NEW VARIABLES (USED IN pdf calculation to reduce number of calls)
		vector<double> Wa(nAlleles*nRep0,0); //prepare data vector ya/scale				
		for (aa = 0; aa < nAlleles; aa++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
			for(kk = 0; kk < nRep0; kk++) { //traverse each replicates (observed alleles indicated by PH)
				jj = aa*nRep0 + kk; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
				Wa[jj] = peaksLong[startIndMarker_nAllelesReps0 + jj]/scale0; //obtain updated
			}
		}		
		
		//CALCULATING LARGE SUM:
		//double bigsum = 0.0;  //total sum over all genotypes

		#pragma omp parallel for //reduction(+:bigsum) //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
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
			int cind1,cind2; //cumulative allele indexing for PHs (for traversing over all replicates)
			double peak; //obtaining observeed peak height
			for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
				for(r = 0; r < nRep0; r++) { //traverse each replicates (observed alleles indicated by PH)
					cind1 = a*nRep0 + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
					cind2 = startIndMarker_nAllelesReps0 + cind1; //get marker corrected index (since long vector)
					peak = peaksLong[cind2]; //obtain coverage/peak
					
					//IF NOT DROPOUT
					if( peak > smalltol ){ //if PH>0 (this is same as PH>=AT since threshold has already been applied)
						if(shapev2[a] > smalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)							
							logevidProb += -lgamma(shapev2[a]) - shapev2[a]*const0 + (shapev2[a]-1)* peaksLong2[cind2] - Wa[cind1]; //gamma/lgamma is found in cmath														
							//logevidProb += -lgamma(shapev2[a]) - shapev2[a]*const2 + (shapev2[a]-1)*log(peak) - peak*const1; //gamma/lgamma is found in cmath							
						} else { //If contribution and PH>0 	//DROPIN SET (C)	
							logevidProb += noiseSizeWeightLong[ cind2 ]; //likelihood for noise coverage
							nDropin[r] += 1; //count dropin for particular replicate															   
						}							
						
					//OTHERWISE IT IS DROPOUT
					} else if(shapev2[a] > smalltol) { //IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)
						logevidProb +=  pgamma(AT0,shapev2[a],scale0,1,1); //Add log(dropout-probability)
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
		   //bigsum += val;
		   bigsumVEC[startIndMarker_nJointCombs0 + iter] = val; //insert element
		} //end for each combination iterations (bigsum)
		//logLik[0] += log(bigsum); //calculate logLik by adding log(innerSUM)	

	} //end for each marker
} //end main function


/////////////////////////////////
//IMPLEMENTATION OF THE NB MODEL/
/////////////////////////////////


void loglikPrediction_allcomb2_NB(double *bigsumVEC, int *nJointCombs, int *NOC, int *NOK, 
	double *mixprop, double *mu, double *omega, double *beta, double *theta_Am, //eps is stutter proportion vector (1 element for each kind)
	double *AT, double *fst, double *nNoiseParam, double *noiseSizeWeightLong, //Noise parameters are given
   	int *nMarkers, int *nRep, int *nAllelesVEC, int *startIndMarker_nAlleles, int *startIndMarker_nAllelesReps,
	double *peaksLong, double *peaksLong2, double *freqsLong, double *nTypedLong, double *maTypedLong, double *basepairLong,
	int *nGenos, int *outG1allele, int *outG1contr, int *startIndMarker_outG1allele, int *startIndMarker_outG1contr,
	int *nStutters, int *stuttFromIndVEC, int *stuttToIndVEC, double *stuttExpVEC, int *startIndMarker_nStutters, 	
	int *knownGind, int *startIndMarker_nJointCombs) {
	
	
	int numThreads = thread::hardware_concurrency();
	int useThreads = numThreads; //number of threads to use
	omp_set_num_threads(useThreads);  //set number of threads to use 
	
	//parameter values (per marker):
	double mu0 = mu[0]; //assume constant
	double omega0 = omega[0]; //assume constant
	double beta0 = beta[0]; //assume constant
	double size0 = mu0/(mu0*omega0*omega0 - 1); //Obtain size param from mean and cv
	double const0 = size0*log(size0) -  lgamma(size0); //constant1 - constant2

	//Calculating the loglik (over all markers), which is returned
	const double smalltol = 1.0e-30; //a tiny number > 0 (avoiding zero roundoff errors)
	for(int locind=0; locind< *nMarkers; locind++) {	//for each marker:
	
		//OBTAIN CONSTANTS FOR SPECIFIC MARKER:
	
		//default settings
		double fst0 = fst[locind]; //theta-correction param
		double AT0 = AT[locind]; //analytical threshold
	    double nNoiseParam0 = nNoiseParam[locind]; //Number of noise count param (geometric)
		
		//SHAPE IS NOW EQUIVALENT TO EXPECTATION!
		double shape1 = mu0*theta_Am[locind]; //scale mu param with marker efficiency (theta_Am): 
		
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
		int startIndMarker_nJointCombs0 = startIndMarker_nJointCombs[locind]; //obtain start index
		int startIndMarker_nStutters0 = startIndMarker_nStutters[locind];
		int startIndMarker_outG1allele0 = startIndMarker_outG1allele[locind];
		int startIndMarker_outG1contr0 = startIndMarker_outG1contr[locind];
				
		//Prepare vector for known contributors (and also unknown): Need to know positions!
		vector<int> GindKnown(NOK0,0); //genotype index in vector 
		vector<int> kindKnown(NOK0,0); //contributor index in vector
		vector<int> kindUnknown(NOU,0); //contributor index in vector
		
		int aa, kk; //indices for alleles, contributors and replicates
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

		//PREPARE NEW VARIABLES (USED IN pdf calculation to reduce number of calls)
		vector<double> Wa(nAlleles*nRep0,0); //prepare data vector ya + omega
		vector<double> Za(nAlleles*nRep0,0); //prepare data vector lgamma(ya + omega)=lgamma(Wa)
		for (aa = 0; aa < nAlleles; aa++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
			for(kk = 0; kk < nRep0; kk++) { //traverse each replicates (observed alleles indicated by PH)
				jj = aa*nRep0 + kk; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
				Wa[jj] = peaksLong[startIndMarker_nAllelesReps0 + jj] + size0; //obtain updated
				Za[jj] = lgamma(Wa[jj]); //obtain updated
			}
		}					
		
		//CALCULATING LARGE SUM:
		//double bigsum = 0.0;  //total sum over all genotypes
		#pragma omp parallel for //reduction(+:bigsum) 
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
			int cind1,cind2; //cumulative allele indexing for PHs (for traversing over all replicates)
			double peak; //obtaining observeed peak height
			for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
				for(r = 0; r < nRep0; r++) { //traverse each replicates (observed alleles indicated by PH)
					cind1 = a*nRep0 + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
					cind2 = startIndMarker_nAllelesReps0 + cind1; //get marker corrected index (since long vector)
					peak = peaksLong[cind2]; //obtain coverage/peak
					
					//IF NOT DROPOUT
					if( peak > smalltol ){ //if PH>0 (this is same as PH>=AT since threshold has already been applied)
						if(shapev2[a] > smalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)				
							logevidProb += peak*log(shapev2[a]) - Wa[cind1]*log(shapev2[a]+size0) + const0 + Za[cind1] - peaksLong2[cind2]; 	//pdf of NegBinom						
						} else { //If contribution and PH>0 	//DROPIN SET (C)	
							logevidProb += noiseSizeWeightLong[ cind2 ]; //likelihood for noise coverage
							nDropin[r] += 1; //count dropin for particular replicate															   
						}							
						
					//OTHERWISE IT IS DROPOUT
					} else if(shapev2[a] > smalltol) { //IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)
						logevidProb += pnbinom_mu(AT0-1,size0,shapev2[a],1,1); //Add log(dropout-probability)
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
		   //bigsum += val;
		   bigsumVEC[startIndMarker_nJointCombs0 + iter] = val; //insert element
		} //end for each combination iterations (bigsum)
		//logLik[0] += log(bigsum); //calculate logLik by adding log(innerSUM)	

	} //end for each marker
} //end main function


} //end external