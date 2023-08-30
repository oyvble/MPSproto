#TODO:
- Optional pre-processing step for MPS sequence data: BayesHammer for sequence correction of noise.
- Optional to have a dynamic AT (per-sample)?


MPSproto v0.9.3 (Release date: 2023-08-23)
=============================================
- Fixed bugs:
	- The degradation model was turned on even when the user had selected "No" when having kit selected (gui updated).
	- The program chashed for calculation under Hd when a reference had missing markers (getUpperLR updated).
    	
- Updated to handle missing evidence markers (prepareData_prediction,prepareC_prediction).

MPSproto v0.9.2 (Release date: 2023-08-23)
=============================================
- Avoid program to crash when using degradation model and having kit selected, and "CE:" prefix format is not used.
- Improved robusness of generation of noise in genMPSevidence.
- An error is given if the number of alleles and reads are not the same.

MPSproto v0.9.1 (Release date: 2023-05-11)
=============================================
- Fixed bug when assigning values in toolbar in GUI. 
-> Adding "ncolumns=1" at gui.R:L279.

MPSproto v0.9.0 (Release date: 2022-08-23)
=============================================
- Including a GUI for interpretation (gui.R)

MPSproto v0.8.1 (Release date: 2022-07-29)
=============================================
- Fixed bug in getStutterData causing stutter from homozygous genotypes to be removed.

MPSproto v0.8.0 (Release date: 2022-03-15)
=============================================
- Implementing Negative Binomial (NB) as an model option in inferEvidence function: model=("GA,"NB")
	- Following functions changed: calcLogLikC_prediction, prepareC_prediction, prefit_prediction, validMLE, calcLogLikR_prediction
	- C++ now utilizing the function pnbinom_mu for cumulative calculation of NB
- genMPSevidence: Degradation not possible when generating evidence for MPS
- fitgammamodel renamed as fitSUMmodel which fits all implemented models (GA/NB).

- Updated functions for stutter calibration:
	- getStutterData function is more comprehensive (also take minStuttOccurence argument)
	- getFilteredData replaced by getStutterData (deprecated)
	- inferStuttRegModel renamed to inferStutterModel

MPSproto v0.7.1 (Release date: 2022-03-02)
=============================================
- Removed function getComplement (not used)

MPSproto v0.7.0 (Release date: 2022-02-02)
=============================================
- Included tutorial and scripts for calibrating data from paper

MPSproto v0.6 (Release date: 2022-02-14)
=============================================
- Now also supporting CE format (Stutter model for BLMM can be included): getStutterIndex updated
- Take into account that getStutterIndex returns NULL if none of the alleles can have a certain stutter type 
- Added functions:
	- getStutterTypeRule: Get motif change and motif rank for given stutter type
	- genMPSevidence: Used to generate MPS mixture based on a calibrated model

MPSproto v0.5 (Release date: 2022-01-12)
=============================================
- Updated inferStuttRegModel function to also include validation plots. Also supports precision expression.

MPSproto v0.4 (Release date: 2021-12-03)
=============================================
- n0 stutters are defined differently (getStutterIndex changed).
- normalising/minFreq argument can be modified by including it into the hypotheses object

MPSproto v0.3 (Release date: 2021-11-12)
=============================================
- Including inferEvidence2 function to also estimate marker efficiency, penalized with it's defined prior.
- Increase the theoretical maximum for validMLE (take into account marker efficiency)
- Including helpfunctions for easy extraction of Param, LR, AIC.
- Fixed bug in plotMPS: Marker specific AT/ST argument was not shown.

MPSproto v0.2 (Release date: 2021-11-03)
=============================================
- Changing NoiseSize model to be discrete instead of continuous (still Pareto).
- inferMarkerEfficiency function can now directly perform MCMC if mleObj argument is provided.
- Including logliki function to obtain per-marker logLikelihood.
- Including plotTopMPS for checking model (also gives  stutter contribution per contributors)
- Including getStutterExpectations to obtain expected stutter proportions between each alleles in evidence data

MPSproto v0.1 (Release date: 2021-10-19)
=============================================
- init release

