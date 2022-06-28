#library(MPSproto)

c = list()
c$NOC = as.integer(2)
c$NOK =  as.integer(2)
par = list()
par$mx = c(0.06087969, 0.93912031)
par$mu = 10179.32 
par$omega = 0.1780542
par$beta = 0.7984095
theta_Am = 2.069644

c$AT = 10
c$fst = 0.01
c$nNoiseParam =  0.4038462
c$noiseSizeWeight = c(-8.67574223651937,-9.6941689583176,-9.91539835448796,-13.260597607381,-4.11007264862848,-11.2438569927287,-13.2195753840045,-7.87417170064168,0)
c$nLocs = as.integer(1)
c$nSamples = as.integer(1)
c$nAlleles = as.integer(9 )
#c$nAlleles2 = 11
c$startIndMarker_nAlleles = as.integer(0)
c$startIndMarker_nAllelesReps = as.integer(0)
c$peakHeights = c(500, 1044, 1225, 13733, 18, 3199, 13332, 280, 0)
c$peakHeights2 = log(c$peakHeights)
c$peakHeights2[is.infinite(c$peakHeights2)] = 0 #will not be used

c$freq = c( 0.10679612, 0.23023578, 0.14840499, 0.09292649, 0.01109570, 0.10957004, 0.13869626, 0.01664355, 0.14563107)
c$nTyped = 2
c$maTyped =  c(0, 1, 0, 1, 0, 0, 2, 0, 0)
c$basepairs = c(0.25, 0.29, 0.33, 0.37, 0.41, 0.33, 0.37, 0.41, 0.57)
c$nGenos = as.integer(45)
c$nJointCombs = as.integer(c$nGenos^(c$NOC - c$NOK))
  
c$outG1allele = as.integer(c(0,0,0,1,0,2,0,3,0,4,0,5,0,6,0,7,0,8,1,1,1,2,1,3,1,4,1,5,1,6,1,7,1,8,2,2,2,3,2,4,2,5,2,6,2,7,2,8,3,3,3,4,3,5,3,6,3,7,3,8,4,4,4,5,4,6,4,7,4,8,5,5,5,6,5,7,5,8,6,6,6,7,6,8,7,7,7,8,8,8))
c$outG1contr = as.integer(c(2,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0,2,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,1,0,0,2,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,2,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1,0,0,0,1,0,0,0,0,0,2,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0,0,1,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,2))

c$nStutters = as.integer(8)
c$stuttFromInd = as.integer(c(0,1,2,3,4,5,6,7))
c$stuttToInd = as.integer(c(9,0,1,2,3,10,5,6))
c$stuttExp = rep(0.08341893,c$nStutters)
c$startIndMarker_nStutters <- c$startIndMarker_outG1allele <- c$startIndMarker_outG1contr <- as.integer(0)

c$knownGind = as.integer(c(14,27))

#GAMMA MODEL
calc = .C("loglikPrediction_allcomb_GA",as.numeric(0), c$nJointCombs,as.integer(c$NOC), as.integer(c$NOK),
          as.numeric(par$mx),  as.numeric(par$mu), as.numeric(par$omega), as.numeric(par$beta), as.numeric(theta_Am),
          as.numeric(c$AT),as.numeric(c$fst),as.numeric(c$nNoiseParam),as.numeric(c$noiseSizeWeight),
          c$nLocs, c$nSamples, c$nAlleles, c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps,
          c$peakHeights, c$peakHeights2, c$freq, c$nTyped, c$maTyped, c$basepairs,
          c$nGenos, c$outG1allele, c$outG1contr, c$startIndMarker_outG1allele, c$startIndMarker_outG1contr,
          c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttExp , c$startIndMarker_nStutters,
          c$knownGind)

expect_equal(round(calc[[1]],3),-66.833)
