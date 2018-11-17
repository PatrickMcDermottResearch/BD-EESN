sourceCpp("Functions/ensembleESNCPPFuncts_001.cpp")


######### QESN Model
embedInd=TRUE
quadInd=FALSE
featureLinkInd=TRUE


######ensemble length
ensembleLen=500

#####Number of Layers
numLayers=3


################################################
########## Type of Dim. Red. ###################
################################################
#####type
dimRedESN="EOF"
# dimRedESN="Fourier"

#####scaling ind.
scaleByHInd=TRUE

if(dimRedESN=="EOF"){
  scaleByHInd=FALSE
}


######################################
##### Layer One ######################
######################################

#####Sparse Percentage
###### W sparseness
piWESNOne=.10

###### U sparseness
piUESNOne=.10

#####wWidth
wWidthOne=.10
uWidthOne=.10


######################################
##### Layer Two ######################
######################################

#####Sparse Percentage
###### W sparseness
piWESNTwo=.10

###### U sparseness
piUESNTwo=.10

#####wWidth
wWidthTwo=.10
uWidthTwo=.10

curMTwo=0

#################################
### GA Parameters ###############
#################################
#######MSeq
mMin=2
mMax=4

tauEmb=6

####ridge parameters
ridgeSeq=c(.01,.05)

######Red dim
redDimMin=6
redDimMax=8

######Second layer: number of hidden units
nhTwoMin=25
nhTwoMax=35
# nhTwoMax=120

#######delta
deltaMin=.10
deltaMax=1

deltaMinBoundRep=rep(deltaMin,numLayers)
deltamMaxBoundRep=rep(deltaMax,numLayers)

########curNhOne
curNhOne=84
hMatdimOne=curNhOne


########GA parameters
popSize=20
numGenerations=40
suggestions=c(3,.05,6,25,1,.5,1)


gaObj=ga(type="real-valued",fitness = function(x)-smGAFunction(x[1],x[2],x[3],x[4],c(x[5],x[6],x[7])),suggestions=suggestions,lower=c(mMin,ridgeSeq[1],redDimMin,nhTwoMin,deltaMinBoundRep),upper=c(mMax,ridgeSeq[2],redDimMax,nhTwoMax,deltamMaxBoundRep),popSize = popSize,maxiter=numGenerations,monitor = TRUE )


bestParSet=summary(gaObj)$solution
bestParSet[1]=as.integer(bestParSet[1])
bestParSet[3]=as.integer(bestParSet[3])
bestParSet[4]=as.integer(bestParSet[4])


cat("Best Parameter Set: m: ",bestParSet[1], " ridge: ", bestParSet[2], " dim. red. ", bestParSet[3], " nhTwo:", bestParSet[4], " delta 1: ",bestParSet[5], " delta 2:  " , bestParSet[6], " delta 3:  " , bestParSet[7], "\n")
