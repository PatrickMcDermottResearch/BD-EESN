sourceCpp("Functions/ensembleESNCPPFuncts_001.cpp")

######### QESN Model
embedInd=TRUE
quadInd=TRUE
featureLinkInd=FALSE




######ensemble length
ensembleLen=500

numLayers=7
# numLayers=6
# numLayers=5


################################################
########## Type of Dim. Red. ###################
################################################
#####type
dimRedESN="EOF"

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
mMin=0
mMax=6

tauEmb=3

####ridge parameters
ridgeMin=.0001
ridgeMax=.005
######Red dim
# redDimMin=20
redDimMin=6
redDimMax=20


######Second layer: number of hidden units
nhTwoMin=60
# nhTwoMax=35
nhTwoMax=75
#######delta
deltaMin=.10
deltaMax=1

deltaMinBoundRep=rep(deltaMin,numLayers)
deltamMaxBoundRep=rep(deltaMax,numLayers)

####curNhOne
curNhOne=84
hMatdimOne=curNhOne

isLogInput=TRUE
isLogOutput=TRUE
isExp=TRUE


########GA parameters
popSize=20
numGenerations=40
suggestions=c(4,.001,20,75,rep(1,numLayers))



gaObj=ga(type="real-valued",fitness = function(x)-deepGAFunction(x[1],x[2],x[3],x[4],c(x[5],x[6],x[7],x[8],x[9],x[10],x[11]),isLogInput,isLogOutput,isExp),suggestions=suggestions,lower=c(mMin,ridgeMin,redDimMin,nhTwoMin,deltaMinBoundRep),upper=c(mMax,ridgeMax,redDimMax,nhTwoMax,deltamMaxBoundRep),popSize = popSize,maxiter=numGenerations,monitor = TRUE )


bestParSet=summary(gaObj)$solution
bestParSet[1]=as.integer(bestParSet[1])
bestParSet[3]=as.integer(bestParSet[3])
bestParSet[4]=as.integer(bestParSet[4])


cat("Best Parameter Set: m: ",bestParSet[1], " ridge: ", bestParSet[2], " dim. red. ", bestParSet[3], " nhTwo:", bestParSet[4], " delta 1: ",bestParSet[5], " delta 2:  " , bestParSet[6], "\n")


