#######Read in function file
source("/Users/PatrickMcDermott/Desktop/Dissertation_Code/RNN_Code/Master_Bayes_RNN/masterBRNNFunctions_001.R")
############# Libaries 
# library(shapes)
library(vegan)
library(pscl)
library(animation)
# library(mvnfast)
library(RSNNS)
library(randomForest)
library(fields)
library(truncnorm)
library(RcppTN)
library(matrixStats)
library(dtw)
library(verification)
library(deSolve)
library(glmnet)
library(foreach)
library(doParallel)
library(ncdf4)
library(abind)
library(NMF)
library(sp)
library(maptools)
library(shapefiles)
library(scoringRules)
library(GA)
# library(devtools)
# install_github(repo = "RcppTN",
#                username = "olmjo",
#                subdir = "pkg",
#                ref = "development"
# )

par(mar=c(5.1,4.1,4.1,2.1))


######################################
##### Select dataset  ################
######################################
# dataSet="SST"
# dataSet="SSTBerliner"
# dataSet="SSTEnsoCycle_2016"
# dataSet="Lorenz3"
# dataSet="Lorenz96"
dataSet="multiLorenz"
# dataSet="MJO_RMM1_RMM2"
# dataSet="unemployment"

#############set wd 
############# wd for SST Data
setwd("/Users/PatrickMcDermott/Documents/Data_Sets/SST")

#############set wd 
############# wd for MJO Data
# setwd("/Users/PatrickMcDermott/Documents/Data_Sets/MJO_Data")


#######set number of locs (or EOFs) 
# numLocs=10

if(dataSet=="MJO_RMM1_RMM2"){
  numLocs=2
  setwd("/Users/PatrickMcDermott/Documents/Data_Sets/MJO_Data")
}

if(dataSet=="SSTBerliner"  || dataSet=="SSTEnsoCycle_2016" ){
  numLocs=10
  setwd("/Users/PatrickMcDermott/Documents/Data_Sets/SST")
}

if(dataSet =="Lorenz96"){
  numLocs=40
  # numLocs=24
  # numLocs=10
}

if(dataSet =="multiLorenz"){
  numLocs=18
}

if(dataSet=="unemployment"){
  numLocs=12
}


############################
## Turn Data Types Off #####
############################
dataTypeSST=FALSE
dataTypeLorenz=FALSE
dataTypeMultiLorenz=FALSE
dataTypeMJO=FALSE
dataTypeUnemploy=FALSE


####################################################
############## Indicators for Data #################
####################################################
####################################
######## All Location Indicator ####
####################################
# allLocsSSTInd=TRUE
allLocsSSTInd=FALSE

set.seed(315)
# set.seed(180)

#########################################################################
########## Reduced Dimension Multi-Scale Lorenz-40 ######################
#########################################################################
if(dataSet=="multiLorenz"){
  
  
  #####turn datatype on 
  dataTypeMultiLorenz=TRUE
  
  ####Set Embedding Indicator 
  embedInd=FALSE
  
  #######set tau
  # tau=6
  # tau=3
  # tau=2
  tau=3
  
  numFullDimLocs=numLocs
  # truncLength=511+350
  # tempTrainLen=435-tau+1+350
  
  truncLength=511
  tempTrainLen=435-tau+1
  
  
  #############################################
  # Solve with external step size delta=.005
  # (internal step size dt=.005)
  #Larger delta more nonlinear
  #############################################
  # largeScaleVar=7.5
  # largeScaleVar=.05
  largeScaleVar=sqrt(.25)
  
  ###THIS IS THE BEST VALUE!!! deltaLorSim = .08
  deltaLorSim = .105
  # deltaLorSim = .08
  # deltaLorSim = .05
  
  # thetaLorSim = 18
  # thetaLorSim = 15
  # thetaLorSim = 18
  thetaLorSim = 10
  
  ######smaller dt, better approximation?
  dt = .0001
  # dt = .0005
  # TT = 500
  TT = 750
  
  lorenzM = deltaLorSim/dt
  rawNumTimeLor40 = TT/deltaLorSim + 1
  
  ######number of small scale locations
  # numSmScaleLocs=4*numLocs
  # numSmScaleLocs=16
  numSmScaleLocs=20
  
  numSmallScaleLocsTot=numFullDimLocs*numSmScaleLocs
  numTotLocs=numSmallScaleLocsTot+numFullDimLocs
  
  ################################################################
  ############## Wilks Parameterization ##########################
  ################################################################
  # hLorenzMult=1
  hLorenzMult=.35
  ### seperates two time scales
  cLorenzMult=10
  bLorenzMult=10
  
  ######################################################################
  ############## Alternative Parameterization ##########################
  ######################################################################
  hx=-1.90
  hy=1
  epsilonSclae=.045
  
  
  #######Large scale variables 
  rawDataTempLarge = matrix(NA,rawNumTimeLor40,numFullDimLocs)
  ######initial conditions
  rawDataTempLarge[1,] = rnorm(numFullDimLocs)
  
  #######Small scale variables 
  rawDataSmall = matrix(NA,rawNumTimeLor40,numSmallScaleLocsTot)
  rawDataSmall[1,]=rnorm(numFullDimLocs*numSmScaleLocs)
  
  isError=FALSE
  # isError=TRUE
  
  for (i in 2:rawNumTimeLor40){
    
    ####alternative parameterization
    multiLorObj=multiLorenz40CppAlt(rawDataTempLarge[i-1,],rawDataSmall[i-1,],thetaLorSim,dt,numFullDimLocs,lorenzM,numSmScaleLocs,hx,hy,epsilonSclae,largeScaleVar,numSmallScaleLocsTot,isError)
    # multiLorObj=multiLorenz40Cpp(rawDataTempLarge[i-1,],rawDataSmall[i-1,],thetaLorSim,dt,numFullDimLocs,lorenzM,numSmScaleLocs,hLorenzMult,cLorenzMult,bLorenzMult,largeScaleVar,numSmallScaleLocsTot)
    # cat(" it ", i , " large vals ", multiLorObj$xx , "\n")
    rawDataTempLarge[i,]=multiLorObj$xx
    rawDataSmall[i,]=multiLorObj$yy
    
  }
  
  strLorenz40=2500
  endLorenz40=3500
  
  ######large scale variation
  # sigmaEpsTr=2.5
  
  # sigmaEpsTr=.05
  # sigmaEpsTr=.50
  # sigmaEpsTr=sqrt(.10)
  sigmaEpsTr=sqrt(.25)
  # sigmaEpsTr=2.5
  
  
  ######small scale variation
  # sigmaEpsSmScale=.1
  # sigmaEpsSmScale=1
  
  ######for now we will only use the large scale variables 
  ###### with data model error
  rawDataLorenz40Lrg=matrix(rlnorm(length(strLorenz40:endLorenz40)*numFullDimLocs,abs(as.vector(rawDataTempLarge[strLorenz40:endLorenz40,]/2)),sigmaEpsTr),nrow=length(strLorenz40:endLorenz40),ncol=numFullDimLocs)
  # rawDataLorenz40Sm=rawDataSmall[strLorenz40:endLorenz40,]+matrix(rnorm(length(strLorenz40:endLorenz40)*numFullDimLocs*numSmScaleLocs,0,sigmaEpsSmScale),nrow=length(strLorenz40:endLorenz40),ncol=numFullDimLocs*numSmScaleLocs)
  
  ###### without data model error
  # rawDataLorenz40Lrg=rawDataTempLarge[strLorenz40:endLorenz40,]
  # rawDataLorenz40Sm=rawDataSmall[strLorenz40:endLorenz40,]
  
  rawDataLorenz40Comb=cbind(rawDataLorenz40Lrg,rawDataLorenz40Sm)
  
  rawData=rawDataLorenz40Lrg
  # rawData=rawDataLorenz40Sm[,1:numLocs]
  
  
  #############################################
  ##### Create Train and Test Data ############
  #############################################
  rawDataInput=rawData
  inputInSamplSeq=1:(tempTrainLen-tau+1)
  
  ###### Training Data
  xTrain=rawData[inputInSamplSeq,]
  yTrain=rawData[(tau+1):(tempTrainLen+1),]
  yTrainFullDim=log(yTrain)
  
  ######Plot Objects
  plotLorenzYInsamp=(rawDataSmall[strLorenz40:endLorenz40,])[(tau+1):(tempTrainLen+1),]
  plotLorenzXInsamp=(rawDataTempLarge[strLorenz40:endLorenz40,])[(tau+1):(tempTrainLen+1),]
  
  ###### Test Data
  xTestIndex=(tempTrainLen+1):(truncLength-tau)
  yTestIndex=(tempTrainLen+tau+1):truncLength
  
  xTest=rawData[xTestIndex,]
  yTest=rawData[yTestIndex,]
  
  testLen=nrow(xTest)
  trainLen=nrow(xTrain)
  
  ####inbetween indexes
  inbetweenIndexes=(max(inputInSamplSeq)+1):(min(xTestIndex)-1)
  lenInbetweenSeq=length(inbetweenIndexes)
  
  #############################################
  ##### Create Validation Data Index ##########
  #############################################
  
  validLenDif=100
  validLen=trainLen-validLenDif
  
  xValTestIndex=(tempTrainLen+1-validLenDif):(truncLength-tau-validLenDif)
  yValTestIndex=(tempTrainLen+tau+1-validLenDif):(truncLength-validLenDif)
  
  validTestLen=length(xValTestIndex)
  
  validInbetweenIndexes=(validLen+1):(min(xValTestIndex)-1)
  validInbetweenLen=length(validInbetweenIndexes)
  

}




for(i in 1:6){
  
  dev.new()
  cairo_ps("/Users/PatrickMcDermott/Desktop/Figure2.eps",height = 13,width=16)
  split.screen( rbind(c(0, 1,.45,1), c(0,1,.01,.45)))
  split.screen(c(2,3), screen=1)-> ind
  
  
  for(j in 1:3){
    screen( ind[j])
    par(mar=c(2.35,2.25,2,1.5))
    par(mgp=c(1.45,.5,0))
    plot(yTrain[1:100,((i-1)*3+j)],type='l',main=paste("(a) z (Large-Scale) Realization; Location: ",((i-1)*3+j)),xlab="Period",ylab="State",lwd=2,cex.axis=.95,cex.lab=.95,cex.main=1.25)
  }
  for(j in 1:3){
    screen( ind[3+j])
    par(mar=c(2.35,2.25,2,1.5))
    par(mgp=c(1.25,.5,0))
    plot(plotLorenzXInsamp[1:100,((i-1)*3+j)],type='l',main=paste("(b)  x (Large-Scale)  Realization; Location: ",((i-1)*3+j)),xlab="Period",ylab="State",lwd=2,cex.axis=.95,cex.lab=.95,cex.main=1.24,col="deepskyblue")
  }
  
  # screen(2)
  # par(mar=c(0,3.5,.25,6))
  split.screen(c(2,6), screen=2)-> indTwo
  
  
  for(plotIt in 1:2){
    for(j in 1:3){
      
      for(k in 1:2){
        screen( indTwo[(j-1)*2+k+(plotIt-1)*6 ])
        par(mar=c(1.5,1.25,2,1.25))
        par(mgp=c(1.25,.5,0))
        tempLorenzPlotIndex=((i-1)*3+j)
        
        if(j==1 & k==1){
        titleDeepLorenzPlot=paste("(c) y ; Large Loc:",((i-1)*3+j), "; Small Loc:", (k+(plotIt-1)*2) )
        tilteSize=.665
        }else{
          titleDeepLorenzPlot=paste("y (Small); Large Loc:",((i-1)*3+j), "; Small Loc:", (k+(plotIt-1)*2) )
          tilteSize=.665
        }
        
        plot(plotLorenzYInsamp[1:100,(numSmScaleLocs*(tempLorenzPlotIndex-1)+k+(plotIt-1)*2)],type='l',main=titleDeepLorenzPlot,xlab="Period",ylab="State",lwd=2,col="red",cex.axis=.85,cex.lab=.85,cex.main=.85)
      }
    }
  }
  close.screen( all=TRUE)
  dev.off()
}


