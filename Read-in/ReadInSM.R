#######Read in function file
source("Functions/masterBRNNFunctions_001.R")
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
library(GA)
dataSet="SoilMoisture"



#########################################################
############ Soil Moisture  #############################
#########################################################
if(dataSet=="SoilMoisture"){
  ####Set Embedding Indicator
  embedInd=TRUE
  
  ######tau
  tau=6
  # tau=12
  
  ##########################################
  ########### Get Soil Moisture Data #######
  ##########################################
  
  setwd("Datasets")
  rawNetSM=nc_open('SoilMoisture_data.nc')
  
  
  rawSM=ncvar_get(rawNetSM,"w")
  
  smFirstYear=1948
  strYrSM=1948
  
  #########Set out-of-sample years so we can calculate EOFs
  # strOutSampleYear=2014
  strOutSampleYear=2011
  tempTrainLen=length(strYrSM:strOutSampleYear)*12+5-tau
  inputInSamplSeq=1:(tempTrainLen-tau)
  
  tempReducedSM=rawSM[,,((strYrSM-smFirstYear)*12+1):dim(rawSM)[3]]
  
  numXLocs=dim(rawSM)[1]
  numYLocs=dim(rawSM)[2]
  
  
  smLat=seq(35.75,48.75,by=.5)
  smLon=seq(-101.75,-80.25,by=.5)
  
  lonLatMat=matrix(NA,numXLocs*numYLocs,2)
  
  for(i in 1:numYLocs){
    lonLatMat[((i-1)*numXLocs+1):(i*numXLocs),]=cbind(smLon,rep(smLat[i],numXLocs))
  }
  
  
  withNaNReducedSM=matrix(tempReducedSM,nrow=dim(tempReducedSM)[1]*dim(tempReducedSM)[2],ncol=dim(tempReducedSM)[3])
  
  
  nanIndex=unique(sapply(1:ncol(withNaNReducedSM),function(funVar)which(is.nan(withNaNReducedSM[,funVar]))))[[2]]
  
  reducedSM=withNaNReducedSM[-nanIndex,]
  
  numFullDimLocs=nrow(reducedSM)
  
  numPlotLocs=nrow(withNaNReducedSM)
  
  
  climateSM=matrix(NA,12,nrow(reducedSM))
  sdClimatologySM=matrix(NA,12,nrow(reducedSM))
  for(i in 1:12){
    climateSM[i,]=rowMeans(reducedSM[ ,seq(i,max(inputInSamplSeq),by=12)])
    sdClimatologySM[i,]=rowSds(reducedSM[ ,seq(i,max(inputInSamplSeq),by=12)])
  }
  
  
  
  anaomSM=matrix(NA,nrow(reducedSM),ncol(reducedSM))
  for(i in 1:ncol(reducedSM)){
    
    curMonthIndex=i%%12
    
    if(curMonthIndex==0){
      curMonthIndex=12
    }
    anaomSM[,i]=reducedSM[,i]-climateSM[curMonthIndex,]
  }
  
  # tempPlotVec=rep(NA,numPlotLocs)
  # tempPlotVec[-nanIndex]=anaomSM[,1]
  # image.plot(matrix(tempPlotVec,nrow=numXLocs,ncol=numYLocs))
  
  
  ##########Calculate SM EOFs
  # numSMEOFs=10
  numSMEOFs=15
  numLocs=numSMEOFs
  
  smEOFObj=calcEOFS(anaomSM,numSMEOFs,1:(max(inputInSamplSeq)+tau))
  
  phiSM=t(smEOFObj$phiMat)
  rawSMEOFs=t(smEOFObj$basisCoefs)
  leftOverPhiMat=t(smEOFObj$leftOverPhiMat)
  allSMEigs=smEOFObj$allEigVals
  
  #######################################
  ######## Get SST Data  ################
  #######################################
  
  srtYearSST=strYrSM
  #####last year for climate avg
  bgyrClimateAvg=1981
  endyrClimateAvg=2010
  dataSetType="Jan2017"
  
  SSTanom.data.obj=get.SST.data(srtYearSST,bgyrClimateAvg,endyrClimateAvg,dataSetType)
  SSTanom.dataRaw=SSTanom.data.obj$sst.anom
  sst.anom.land.sea=SSTanom.data.obj$sst.anom.land.sea
  sst.lat.long=SSTanom.data.obj$SSTlonlatTemp
  sst.long=SSTanom.data.obj$sst.long
  sst.lat=SSTanom.data.obj$sst.lat
  
  allLatLongSST=matrix(NA,nrow(sst.anom.land.sea),2)
  for(i in 1:length(sst.long)){
    allLatLongSST[((i-1)*length(sst.lat)+1):(i*length(sst.lat)),]=cbind(rep(sst.long[i],length(sst.lat)),as.vector(sst.lat))
  }
  
  
  #######changes spatial domain
  spatDomainIndex=which(sst.lat.long[,2]<=29 & sst.lat.long[,2]>=-29  )
  SSTanom.data=SSTanom.dataRaw[spatDomainIndex,]
  
  ######change spatial domain for plotting 
  SSTAllSpatDomainIndex=which(allLatLongSST[,2]<=29 & allLatLongSST[,2]>=-29  )
  SSTTruncAllLatLong=allLatLongSST[SSTAllSpatDomainIndex,]
  
  allSSTLongTrunc=sort(unique(SSTTruncAllLatLong[,1]))
  allSSTLatTrunc=sort(unique(SSTTruncAllLatLong[,2]))
  
  truncNxPlot=length(allSSTLongTrunc)
  truncNyPlot=length(allSSTLatTrunc)
  
  ######Indexs of sea locs (i.e., non NA values)
  plotIndexSST=which(!is.na(sst.anom.land.sea[SSTAllSpatDomainIndex,1]))
  
  
  ####### calculate SST EOFs
  numSSTEOFs=5
  
  sstEOFObj=calcEOFS(SSTanom.data,numSSTEOFs,inputInSamplSeq)
  rawSSTEOFs=sstEOFObj$basisCoefs
  sstPhi=sstEOFObj$phiMat
  
  
  #############################################
  ##### Create Train and Test Data ############
  #############################################
  
  rawDataInput=rawSSTEOFs
  numRawInputVars=ncol(rawDataInput)
  rawDataOutput=rawSMEOFs
  
  
  ###### Training Data
  xTrain=rawDataInput[inputInSamplSeq,]
  yTrain=t(rawDataOutput[,inputInSamplSeq+tau])
  yTrainFullDim=t(anaomSM[,inputInSamplSeq+tau])
  
  
  numRedPrecipPeriods=nrow(rawDataInput)+tau
  # numRedPrecipPeriods=nrow(rawDataInput)-tau
  
  ###### Test Data
  xTestIndex=(tempTrainLen):(numRedPrecipPeriods-tau)
  yTestIndex=(tempTrainLen+tau):(numRedPrecipPeriods)
  
  # xTest=t(rawDataInput[xTestIndex,])
  # yTestALLLocs=t(landPrecipData[,yTestIndex])
  
  xTest=rawDataInput[xTestIndex,]
  yTestALLLocs=t(anaomSM[,yTestIndex])
  
  testLen=nrow(xTest)
  trainLen=nrow(xTrain)
  
  ####inbetween indexes
  inbetweenIndexes=(max(inputInSamplSeq)+1):(min(xTestIndex)-1)
  inbetweenLen=length(inbetweenIndexes)
  
  ##### find months and years for out of sample forecasts
  monthSeq=sapply(1:testLen,function(funVar) findMonth( (yTestIndex)[funVar]  ) )
  yearSeq=sapply(1:testLen,function(funVar) findYear(strYrSM, (yTestIndex)[funVar]  ) )
  yrMonthLabels=sapply(1:testLen,function(funVar) yrMonthLabel( (yTestIndex)[funVar],yearSeq[funVar]  ) )
  
  mayIndex=which(monthSeq=="May")
  
  #############################################
  ##### Create Validation Data Index ##########
  #############################################
  
  
  validStrOutSampleYear=2005
  # validStrOutSampleYear=2012
  
  tempValidTrainLen=length(strYrSM:validStrOutSampleYear)*12+2-tau
  
  inputValidInSamplSeq=1:(tempValidTrainLen-tau)
  
  validLen=length(inputValidInSamplSeq)
  
  xValTestIndex=tempValidTrainLen:(yTestIndex[1]-1-tau)
  yValTestIndex=xValTestIndex+tau
  validTestLen=length(xValTestIndex)
  
  yValidALLLocs=t(anaomSM[,yValTestIndex])
  
  
  maySeq=seq(5,ncol(rawSMEOFs),by=12)
  
  validMayIndex=na.omit(match(maySeq,yValTestIndex))
  validInbetweenIndexes=(validLen+1):(min(xValTestIndex)-1)
  validInbetweenLen=length(validInbetweenIndexes)
  
  
  
  
  #############################################
  ##### Climate Means #########################
  #############################################
  
  
  outSampMonMeans=matrix(NA,testLen,numFullDimLocs)
  # outSampMonLw95=matrix(NA,testLen,numFullDimLocs)
  # outSampMonUp95=matrix(NA,testLen,numFullDimLocs)
  # outSampMonSD=matrix(NA,testLen,numFullDimLocs)
  for(i in 1:testLen){
    
    curMonthIndex=yTestIndex[i]%%12
    
    if(curMonthIndex==0){
      curMonthIndex=12
    }
    outSampMonMeans[i,]=climateSM[curMonthIndex,]
    # outSampMonLw95[i,]=inSampMonthLw95[curMonthIndex,]
    # outSampMonUp95[i,]=inSampMonthUp95[curMonthIndex,]
    # outSampMonSD[i,]=inSampMonthSD[curMonthIndex,]
    
  }
  
}