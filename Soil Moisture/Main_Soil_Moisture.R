#Read in funciton file 
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
library(scoringRules)
library(GA)

sourceCpp("Functions/ensembleESNCPPFuncts_001.cpp")
######### QESN Model
embedInd=TRUE
quadInd=FALSE
featureLinkInd=TRUE




####################################################
##### Process Results for CV #######################
####################################################
####find min MSE
####use Nino3.4 as the metric for SST
# curM=0
curM=3
tauEmb=6

#######################################################
##### Create Current Dataset ##########################
#######################################################


###Variables for out-of-sample test
curTrainLen=trainLen
curXTestIndex=xTestIndex
curYTestLen=length(yTestIndex)
curInbetweenIndexes=inbetweenIndexes
curInbetweenLen=inbetweenLen
curTestLen=testLen
curYTestAllLocs=yTestALLLocs
curMayIndex=mayIndex

createESNDataObj=createSMEmbedRNNData(curTrainLen,curM,tauEmb,yTrain,rawDataInput,numRawInputVars,curXTestIndex,curTestLen,curInbetweenIndexes,curInbetweenLen,numLocs)

lenInSampleEmb=createESNDataObj$lenInSampleEmb
addScaleMat=createESNDataObj$addScaleMat
scaleFactorMean=addScaleMat[1,1]
scaleFactor=createESNDataObj$curInSampSD
designMatrix=createESNDataObj$designMatrix
designMatrixOutSample=createESNDataObj$designMatrixOutSample

designMatrixInbetween=createESNDataObj$inbetweenDesignMatrix
yScaledTrain=createESNDataObj$curYInSamp

################################################
################################################
##########  Simulate Reservoirs ################
################################################
################################################


######ensemble length
ensembleLen=100
layerOneDimRed=6

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

curNhOne=84
hMatdimOne=curNhOne
numLayers=2
# numLayers=3


curNhTw0=25
curMTwo=0

if(numLayers==2){
  deltaVec=c(0.9990464,0.9949445)
}
if(numLayers==3){
  deltaVec=c(0.9997245,0.5050391,0.9812467)
}
if(numLayers==2){
curRidge=0.04720978
}
if(numLayers==3){
  curRidge=0.04651549
}

##########first layer
setParObjFirst=setParsEESN(curRidge ,curNhOne,numLocs,curM)
sampVecESNOne=setParObjFirst$sampVecESN
stratValuesXtempOne=setParObjFirst$stratValuesXtemp


setParObjTwo=setParsEESN(curRidge ,curNhTw0,layerOneDimRed,curMTwo)
sampVecESNTwo=setParObjTwo$sampVecESN
stratValuesXtempTwo=setParObjTwo$stratValuesXtemp

if(featureLinkInd){
  if(quadInd){
    numVPars=2*curNhTw0+(numLayers-1)*layerOneDimRed
    ridgeMat=diag(curRidge,numVPars)
  }else{
    numVPars=curNhTw0+(numLayers-1)*layerOneDimRed
    ridgeMat=diag(curRidge,numVPars)
  }
}else{

  if(quadInd){
    numVPars=2*curNhTw0
    ridgeMat=diag(curRidge,2*numVPars)

  }else{
    numVPars=curNhTw0
    ridgeMat=diag(curRidge,numVPars)
  }

}



nColsUOne=ncol(designMatrix)
nColsUTwo=layerOneDimRed+1


#############in-sample and out-sample ESN mats
esnInSampArray=array(NA,c(ensembleLen,lenInSampleEmb,numLocs))
esnOutSampArray=array(NA,c(ensembleLen,curTestLen,numLocs))
esnOutSampArrayAll=array(NA,c(ensembleLen,curTestLen,numFullDimLocs))
##############################################
################ collect HMats ###############
##############################################

# hArrayInSamp=array(NA,c(ensembleLen,2*curNhTw0,lenInSampleEmb))
# hArrayOutSamp=array(NA,c(ensembleLen,2*curNhTw0,curTestLen))

hArrayInSamp=matrix(NA,lenInSampleEmb,numVPars*ensembleLen)
hArrayOutSamp=matrix(NA,curTestLen,numVPars*ensembleLen)

###For plotting purposes only
if(numLayers==2){
plotHOneArrayInSamp=matrix(NA,lenInSampleEmb,curNhOne*ensembleLen)
plotInSampFeatLinks=matrix(NA,lenInSampleEmb,layerOneDimRed*ensembleLen)
}

for(iEnsem in 1:ensembleLen){


  hMatObj=multiLayerDeepESNFor(curNhOne,lenInSampleEmb,designMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUOne,deltaVec,stratValuesXtempOne,curTestLen,dimRedESN,layerOneDimRed,curNhTw0,sampVecESNTwo,uWidthTwo,piWESNTwo,piUESNTwo,wWidthTwo,nColsUTwo,stratValuesXtempTwo,ridgeMat,yScaledTrain,scaleFactor,addScaleMat,designMatrixOutSample,designMatrixInbetween,curInbetweenLen,featureLinkInd,numLayers )

 #####Collect h Mats
  hArrayInSamp[,(1+(iEnsem-1)*(numVPars)):(iEnsem*(numVPars))]=t(hMatObj$inSampHMat)
  hArrayOutSamp[,(1+(iEnsem-1)*(numVPars)):(iEnsem*(numVPars))]=t(hMatObj$outSampHMat)

  ###For plotting purposes only
  if(numLayers==2){
  plotHOneArrayInSamp[,(1+(iEnsem-1)*(curNhOne)):(iEnsem*(curNhOne))]=t(hMatObj$hMatOneInSamp)
  plotInSampFeatLinks[,(1+(iEnsem-1)*(layerOneDimRed)):(iEnsem*(layerOneDimRed))]=t(hMatObj$inSampleFeatureLinks)
  }
 #####For comparison, calculate in-sample and out-sample MSE
  esnVMat=solve(hMatObj$inSampHMat%*%t(hMatObj$inSampHMat)+diag(curRidge,numVPars))%*%hMatObj$inSampHMat%*%yScaledTrain
  esnInSampArray[iEnsem,,]=t(hMatObj$inSampHMat)%*%esnVMat
  esnOutSampArray[iEnsem,,]=(t(hMatObj$outSampHMat)%*%esnVMat)*scaleFactor+scaleFactorMean

  esnOutSampArrayAll[iEnsem,,]=esnOutSampArray[iEnsem,,]%*%t(phiSM)

  print(iEnsem)
}

########ESN in-sample and out-sample MSE
esnInSampEnsemMean=sapply(1:numLocs,function(funVar) colMeans(esnInSampArray[,,funVar]))
esnOutSampEnsemMean=sapply(1:numLocs,function(funVar) colMeans(esnOutSampArray[,,funVar]))



esnInSampMSE=mean((esnInSampEnsemMean-yScaledTrain)^2)

cat("ESN In-Sample MSE ",esnInSampMSE , "\n" )


allLocsESNOutSamp=t(phiSM%*%t(esnOutSampEnsemMean))
esnOutSampMSE=mean((allLocsESNOutSamp[curMayIndex,]-curYTestAllLocs[curMayIndex,])^2)

cat("ESN Out-Sample MSE ",esnOutSampMSE , "\n" )




postMeanAll=matrix(NA,numFullDimLocs,length(curMayIndex))
postSDAll=matrix(NA,numFullDimLocs,length(curMayIndex))
sumCRPSALL=0


for(i in 1:length(curMayIndex)){

  postMeanAll[,i]=apply(esnOutSampArrayAll[,curMayIndex[i],],2,mean)
  postSDAll[,i]=apply(esnOutSampArrayAll[,curMayIndex[i],],2,sd)
  curCRPS=crps(yTestALLLocs[curMayIndex[i],],cbind(postMeanAll[,i],postSDAll[,i]))$CRPS
  sumCRPSALL=sumCRPSALL+curCRPS
  print(i)
}


cat("All Locs CRPS", sumCRPSALL , "\n" )


par(mfrow=c(1,1))
par(mar=c(4.5,4.5,3.5,2.1))
skillScoreVar=1-colMeans((t(postMeanAll) - curYTestAllLocs[curMayIndex,] )^2)/colMeans((curYTestAllLocs[curMayIndex,])^2)

cat("D-EESN: % Skill Score >0 ", length(which(skillScoreVar>0))/numFullDimLocs , "\n")



sourceCpp("Functions/DB-EESNNCPPFunctions.cpp")
library(mvnfast)




##########################################
### Calculate Vars. for Data Model   #####
##########################################

sigmaZMat=matrix(0,numFullDimLocs,numFullDimLocs)
numNonEof=numFullDimLocs-numLocs

for(i in 1:numNonEof){
  sigmaZMat=sigmaZMat+allSMEigs[numLocs+i]*(leftOverPhiMat[,i]%*%t(leftOverPhiMat[,i]))
  print(i)
}

cSigmaZ=.01

sigmaZMat=sigmaZMat+diag(cSigmaZ,numFullDimLocs)

sigmaZInverse=solve(sigmaZMat)

muAlphaVec=matrix(scaleFactorMean,numLocs,1)

############################################
##### Set MCMC Parameters  #################
############################################
# num.Its=100000
# num.Its=50000
# num.Its=35000
num.Its=5000
# burnInPerc=.10
# burnIn=num.Its*burnInPerc
burnIn=1000

numQuadIntTerms=(curNhTw0*(curNhTw0+1) )/2 
# numVPars=2*curNhTw0+1
# numVPars=2*curNhTw0

###############################
####### Find Quad Index #######
###############################
rowCount=1
colCount=1
indexMatCrtH=rep(NA,numQuadIntTerms)

for(i in 1:curNhTw0){
  
  for(j in 1:curNhTw0){
    
    if(i<=j){
      indexMatCrtH[rowCount]=colCount
      rowCount=rowCount+1
    }
    colCount=colCount+1
  }
  
}

###########################################
############### Forecast Objects ##########
###########################################

outMSEVec=matrix(NA,num.Its,numLocs)
outMSEVec[1,]=1000000

thinSegLen=1
thinSeq=seq(1,num.Its,by=thinSegLen)
allPostForecasts=array(0,c(length(thinSeq),numFullDimLocs,testLen))
forcastCounter=1

#########posterior forecasts Out-Sample
postForOutSample=array(0,c(num.Its,testLen,numLocs))
outMSEVec=rep(num.Its)
outMSEVec[1]=1000000


if(featureLinkInd){
  if(quadInd){
    linCoefSeq=1:curNhTw0
    quadCoefSeq=(curNhTw0+1):(2*curNhTw0)
    featCoefSeq=(2*curNhTw0+1):numVPars
  }else{
    linCoefSeq=1:curNhTw0
    featCoefSeq=(curNhTw0+1):numVPars
  }
}else{
  
  if(quadInd){
    linCoefSeq=1:curNhTw0
    quadCoefSeq=(curNhTw0+1):(2*curNhTw0)
    
  }else{
    linCoefSeq=1:curNhTw0
  }
  
}

##################################
##### Set Priors  ################
##################################
######mu
priorMuSigma2=100
################
### V Priors ###
################
if(numLayers==2){
  sigmaV0=5
  sigmaV1=.001
  
  
  sigmaV0Feat=5
  # sigmaV0Feat=2
  sigmaV1Feat=.001
  
}

if(numLayers==3){
  # sigmaV0=3
  # sigmaV1=.001
  # 
  # 
  # sigmaV0Feat=3
  # sigmaV1Feat=.001
  
  sigmaV0=4
  sigmaV1=.001
  
  
  sigmaV0Feat=4
  sigmaV1Feat=.001
  
}

# sigmaV0Quad=4
# sigmaV1Quad=.01

# sigmaV1Vec=1/c(rep(sigmaV1,curNhTw0))
# sigmaV1Vec=1/c(rep(sigmaV1,curNhTw0),rep(sigmaV1Quad,curNhTw0))
sigmaV1Vec=1/c(rep(sigmaV1,curNhTw0),rep(sigmaV1Feat,layerOneDimRed))


curGammaOneVals=rep(NA,numVPars)
curGammaZeroVals=rep(NA,numVPars)

##### V parameters 
vMatSample=matrix(NA,numLocs,numVPars*ensembleLen)

######################################
####### Get initial V values #########
######################################
for(i in 1:ensembleLen){
  curHMat=matrix(t(hArrayInSamp[,(1+(i-1)*(numVPars)):(i*(numVPars))]),nrow=numVPars,ncol=lenInSampleEmb)
  vMatSample[,(1+(i-1)*(numVPars)):(i*(numVPars))]=t(solve(curHMat%*%t(curHMat)+diag(.001,numVPars))%*%curHMat%*%yScaledTrain)
}


##### Gamma V
gammaVMat=array(NA,c(ensembleLen,numVPars,numLocs))
countVMat=array(0,c(ensembleLen,numVPars,numLocs))



###### piV
if(numLayers==2){
  piV=.25
  # piVQuad=.5
  piVFeat=.25
}

if(numLayers==3){
  piV=.10
  # piVQuad=.5
  piVFeat=.10
  
}

gammaVStrtLin=rbinom(curNhTw0*numLocs*ensembleLen,1,piV)
gammaVMat[,linCoefSeq,]=array(gammaVStrtLin,c(ensembleLen,curNhTw0,numLocs))

# gammaVStrtQuad=rbinom(curNhTw0*numLocs,1,piVQuad)
# gammaVMat[,quadCoefSeq,]=array(gammaVStrtQuad,c(ensembleLen,curNhTw0,numLocs))

gammaVStrtFeat=rbinom(layerOneDimRed*numLocs,1,piVFeat)
gammaVMat[,featCoefSeq,]=array(gammaVStrtFeat,c(ensembleLen,layerOneDimRed,numLocs))

oneIndexGammaV=rep(NA,numVPars)


#####Sigma Epsilon
sigmaEps=rep(NA,num.Its)
sigmaEps[1]=.85


alphaEps=1
betaEps=1


# yTimesNumEnsem=ensembleLen*yScaledTrain


#######################################
####### Basis Coefficient Parameters ##
#######################################
inSampAlphaMat=matrix(NA,nrow(yScaledTrain),numLocs)

outSampAlphaMat=matrix(NA,testLen,numLocs)


###########For plotting
tempAlphaPlotVec=matrix(NA,num.Its,numLocs)

phiTimesPhi=(scaleFactor^2)*t(phiSM)%*%sigmaZInverse%*%phiSM
phiTimesMu=phiSM%*%muAlphaVec
tempPhiTimesZ=scaleFactor*t(phiSM)%*%sigmaZInverse

reducedYTrainFullDim=yTrainFullDim[(curM*tauEmb+1):trainLen,]

phiTimesZ=sapply(1:nrow(yScaledTrain),function(funVar)  tempPhiTimesZ%*%(reducedYTrainFullDim[funVar,]-phiTimesMu) )


#######################################
####### Get Initial In-Samp Preds  ####
#######################################
inSampForMat=matrix(0,lenInSampleEmb,numLocs)
for(j in 1:ensembleLen){
  curEnsemHMat=matrix(t(hArrayInSamp[,(1+(j-1)*(numVPars)):(j*(numVPars))]),nrow=numVPars,ncol=lenInSampleEmb)
  for(i in 1:numLocs){
    curIndexV=(1+(j-1)*(numVPars)):(j*(numVPars))
    inSampForMat[,i]=inSampForMat[,i]+vMatSample[i,curIndexV]%*%curEnsemHMat
  }
  
}

inSampForMat=inSampForMat/ensembleLen

set.seed(887892)

for(it.Count in 2:num.Its){
  
  
  ############################################################
  ######### Sample In-sample Basis Coefs #####################
  ############################################################
  sigmaEpsDiagMat=diag(1/sigmaEps[it.Count-1],numLocs)
  curInSampBasisMean=inSampForMat
  
  for(i in 1:nrow(yScaledTrain)){
    AvalAlpha=phiTimesPhi+sigmaEpsDiagMat
    BvalAlpha=phiTimesZ[,i]+curInSampBasisMean[i,]/sigmaEps[it.Count-1]
    inverseAVal=solve(AvalAlpha)
    
    inSampAlphaMat[i,]= rmvn(1,BvalAlpha%*%inverseAVal,inverseAVal,isChol = F)
    
  }
  
  ####inSampPreds
  inSampForMat=matrix(0,lenInSampleEmb,numLocs)
  outSampForMat=matrix(0,testLen,numLocs)
  
  yTimesNumEnsem=ensembleLen*inSampAlphaMat
  
  #####################################################
  ######### Sample V Parameters  ######################
  #####################################################
  
  
  for(j in 1:ensembleLen){
    
    curEnsemHMat=matrix(t(hArrayInSamp[,(1+(j-1)*(numVPars)):(j*(numVPars))]),nrow=numVPars,ncol=lenInSampleEmb)
    hPrimeMat=curEnsemHMat%*%t(curEnsemHMat)
    
    curOutSampHMat=matrix(t(hArrayOutSamp[,(1+(j-1)*(numVPars)):(j*(numVPars))]),nrow=numVPars,ncol=testLen)
    
    curBtilVTemp=multiDBayesEESNHelper( numVPars, vMatSample, hArrayInSamp,  (1:ensembleLen)[-j] )
    
    
    curIndexV=(1+(j-1)*(numVPars)):(j*(numVPars))
    
    for(i in 1:numLocs){
      
      # BtilVTemp=(curEnsemHMat%*%(inSampAlphaMat[,i]- rowSums( sapply((1:ensembleLen)[-j],function(funVar) t(curVPars[(1+(funVar-1)*(numVPars)):(funVar*(numVPars))]%*%matrix(t(hArrayInSamp[,(1+(funVar-1)*(numVPars)):(funVar*(numVPars))]),nrow=numVPars,ncol=lenInSampleEmb)) ))/ensembleLen) )/sigmaEps[it.Count-1]
      BtilVTemp=(curEnsemHMat%*%(yTimesNumEnsem[,i]-curBtilVTemp[,i]) )/(sigmaEps[it.Count-1]*ensembleLen^2)
      
      
      oneIndexGammaV=which(gammaVMat[j,,i]==1)
      
      sigmaGamMatTemp=diag(sigmaV1Vec,numVPars)
      
      
      diag(sigmaGamMatTemp)[oneIndexGammaV[which(oneIndexGammaV<=curNhTw0)]]=1/sigmaV0
      diag(sigmaGamMatTemp)[oneIndexGammaV[-which(oneIndexGammaV<=curNhTw0)]]=1/sigmaV0Feat
      
      sigmaGamMat=diag(diag(sigmaGamMatTemp))
      
      meanV=hPrimeMat/(sigmaEps[it.Count-1]*ensembleLen^2)+sigmaGamMat
      invMeanV=chol2inv(chol(meanV))
      meanPostV=invMeanV%*%BtilVTemp
      
      vMatSample[i,curIndexV]=as.vector(rmvn(1,meanPostV,chol(invMeanV),isChol=T))
      
      
      #################################################
      ######### Sample Gamma V Parameters  ############
      #################################################
      
      curGammaOneVals[linCoefSeq]=piV*dnorm(vMatSample[i,curIndexV[linCoefSeq]],0,sqrt(sigmaV0))
      curGammaOneVals[featCoefSeq]=piVFeat*dnorm(vMatSample[i,  curIndexV[featCoefSeq]],0,sqrt(sigmaV0Feat))
      
      curGammaZeroVals[linCoefSeq]=(1-piV)*dnorm(vMatSample[i,curIndexV[linCoefSeq]],0,sqrt(sigmaV1))
      curGammaZeroVals[featCoefSeq]=(1-piVFeat)*dnorm(vMatSample[i,  curIndexV[featCoefSeq]],0,sqrt(sigmaV1Feat))
      
      gammaVMat[j,,i]=rbinom(numVPars,1,curGammaOneVals/(curGammaOneVals+curGammaZeroVals))
      countVMat[j,,i]=countVMat[j,,i]+gammaVMat[j,,i]
      
    }
    inSampForMat=inSampForMat+t(vMatSample[,curIndexV]%*%curEnsemHMat)
    outSampForMat=outSampForMat+t(vMatSample[,curIndexV]%*%curOutSampHMat)
    
  }
  
  
  
  inSampForMat=inSampForMat/ensembleLen
  curOutSampFor=outSampForMat/ensembleLen
  
  ##################################################### 
  ######### Sample Sigma_epsilon ###################### 
  #####################################################
  alphaTilEps=(lenInSampleEmb*numLocs)/2+alphaEps
  
  sumDiff=sum((inSampForMat-inSampAlphaMat)^2)
  betaTilEps=.5*(sumDiff)+betaEps
  sigmaEps[it.Count]=1/rgamma(1,alphaTilEps,betaTilEps)
  
  
  
  #####################################################
  ######### Post-Processing  ##########################
  #####################################################
  # tempAllForMat=sapply(1:testLen,function(funVar) phiSM%*%as.matrix(rnorm(numLocs, curOutSampFor[funVar,],sqrt(sigmaEps[it.Count]) ) *scaleFactor +scaleFactorMean) )
  
  
  ###############################################
  ### Calculate Out-Sample Basis Coeffs  ########
  ###############################################
  sigmaEpsDiagMat=diag(1/sigmaEps[it.Count],numLocs)
  
  tempAllForMat=matrix(NA,numFullDimLocs,testLen)
  
  
  for(i in 1:testLen){
    AvalAlphaOut=sigmaEpsDiagMat
    BvalAlphaOut=curOutSampFor[i,]/sigmaEps[it.Count]
    inverseAvalAlphaOut=solve(AvalAlphaOut)
    
    outSampAlphaMat[i,]= rmvn(1,BvalAlphaOut%*%inverseAvalAlphaOut,inverseAvalAlphaOut,isChol = F)
    
    # tempAllForMat[,i]=phiSM%*%(outSampAlphaMat[i,]*scaleFactor+scaleFactorMean)
    tempAllForMat[,i]=rmvn(1,phiSM%*%(outSampAlphaMat[i,]*scaleFactor+scaleFactorMean),sigmaZMat,isChol = F)
  }
  
  
  tempAllForMSE=mean((tempAllForMat[,curMayIndex]-t(yTestALLLocs)[,curMayIndex])^2)
  
  outMSEVec[it.Count]=tempAllForMSE
  
  cat("It ", it.Count, " sigmaEps ", sigmaEps[it.Count], " MSE ", tempAllForMSE, "\n" )
  
  if((it.Count%%thinSegLen)==0){
    # allPostForecasts[forcastCounter,,]=sapply(1:testLen,function(funVar) t(phiPrecip)%*%(rnorm(numLocs, temppostForOutSample[funVar,],sqrt(sigmaEps[it.Count]) ) *scaleFactor +scaleFactorMean) )
    allPostForecasts[forcastCounter,,] =tempAllForMat
    forcastCounter=forcastCounter+1
  }
  
}

save.image("/Users/PatrickMcDermott/Desktop/Two_Layer_SM_B_DEESN_Environ_8_28_2018.RData")






#############################################
#############  Thin the Data  ###############
#############################################
postThinSeq=min(which(thinSeq>burnIn)):max(which(thinSeq<num.Its))
postBurnForcasts=allPostForecasts[postThinSeq,,]
lwCI=.025



#########################################################
###### Posterior Summaries For All Data Types ###########
#########################################################


postMeanAll=matrix(NA,numFullDimLocs,length(curMayIndex))
postSDAll=matrix(NA,numFullDimLocs,length(curMayIndex))
sumCRPSALL=0


for(i in 1:length(curMayIndex)){
  
  postMeanAll[,i]=apply(postBurnForcasts[,,curMayIndex[i]],2,mean)
  postSDAll[,i]=apply(postBurnForcasts[,,curMayIndex[i]],2,sd)
  curCRPS=crps(yTestALLLocs[curMayIndex[i],],cbind(postMeanAll[,i],postSDAll[,i]))$CRPS
  sumCRPSALL=sumCRPSALL+curCRPS
  print(i)
}


allPostForMSE=mean((postMeanAll-t(yTestALLLocs)[,curMayIndex])^2)
cat("All Locs MSE", allPostForMSE , "\n" )

cat("All Locs CRPS", sumCRPSALL , "\n" )




############################################
###### Add back on Anomalies ###############
############################################
posSupportBurnForMat=array(NA,c(dim(postBurnForcasts)[1],dim(postBurnForcasts)[2],dim(postBurnForcasts)[3]))

for(i in 1:dim(postBurnForcasts)[1]){
  posSupportBurnForMat[i,,]=postBurnForcasts[i,,]  +t(outSampMonMeans)
}

posSupportSd=matrix(NA,numFullDimLocs,testLen)
scaledSupportSd=matrix(NA,numFullDimLocs,testLen)
lwPosSupport=matrix(NA,numFullDimLocs,testLen)
upPosSupport=matrix(NA,numFullDimLocs,testLen)


for(i in 1:testLen){
  posSupportSd[,i]=apply(posSupportBurnForMat[,,i],2,sd)
  scaledSupportSd[,i]=apply(t(scale(t(postBurnForcasts[,,i]))),2,sd)
  lwPosSupport[,i]=apply(posSupportBurnForMat[,,i],2,function(x)quantile(x,lwCI))
  upPosSupport[,i]=apply(posSupportBurnForMat[,,i],2,function(x)quantile(x,(1-lwCI)))
  print(i)
}



###################################################
############ SM Spatial Plots #####################
###################################################
####Width:400 Height:700
par(mfrow=c(3,1))
par(mar=c(4.5,6.5,3.5,2.1))
plotPeriod=6


###################################################
############ Anomaly Spatial Plots ################
###################################################
tempAnomSDPlotMax=rep(NA,length(curMayIndex))
tempAnomSDPlotMin=rep(NA,length(curMayIndex))
for(i in 1:length(curMayIndex)){
  tempPlotVecSDESN=rep(NA,numPlotLocs)
  tempSDPlotVals=scaledSupportSd[,curMayIndex[i]]
  sdNonOutliers=which(tempSDPlotVals<quantile(tempSDPlotVals,.98))
  tempSDPlotVals[-sdNonOutliers]=NA
  tempPlotVecSDESN[-nanIndex]=tempSDPlotVals
  tempAnomSDPlotMax[i]=max(na.omit(tempSDPlotVals))
  tempAnomSDPlotMin[i]=min(na.omit(tempSDPlotVals))
}


###Observed
plotVecObs=rep(NA,numPlotLocs)
plotVecObs[-nanIndex]=scale(curYTestAllLocs[curMayIndex[plotPeriod],])
# zlim=c(min(plotVecObs[-nanIndex]),max(plotVecObs[-nanIndex]))
zlimTemp=c(-max(plotVecObs[-nanIndex]),max(plotVecObs[-nanIndex]))

obsPlotObj=image.nan(plotVecObs,zlimTemp,tim.colors(n = 64))
image.plot(smLon,smLat,matrix(obsPlotObj$z,nrow=numXLocs,ncol=numYLocs),xlab="Longitude (deg)",ylab="Lattitude (deg)",zlim=obsPlotObj$zlim,col=obsPlotObj$col,main=paste("True Soil Moisture Anomalies ",monthSeq[curMayIndex[plotPeriod]], " ",yearSeq[curMayIndex[plotPeriod]] ))
map("state",add=T)

###Forecast
plotVecForESN=rep(NA,numPlotLocs)
plotVecForESN[-nanIndex]=scale(postMeanAll[,plotPeriod])
# mean((plotVecForESN[-nanIndex]-plotVecObs[-nanIndex])^2)

forPlotObj=image.nan(plotVecForESN,zlimTemp,tim.colors(n = 64))
image.plot(smLon,smLat,matrix(forPlotObj$z,nrow=numXLocs,ncol=numYLocs),zlim=forPlotObj$zlim,col=forPlotObj$col,xlab="Longitude (deg)",ylab="Lattitude (deg)",main=paste("Bayesian D-EESN Mean Forecast Mean (",numLayers ,"  Layers) ; Lead Time: ",tau, " Months" ))
map("state",add=T)

plotVecSDESN=rep(NA,numPlotLocs)
# tempSDPlotVals=posSupportSd[,curMayIndex[plotPeriod]]
tempSDPlotVals=scaledSupportSd[,curMayIndex[plotPeriod]]
# sdNonOutliers=which(tempSDPlotVals<quantile(tempSDPlotVals,.98))
# tempSDPlotVals[-sdNonOutliers]=NA
plotVecSDESN[-nanIndex]=tempSDPlotVals

# zLimAnomSD=c(min(tempAnomSDPlotMin),max(tempAnomSDPlotMax))
zLimAnomSD=c(.2,max(tempAnomSDPlotMax))

sdPlotObj=image.nan(plotVecSDESN,zLimAnomSD,tim.colors(n = 64))
image.plot(smLon,smLat,matrix(sdPlotObj$z,nrow=numXLocs,ncol=numYLocs),xlab="Longitude (deg)",ylab="Lattitude (deg)",main=paste("Bayesian D-EESN  Standard Deviation " ),zlim=sdPlotObj$zlim,col=sdPlotObj$col)
map("state",add=T)

############################################
######    Skill Score ######################
############################################
####Height:703  Width:703
par(mfrow=c(1,1))
par(mar=c(4.5,4.5,3.5,2.1))
skillScoreVar=1-colMeans((t(postMeanAll) - curYTestAllLocs[curMayIndex,] )^2)/colMeans((curYTestAllLocs[curMayIndex,])^2)

plotVecSS=rep(NA,numPlotLocs)
plotVecSS[-nanIndex]=skillScoreVar
zLimSS=c(-1,1)
ssPlotObj=image.nan(plotVecSS,zLimSS,tim.colors(n = 64))
ssPlotObj$col[which(ssPlotObj$col=="black")]="white"

image.plot(smLon,smLat,matrix(ssPlotObj$z,nrow=numXLocs,ncol=numYLocs),xlab="Longitude (deg)",ylab="Lattitude (deg)",main=paste("Bayesian D-EESN: Skill Score (to Climatology)" ),zlim=ssPlotObj$zlim,col=ssPlotObj$col)
map("state",add=T)
cat("Bayesian D-EESN: % Skill Score >0 ", length(which(skillScoreVar>0))/numFullDimLocs , "\n")



#######delta and sigmaEps
par(mfrow=c(2,2))
par(mar=c(2.1,2.1,2.5,2.1))
plot(sigmaEps[burnIn:num.Its],type='l',main=expression(paste("Trace Plot : ",sigma[epsilon]^2)))


sigmaRunMean=sapply(burnIn:num.Its,function(funVar)mean(sigmaEps[burnIn:funVar]) )
plot(sigmaRunMean,type='l',main=expression(paste("Running Mean : ",sigma[epsilon]^2)))


plot(outMSEVec[burnIn:num.Its],type='l',main="Out-Sample MSE ")


outSampMSERunMean=sapply(burnIn:num.Its,function(funVar)mean(outMSEVec[burnIn:funVar]) )
plot(outSampMSERunMean,type='l',main="Running Mean: Out-Sample MSE ")












############################################
#### Plot All Hidden Units #################
############################################
# for(j in 1:ensembleLen){
j=25  

par(mfrow=c(3,3))
par(mar=c(4.5,4.5,3.5,2.1))
curEnsemHMat=matrix(t(hArrayInSamp[,(1+(j-1)*(numVPars)):(j*(numVPars))]),nrow=numVPars,ncol=lenInSampleEmb)
countLab=1
for(i in 1:curNhTw0){
  if(sum(abs(curEnsemHMat[i,]))>0){
    plot(curEnsemHMat[i,],type='l',lwd=2,xlab="Time Period",ylab="State",main=paste("Hidden-Unit ", countLab,"; Layer 1; Ensemble",j))
    countLab=countLab+1
  }
}

par(mfrow=c(3,2))
par(mar=c(4.5,4.5,3.5,2.1))
tempInSampFLs=matrix(t(plotHOneArrayInSamp[,(1+(j-1)*(layerOneDimRed)):(j*(layerOneDimRed))]),nrow=layerOneDimRed,ncol=lenInSampleEmb)
countLab=1
for(i in 1:layerOneDimRed){
  if(sum(abs(tempInSampFLs[i,]))>0){
    plot(tempInSampFLs[i,],type='l',lwd=2,xlab="Time Period",ylab="State",main=paste("Hidden-Unit ", countLab,"; Dim. Reduced Layer 2; Ensemble",j))
    countLab=countLab+1
  }
}


par(mfrow=c(3,3))
par(mar=c(4.5,4.5,3.5,2.1))
tempInLayerTwo=matrix(t(plotHOneArrayInSamp[,(1+(j-1)*(curNhOne)):(j*(curNhOne))]),nrow=curNhOne,ncol=lenInSampleEmb)
countLab=1
for(i in 1:25){
  if(sum(abs(tempInLayerTwo[i,]))>0){
    plot(tempInLayerTwo[i,],type='l',lwd=2,xlab="Time Period",ylab="State",main=paste("Hidden-Unit ", countLab,"; Layer 2; Ensemble",j))
    countLab=countLab+1
  }
}

# }