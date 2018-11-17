################### File of Functions 
################### libraries
library(Rcpp)
library(compiler)
library(RcppArmadillo)
library(inline)
library(RcppEigen)
library(RcppTN)


####################################################################
############# Read-in Ensemble ESN Functions #######################
####################################################################
sourceCpp("Functions/ensembleESNCPPFuncts_001.cpp")


########################################
##### createLayerOneInput ##############
########################################
createLayerOneInput=function(firstMat,secondMat){
  cbind(firstMat,secondMat)
}
########################################
##### Simulate Lorenz System ###########
########################################

simLorenz=function(totLonezLen,lorenzBurnIn){
  
  sigmaLonez=10
  r=28
  b=8/3
  x=rep(NA,totLonezLen)
  y=rep(NA,totLonezLen)
  z=rep(NA,totLonezLen)
  dt=.01
  x[1]=12
  y[1]=2
  z[1]=9
  for(i in 2:totLonezLen){
    x[i]=x[i-1]+(sigmaLonez*(y[i-1]-x[i-1]))*dt
    y[i]=y[i-1]+(-x[i-1]*z[i-1]+r*x[i]-y[i-1])*dt
    z[i]=z[i-1]+(x[i-1]*y[i-1]-b*z[i-1])*dt
  }
  x=x[lorenzBurnIn:totLonezLen]
  y=y[lorenzBurnIn:totLonezLen]
  z=z[lorenzBurnIn:totLonezLen]
  
  
  return(list(x=x,y=y,z=z))
}




########################################
##### Generate Precip Deep ESN  ########
########################################
precipDeepESNFor=function(curNhOne,lenInSampleEmb,designMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUOne,curDeltaOne,stratValuesXtempOne,curTestLen,dimRedESN,layerOneDimRed,curNhTwo,sampVecESNTwo,uWidthTwo,piWESNTwo,piUESNTwo,wWidthTwo,nColsUTwo,curDeltaTwo,stratValuesXtempTwo,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,designMatrixOutSample,designMatrixInbetween,curInbetweenLen,featureLinkInd,hmatInSampSeq,reducedYTestIndex){
  ############################################################
  ############# In-Sample Expansion ##########################
  ############################################################
  
  genResObj=genResCPP(curNhOne,lenInSampleEmb,designMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUOne,curDeltaOne,stratValuesXtempOne,FALSE,hMatdimOne)
  
  hMatOneInSamp=genResObj$hMat
  xTempLastExpan=genResObj$xTempLast
  uMatOne=genResObj$uMatESNFun
  wMatScaledOne=genResObj$wMatESNFun
  
  
  #######Get inbetween H's
  uInbetweenProdMat=uMatOne%*%t(designMatrixInbetween)
  createInbetweenHOne=createHMat( curNhOne,  curInbetweenLen, wMatScaledOne,  uInbetweenProdMat, xTempLastExpan,FALSE,hMatdimOne)
  finalExpanH=createInbetweenHOne$xTempLast
  inbetweenHMat=createInbetweenHOne$hMat
  
  uProdMat=uMatOne%*%t(designMatrixOutSample)
  createHMatObjOut=createHMat( curNhOne,  curTestLen, wMatScaledOne,  uProdMat, finalExpanH,FALSE,hMatdimOne)
  outSampExpanHMat=createHMatObjOut$hMat
  
  
  ####################################
  ############# EOF Dim. Red. #######
  ###################################
  
  if(dimRedESN=="EOF"){
    hiddenEOFObj=calcEOFS(hMatOneInSamp,layerOneDimRed,1:lenInSampleEmb)
    
    phiMat=hiddenEOFObj$phiMat
    
    resEOFMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%hMatOneInSamp)
    resEOFInbetweenMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%inbetweenHMat)
    resEOFOutSampMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%outSampExpanHMat)
  }
  
  # if(dimRedESN=="EOF"){
  #   # hiddenEOFObj=calcEOFS(hMatOneInSamp,layerOneDimRed,1:lenInSampleEmb)
  #   #
  #   # phiMat=hiddenEOFObj$phiMat
  # 
  # 
  # 
  #   rawDistMat=as.matrix(dist(hMatOneInSamp))
  # 
  #   # epsilonSeq=seq(.01,.015,by=.0005)
  #   # epsilonSeq=seq(.01,.02,by=.001)
  #   epsilonSeq=seq(6,96,by=3)
  #   epsilonLength=length(epsilonSeq)
  #   epsilonIntSeq=1:epsilonLength
  # 
  #   # for(i in 1:epsilonLength){
  # 
  #   # curEpsLE=18
  #   curEpsLE=5
  #   #curEpsLE=3
  #   ######Nearest Neighbor distance
  #   knnMat=graph.knn(rawDistMat,epsilonSeq[curEpsLE])
  #   adjMat=graph.adj(knnMat)
  # 
  #   ######Kernel distance
  #   # adjMat=graph.heat(rawDistMat,epsilonSeq[curEpsLE])
  # 
  #   laplaceMat=graph.laplacian(adjMat)
  #   laplaceEigenObj=eigen(laplaceMat)
  #   phiMat=t(laplaceEigenObj$vectors[,(curNhOne-layerOneDimRed):(curNhOne-1)])
  # 
  # 
  # 
  #   resEOFMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%hMatOneInSamp)
  #   resEOFInbetweenMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%inbetweenHMat)
  #   resEOFOutSampMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%outSampExpanHMat)
  # }
  
  
  #######################################
  ############# Fourier Dim. Red. #######
  #######################################
  
  if(dimRedESN=="Fourier"){
    fullPhiMat=calcFourierAlpha(curNhOne)
    fullAlphaMat=t(fullPhiMat)%*%hMatOneInSamp
    
    allBasisEng=findTopEnergy(fullAlphaMat,curNhOne)
    
    tempTopEngIndex=(sort(allBasisEng,index.return=T,decreasing=T)$ix)[1:(layerOneDimRed/2)]
    
    phiMat=fullPhiMat[,sort(c(tempTopEngIndex*2+1,tempTopEngIndex*2+2))]
    
    resEOFMat =t(t(phiMat)%*%hMatOneInSamp)
    resEOFInbetweenMat=t(t(phiMat)%*%inbetweenHMat)
    resEOFOutSampMat=t(t(phiMat)%*%outSampExpanHMat)
    
  }
  
  
  ################################################
  #### scale by hidden unit ######################
  ################################################
  
  
  if(scaleByHInd){
    eofDesignMatrix=cbind(rep(1,lenInSampleEmb),scale(resEOFMat))
    
    inSampMeans=colMeans(resEOFMat)
    inSampSds=colSds(resEOFMat)
    
    combinInbetween=matrix(NA,curInbetweenLen,layerOneDimRed)
    combinDataOutSamp=matrix(NA,curTestLen,layerOneDimRed)
    for(i in 1:layerOneDimRed){
      combinInbetween[,i]=(resEOFInbetweenMat[,i]-inSampMeans[i])/inSampSds[i]
      combinDataOutSamp[,i]=(resEOFOutSampMat[,i]-inSampMeans[i])/inSampSds[i]
    }
    
  }else{
    ################################################
    #### scale overall #############################
    ################################################
    inSampMeans=mean(resEOFMat)
    inSampSds=sd(resEOFMat)
    
    eofDesignMatrix=cbind(rep(1,lenInSampleEmb),(resEOFMat-inSampMeans)/inSampSds)
    
    combinInbetween=matrix(NA,curInbetweenLen,layerOneDimRed)
    combinDataOutSamp=matrix(NA,curTestLen,layerOneDimRed)
    for(i in 1:layerOneDimRed){
      combinInbetween[,i]=(resEOFInbetweenMat[,i]-inSampMeans)/inSampSds
      combinDataOutSamp[,i]=(resEOFOutSampMat[,i]-inSampMeans)/inSampSds
    }
    
  }
  
  
  eofInbetweenDesignMat=cbind(rep(1,curInbetweenLen),combinInbetween)
  eofDesignMatrixOut=cbind(rep(1,curTestLen),combinDataOutSamp)
  
  
  
  if(quadInd){
    hMatdimTwo=2*curNhTwo
  }else{
    hMatdimTwo=curNhTwo
  }
  
  reservoirObj=genResCPP(curNhTwo,lenInSampleEmb,eofDesignMatrix,sampVecESNTwo,wWidthTwo,piWESNTwo,piUESNTwo,uWidthTwo,nColsUTwo,curDeltaTwo,stratValuesXtempTwo,quadInd,hMatdimTwo)
  
  inSampHMat=reservoirObj$hMat[,hmatInSampSeq]
  
  finDesignMat=inSampHMat
  
  if(featureLinkInd){
    finDesignMat=rbind(inSampHMat,tanh(t( eofDesignMatrix[,2:ncol(eofDesignMatrix)] ))[,hmatInSampSeq] )
  }
  
  # Ridge Regression to get vMat
  vMatESNFL=t(altYInSamp)%*%t(finDesignMat)%*%solve(finDesignMat%*%t(finDesignMat)+ridgeMat)
  uTempESN=reservoirObj$uMatESNFun
  
  # inSampMSE=mean((t(hMat)%*%vMatESNFL-altYInSamp)^2)
  # print(inSampMSE)
  # print(mean((t(altYInSamp)- vMatESNFL%*%finDesignMat)^2))
  
  # print(median(abs(inSampHMat)))
  # print(max(vMatESNFL))
  
  ######inbetween
  uTwoInbetweenProdMat=uTempESN%*%t(eofInbetweenDesignMat)
  inbetwenHMatTwo= createHMat( curNhTwo,  curInbetweenLen, reservoirObj$wMatESNFun,  uTwoInbetweenProdMat,reservoirObj$xTempLast,quadInd,hMatdimTwo)
  finalHMatTwoVec=inbetwenHMatTwo$xTempLast
  
  
  ######out-of-sample
  uProdMatOutSamp=uTempESN%*%t(eofDesignMatrixOut)
  outSampHMatObj= createHMat( curNhTwo,  curTestLen, reservoirObj$wMatESNFun,  uProdMatOutSamp,finalHMatTwoVec,quadInd,hMatdimTwo)
  
  # outSampHMat=outSampHMatObj$hMat[,reducedYTestIndex-tau]
  outSampHMat=outSampHMatObj$hMat[,reducedYTestIndex]
  
  
  finOutSampDesignMat=outSampHMat
  if(featureLinkInd){
    finOutSampDesignMat=rbind(outSampHMat,tanh(t(combinDataOutSamp)[,reducedYTestIndex-tau]))
  }
  # dim(finOutSampDesignMat)
  # dim(vMatESNFL)
  
  scaledOutSampleFor=vMatESNFL%*%finOutSampDesignMat
  
  
  return(list(outSampFor=t(scaledOutSampleFor),inSampHMat=inSampHMat,outSampHMat=outSampHMat))
}


########################################
##### Generate Deep ESN  ###############
########################################

deepESNFor=function(curNhOne,lenInSampleEmb,designMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUOne,curDeltaOne,stratValuesXtempOne,curTestLen,dimRedESN,layerOneDimRed,curNhTwo,sampVecESNTwo,uWidthTwo,piWESNTwo,piUESNTwo,wWidthTwo,nColsUTwo,curDeltaTwo,stratValuesXtempTwo,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,designMatrixOutSample,designMatrixInbetween,curInbetweenLen,featureLinkInd){
  ############################################################
  ############# In-Sample Expansion ##########################
  ############################################################
  
  genResObj=genResCPP(curNhOne,lenInSampleEmb,designMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUOne,curDeltaOne,stratValuesXtempOne,FALSE,hMatdimOne)
  
  hMatOneInSamp=genResObj$hMat
  xTempLastExpan=genResObj$xTempLast
  uMatOne=genResObj$uMatESNFun
  wMatScaledOne=genResObj$wMatESNFun
  
  
  #######Get inbetween H's
  uInbetweenProdMat=uMatOne%*%t(designMatrixInbetween)
  createInbetweenHOne=createHMat( curNhOne,  curInbetweenLen, wMatScaledOne,  uInbetweenProdMat, xTempLastExpan,FALSE,hMatdimOne)
  finalExpanH=createInbetweenHOne$xTempLast
  inbetweenHMat=createInbetweenHOne$hMat
  
  uProdMat=uMatOne%*%t(designMatrixOutSample)
  createHMatObjOut=createHMat( curNhOne,  curTestLen, wMatScaledOne,  uProdMat, finalExpanH,FALSE,hMatdimOne)
  outSampExpanHMat=createHMatObjOut$hMat
  
  
  ####################################
  ############# EOF Dim. Red. #######
  ###################################
  
  if(dimRedESN=="EOF"){
    hiddenEOFObj=calcEOFS(hMatOneInSamp,layerOneDimRed,1:lenInSampleEmb)
    
    phiMat=hiddenEOFObj$phiMat
    
    resEOFMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%hMatOneInSamp)
    resEOFInbetweenMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%inbetweenHMat)
    resEOFOutSampMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%outSampExpanHMat)
  }
  
  
  #######################################
  ############# Fourier Dim. Red. #######
  #######################################
  
  if(dimRedESN=="Fourier"){
    fullPhiMat=calcFourierAlpha(curNhOne)
    fullAlphaMat=t(fullPhiMat)%*%hMatOneInSamp
    
    allBasisEng=findTopEnergy(fullAlphaMat,curNhOne)
    # par(mfrow=c(1,1))
    # plot(allBasisEng,type='l',main="Engery for each Basis")
    
    tempTopEngIndex=(sort(allBasisEng,index.return=T,decreasing=T)$ix)[1:(layerOneDimRed/2)]
    
    phiMat=fullPhiMat[,sort(c(tempTopEngIndex*2+1,tempTopEngIndex*2+2))]
    
    resEOFMat =t(t(phiMat)%*%hMatOneInSamp)
    resEOFInbetweenMat=t(t(phiMat)%*%inbetweenHMat)
    resEOFOutSampMat=t(t(phiMat)%*%outSampExpanHMat)
  }
  
  
  ################################################
  #### scale by hidden unit ######################
  ################################################
  
  
  if(scaleByHInd){
    eofDesignMatrix=cbind(rep(1,lenInSampleEmb),scale(resEOFMat))
    
    inSampMeans=colMeans(resEOFMat)
    inSampSds=colSds(resEOFMat)
    
    combinInbetween=matrix(NA,curInbetweenLen,layerOneDimRed)
    combinDataOutSamp=matrix(NA,curTestLen,layerOneDimRed)
    for(i in 1:layerOneDimRed){
      combinInbetween[,i]=(resEOFInbetweenMat[,i]-inSampMeans[i])/inSampSds[i]
      combinDataOutSamp[,i]=(resEOFOutSampMat[,i]-inSampMeans[i])/inSampSds[i]
    }
    
  }else{
    ################################################
    #### scale overall #############################
    ################################################
    inSampMeans=mean(resEOFMat)
    inSampSds=sd(resEOFMat)
    
    eofDesignMatrix=cbind(rep(1,lenInSampleEmb),(resEOFMat-inSampMeans)/inSampSds)
    
    combinInbetween=matrix(NA,curInbetweenLen,layerOneDimRed)
    combinDataOutSamp=matrix(NA,curTestLen,layerOneDimRed)
    for(i in 1:layerOneDimRed){
      combinInbetween[,i]=(resEOFInbetweenMat[,i]-inSampMeans)/inSampSds
      combinDataOutSamp[,i]=(resEOFOutSampMat[,i]-inSampMeans)/inSampSds
    }
    
  }
  
  
  eofInbetweenDesignMat=cbind(rep(1,curInbetweenLen),combinInbetween)
  eofDesignMatrixOut=cbind(rep(1,curTestLen),combinDataOutSamp)
  
  if(featureLinkInd){
    
    
    if(quadInd){
      hMatdimTwo=2*curNhTwo
    }else{
      hMatdimTwo=curNhTwo
    }
    
    reservoirObj=genResCPP(curNhTwo,lenInSampleEmb,eofDesignMatrix,sampVecESNTwo,wWidthTwo,piWESNTwo,piUESNTwo,uWidthTwo,nColsUTwo,curDeltaTwo,stratValuesXtempTwo,quadInd,hMatdimTwo)
    
    inSampHMat=rbind(reservoirObj$hMat,tanh(t( eofDesignMatrix[,2:ncol(eofDesignMatrix)] )))
    
    
    
    # Ridge Regression to get vMat
    vMatESNFL=t(altYInSamp)%*%t(inSampHMat)%*%solve(inSampHMat%*%t(inSampHMat)+ridgeMat)
    uTempESN=reservoirObj$uMatESNFun
    
    # print(mean((t(altYInSamp)- vMatESNFL%*%inSampHMat)^2))
    
    ######inbetween
    uTwoInbetweenProdMat=uTempESN%*%t(eofInbetweenDesignMat)
    inbetwenHMatTwo= createHMat( curNhTwo,  curInbetweenLen, reservoirObj$wMatESNFun,  uTwoInbetweenProdMat,reservoirObj$xTempLast,quadInd,hMatdimTwo)
    finalHMatTwoVec=inbetwenHMatTwo$xTempLast
    
    
    ######out-of-sample
    uProdMatOutSamp=uTempESN%*%t(eofDesignMatrixOut)
    outSampHMatObj= createHMat( curNhTwo,  curTestLen, reservoirObj$wMatESNFun,  uProdMatOutSamp,finalHMatTwoVec,quadInd,hMatdimTwo)
    
    outSampHMat=rbind(outSampHMatObj$hMat,tanh(t(combinDataOutSamp)))
    scaledOutSampleFor=altScaleFactor*(vMatESNFL%*%outSampHMat)+addScaleMat
    
  }else{
    esnForObj=multiObjCalcESNForecastsSpinForward(curNhTwo,lenInSampleEmb,curTestLen,eofDesignMatrix,eofDesignMatrixOut,sampVecESNTwo,uWidthTwo,piWESNTwo,piUESNTwo,wWidthTwo,nColsUTwo,curDeltaTwo,stratValuesXtempTwo,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,quadInd,eofInbetweenDesignMat,curInbetweenLen)
    scaledOutSampleFor=esnForObj$outSampFor
    inSampHMat=esnForObj$hMat
    outSampHMat=esnForObj$hMat
  }
  
  return(list(outSampFor=t(scaledOutSampleFor),inSampHMat=inSampHMat,outSampHMat=outSampHMat))
}



####################################################
##### Generate Multi Layer Deep ESN  ###############
####################################################
#####For ESN that has two or more layers

multiLayerDeepESNFor=function(curNhOne,lenInSampleEmb,designMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUOne,deltaVec,stratValuesXtempOne,curTestLen,dimRedESN,dimRedDim,curNhTwo,sampVecESNTwo,uWidthTwo,piWESNTwo,piUESNTwo,wWidthTwo,nColsUTwo,stratValuesXtempTwo,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,designMatrixOutSample,designMatrixInbetween,curInbetweenLen,featureLinkInd,numLayers){
  ############################################################
  ############# In-Sample Expansion ##########################
  ############################################################
  
  genResObj=genResCPP(curNhOne,lenInSampleEmb,designMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUOne,deltaVec[1],stratValuesXtempOne,FALSE,hMatdimOne)
  
  hMatOneInSamp=genResObj$hMat
  xTempLastExpan=genResObj$xTempLast
  uMatOne=genResObj$uMatESNFun
  wMatScaledOne=genResObj$wMatESNFun
  
  
  #######Get inbetween H's
  uInbetweenProdMat=uMatOne%*%t(designMatrixInbetween)
  createInbetweenHOne=createHMat( curNhOne,  curInbetweenLen, wMatScaledOne,  uInbetweenProdMat, xTempLastExpan,FALSE,hMatdimOne)
  finalExpanH=createInbetweenHOne$xTempLast
  inbetweenHMat=createInbetweenHOne$hMat
  
  uProdMat=uMatOne%*%t(designMatrixOutSample)
  createHMatObjOut=createHMat( curNhOne,  curTestLen, wMatScaledOne,  uProdMat, finalExpanH,FALSE,hMatdimOne)
  outSampExpanHMat=createHMatObjOut$hMat
  
  ####################################
  ############# EOF Dim. Red. #######
  ###################################
  
  curDimRedObj=dimRedDeepESN(dimRedESN,hMatOneInSamp,inbetweenHMat,outSampExpanHMat,dimRedDim,lenInSampleEmb,curNhOne,curInbetweenLen,curTestLen)
  
  curEOFDesignMatrix=curDimRedObj$eofDesignMatrix
  curEOFInbetweenDesignMat=curDimRedObj$eofInbetweenDesignMat
  curEOFDesignMatrixOut=curDimRedObj$eofDesignMatrixOut
  
  
  
  if(numLayers==2){
    
    eofDesignMatrix=curEOFDesignMatrix
    eofInbetweenDesignMat=curEOFInbetweenDesignMat
    eofDesignMatrixOut=curEOFDesignMatrixOut
    
    inSampleFeatureLinks=t( eofDesignMatrix[,2:ncol(eofDesignMatrix)] )
    outSampFeatureLinks=t(eofDesignMatrixOut[,2:ncol(eofDesignMatrixOut)])
    
    
  }else{
    
    inSampleFeatureLinks=matrix(NA,(numLayers-1)*dimRedDim,lenInSampleEmb)
    outSampFeatureLinks=matrix(NA,(numLayers-1)*dimRedDim,curTestLen)
    
    inSampleFeatureLinks[1:dimRedDim,]=t( curEOFDesignMatrix[,2:ncol(curEOFDesignMatrix)] )
    outSampFeatureLinks[1:dimRedDim,]=t( curEOFDesignMatrixOut[,2:ncol(curEOFDesignMatrixOut)] )
    
    
    for(i in 2:(numLayers-1)){
      
      
      #############################################
      ######### Generate Next Reservior ###########
      #############################################
      
      genResObj=genResCPP(curNhOne,lenInSampleEmb,curEOFDesignMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUTwo,deltaVec[i],stratValuesXtempOne,FALSE,hMatdimOne)
      
      hMatOneInSamp=genResObj$hMat
      xTempLastExpan=genResObj$xTempLast
      uMatOne=genResObj$uMatESNFun
      wMatScaledOne=genResObj$wMatESNFun
      
      
      #######Get inbetween H's
      uInbetweenProdMat=uMatOne%*%t(curEOFInbetweenDesignMat)
      createInbetweenHOne=createHMat( curNhOne,  curInbetweenLen, wMatScaledOne,  uInbetweenProdMat, xTempLastExpan,FALSE,hMatdimOne)
      finalExpanH=createInbetweenHOne$xTempLast
      inbetweenHMat=createInbetweenHOne$hMat
      
      uProdMat=uMatOne%*%t(curEOFDesignMatrixOut)
      createHMatObjOut=createHMat( curNhOne,  curTestLen, wMatScaledOne,  uProdMat, finalExpanH,FALSE,hMatdimOne)
      outSampExpanHMat=createHMatObjOut$hMat
      
      ####################################
      ############# EOF Dim. Red. #######
      ###################################
      
      curDimRedObj=dimRedDeepESN(dimRedESN,hMatOneInSamp,inbetweenHMat,outSampExpanHMat,dimRedDim,lenInSampleEmb,curNhOne,curInbetweenLen,curTestLen)
      
      curEOFDesignMatrix=curDimRedObj$eofDesignMatrix
      curEOFInbetweenDesignMat=curDimRedObj$eofInbetweenDesignMat
      curEOFDesignMatrixOut=curDimRedObj$eofDesignMatrixOut
      
      inSampleFeatureLinks[((i-1)*dimRedDim+1):(i*dimRedDim),]=t( curEOFDesignMatrix[,2:ncol(curEOFDesignMatrix)] )
      outSampFeatureLinks[((i-1)*dimRedDim+1):(i*dimRedDim),]=t( curEOFDesignMatrixOut[,2:ncol(curEOFDesignMatrixOut)] )
      
      
      if(i==((numLayers-1))){
        
        eofDesignMatrix=curEOFDesignMatrix
        eofInbetweenDesignMat=curEOFInbetweenDesignMat
        eofDesignMatrixOut=curEOFDesignMatrixOut
        
      }
      
      
    }
    
    
  }
  
  
  
  if(featureLinkInd){
    
    
    if(quadInd){
      hMatdimTwo=2*curNhTwo
    }else{
      hMatdimTwo=curNhTwo
    }
    
    
    reservoirObj=genResCPP(curNhTwo,lenInSampleEmb,eofDesignMatrix,sampVecESNTwo,wWidthTwo,piWESNTwo,piUESNTwo,uWidthTwo,nColsUTwo,deltaVec[numLayers],stratValuesXtempTwo,quadInd,hMatdimTwo)
    
    # inSampHMat=rbind(inSampHMat,tanh(t( eofDesignMatrix[,2:ncol(eofDesignMatrix)] )))
    inSampHMat=rbind(reservoirObj$hMat,tanh(inSampleFeatureLinks))
    
    
    # Ridge Regression to get vMat
    vMatESNFL=t(altYInSamp)%*%t(inSampHMat)%*%solve(inSampHMat%*%t(inSampHMat)+ridgeMat)
    uTempESN=reservoirObj$uMatESNFun
    
    # print(mean((t(altYInSamp)- vMatESNFL%*%inSampHMat)^2))
    
    ######inbetween
    uTwoInbetweenProdMat=uTempESN%*%t(eofInbetweenDesignMat)
    inbetwenHMatTwo= createHMat( curNhTwo,  curInbetweenLen, reservoirObj$wMatESNFun,  uTwoInbetweenProdMat,reservoirObj$xTempLast,quadInd,hMatdimTwo)
    finalHMatTwoVec=inbetwenHMatTwo$xTempLast
    
    
    ######out-of-sample
    uProdMatOutSamp=uTempESN%*%t(eofDesignMatrixOut)
    outSampHMatObj= createHMat( curNhTwo,  curTestLen, reservoirObj$wMatESNFun,  uProdMatOutSamp,finalHMatTwoVec,quadInd,hMatdimTwo)
    
    outSampHMat=rbind(outSampHMatObj$hMat,tanh(outSampFeatureLinks))
    scaledOutSampleFor=altScaleFactor*(vMatESNFL%*%outSampHMat)+addScaleMat
    
  }else{
    esnForObj=multiObjCalcESNForecastsSpinForward(curNhTwo,lenInSampleEmb,curTestLen,eofDesignMatrix,eofDesignMatrixOut,sampVecESNTwo,uWidthTwo,piWESNTwo,piUESNTwo,wWidthTwo,nColsUTwo,deltaVec[numLayers],stratValuesXtempTwo,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,quadInd,eofInbetweenDesignMat,curInbetweenLen)
    scaledOutSampleFor=esnForObj$outSampFor
    inSampHMat=esnForObj$hMat
    outSampHMat=esnForObj$hMatOutSamp
  }
  
  return(list(outSampFor=t(scaledOutSampleFor),inSampHMat=inSampHMat,outSampHMat=outSampHMat,hMatOneInSamp=hMatOneInSamp,inSampleFeatureLinks=inSampleFeatureLinks,eofDesignMatrix=eofDesignMatrix[,2:(dimRedDim+1)],eofInbetweenDesignMat=eofInbetweenDesignMat[,2:(dimRedDim+1)],eofDesignMatrixOut=eofDesignMatrixOut[,2:(dimRedDim+1)]))
}



####################################################
##### Dim. Red. Deep ESN  ##########################
####################################################


dimRedDeepESN=function(dimRedESN,hMatOneInSamp,inbetweenHMat,outSampExpanHMat,dimRedDim,lenInSampleEmb,curNhOne,curInbetweenLen,curTestLen){
  ####################################
  ############# EOF Dim. Red. #######
  ###################################
  
  if(dimRedESN=="EOF"){
    hiddenEOFObj=calcEOFS(hMatOneInSamp,dimRedDim,1:lenInSampleEmb)
    
    phiMat=hiddenEOFObj$phiMat
    
    resEOFMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%hMatOneInSamp)
    resEOFInbetweenMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%inbetweenHMat)
    resEOFOutSampMat=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%outSampExpanHMat)
  }
  
  
  #######################################
  ############# Fourier Dim. Red. #######
  #######################################
  
  if(dimRedESN=="Fourier"){
    fullPhiMat=calcFourierAlpha(curNhOne)
    fullAlphaMat=t(fullPhiMat)%*%hMatOneInSamp
    
    allBasisEng=findTopEnergy(fullAlphaMat,curNhOne)
    # par(mfrow=c(1,1))
    # plot(allBasisEng,type='l',main="Engery for each Basis")
    
    tempTopEngIndex=(sort(allBasisEng,index.return=T,decreasing=T)$ix)[1:(dimRedDim/2)]
    
    phiMat=fullPhiMat[,sort(c(tempTopEngIndex*2+1,tempTopEngIndex*2+2))]
    
    resEOFMat =t(t(phiMat)%*%hMatOneInSamp)
    resEOFInbetweenMat=t(t(phiMat)%*%inbetweenHMat)
    resEOFOutSampMat=t(t(phiMat)%*%outSampExpanHMat)
  }
  
  
  ################################################
  #### scale by hidden unit ######################
  ################################################
  
  
  if(scaleByHInd){
    eofDesignMatrix=cbind(rep(1,lenInSampleEmb),scale(resEOFMat))
    
    inSampMeans=colMeans(resEOFMat)
    inSampSds=colSds(resEOFMat)
    
    combinInbetween=matrix(NA,curInbetweenLen,dimRedDim)
    combinDataOutSamp=matrix(NA,curTestLen,dimRedDim)
    for(i in 1:dimRedDim){
      combinInbetween[,i]=(resEOFInbetweenMat[,i]-inSampMeans[i])/inSampSds[i]
      combinDataOutSamp[,i]=(resEOFOutSampMat[,i]-inSampMeans[i])/inSampSds[i]
    }
    
  }else{
    ################################################
    #### scale overall #############################
    ################################################
    inSampMeans=mean(resEOFMat)
    inSampSds=sd(resEOFMat)
    
    eofDesignMatrix=cbind(rep(1,lenInSampleEmb),(resEOFMat-inSampMeans)/inSampSds)
    
    combinInbetween=matrix(NA,curInbetweenLen,dimRedDim)
    combinDataOutSamp=matrix(NA,curTestLen,dimRedDim)
    for(i in 1:dimRedDim){
      combinInbetween[,i]=(resEOFInbetweenMat[,i]-inSampMeans)/inSampSds
      combinDataOutSamp[,i]=(resEOFOutSampMat[,i]-inSampMeans)/inSampSds
    }
    
  }
  
  eofInbetweenDesignMat=cbind(rep(1,curInbetweenLen),combinInbetween)
  eofDesignMatrixOut=cbind(rep(1,curTestLen),combinDataOutSamp)
  
  return(list(eofDesignMatrix=eofDesignMatrix,eofInbetweenDesignMat=eofInbetweenDesignMat,eofDesignMatrixOut=eofDesignMatrixOut))
  
}

################################################################
############### Ensembel ESN Set Parameters ####################
################################################################
setParsEESN=function(regPar,nh,numLocs,m){
  inSize = numLocs*(m+1)
  nColsU=inSize+1
  
  #########sampleVec
  sampVecESN=0:(nh-1)
  
  stratValuesXtemp=rep(0,nh)
  
  
  #####create Ridge Matrix
  reg =regPar
  if(quadInd){
    ridgeMat=reg*diag(2*nh )
  }else{
    ridgeMat=reg*diag(nh)
  }
  return(list(stratValuesXtemp=stratValuesXtemp,sampVecESN=sampVecESN,ridgeMat=ridgeMat,nColsU=nColsU))
}

################################################################
################################################################
############### Ensembel ESN Function ##########################
################################################################
################################################################
ensembleESN=function(nh,reg,delta,sparsePerc,uSparsePerc,aWidth,uWidth,leakRate,inSampleX,inSampleY,outSampleX,determinMseLen,trainLenEcho,testEchoLen,numLocsEcho){
  
  determinEsnMse=rep(NA,determinMseLen)
  
  forMatDetminESN=array(NA,c(determinMseLen,testEchoLen,numLocsEcho))
  
  for(i in 1:determinMseLen){
    
    
    reservoirObj=genRes(nh,inSize,sparsePerc,trainLenEcho,delta,inSampleX,aWidth,uWidth,leakRate,uSparsePerc,quadInd)
    
    wMatNBEcho=reservoirObj$W
    uMatNBEcho=reservoirObj$U
    hMatNBEcho=reservoirObj$H
    
    xTempLast=reservoirObj$xTempLast
    
    
    ###############################################
    ########### In-Sample MSE #####################
    ###############################################
    
    
    ######allows for quadratic and interaction terms 
    if(quadInd){
      vEcho = t(inSampleY) %*% t(hMatNBEcho) %*%chol2inv( chol( hMatNBEcho %*% t(hMatNBEcho) + reg*diag(nh+(nh*(nh+1) )/2 ) ) )
    }else{
      vEcho = t(inSampleY) %*% t(hMatNBEcho) %*%chol2inv( chol( hMatNBEcho %*% t(hMatNBEcho) + reg*diag(2*nh ) ) )
    }
    
    
    
    yMatEcho = matrix(0,numLocsEcho,testEchoLen)
    
    xTemp =xTempLast
    for (t in 1:testEchoLen){
      if(is.matrix(outSampleX)){
        u=as.vector(outSampleX[t,])
      }else{
        u=as.vector(outSampleX[t,,])
      }
      
      xTemp =(1-leakRate)*xTemp + tanh( uMatNBEcho %*% c(1,u) + wMatNBEcho %*% xTemp )
      ######allows for quadratic and interaction terms 
      if(quadInd){
        yTemp = vEcho %*% c(xTemp,(xTemp%*%t(xTemp))[upper.tri(xTemp%*%t(xTemp) ,diag=T)])
      }else{
        yTemp = vEcho %*% c(xTemp,xTemp^2)
      }
      
      yMatEcho[,t] = yTemp*altScaleFactor+altYMean
      u = yTemp
    }
    
    forMatDetminESN[i,,]=t(yMatEcho)
    
    
  }
  
  return(forMatDetminESN)
  
}



############################################################
############### Create Embed Data ##########################
############################################################

createEmbedRNNData=function(curTrainLen,m,tauEmb,yTrain,rawDataInput,numLocs,curXTestIndex,curTestLen,curInBetweenIndexes,curLenInbetween){
  lenInSampleEmb=curTrainLen-(m*tauEmb)
  
  altInSampleXRaw=array(NA,c(lenInSampleEmb,m+1,numLocs))
  
  for(i in 1:lenInSampleEmb){
    #####using laplacian eigenmaps
    
    altInSampleXRaw[i,,]=rawDataInput[seq(i,(m*tauEmb+i),by=tauEmb),]
  }
  
  #######Scale in-sample x and y
  altYInSampRaw=yTrain[(m*tauEmb+1):curTrainLen,]
  altYMean=mean(altYInSampRaw)
  altScaleFactor=sd(altYInSampRaw)
  
  altYInSamp=(altYInSampRaw-altYMean)/altScaleFactor
  
  
  meanXTrainMatrix=mean(rawDataInput[1:curTrainLen,])
  sdXTrainMatrix=sd(rawDataInput[1:curTrainLen,])
  
  
  altInSampleX=(altInSampleXRaw-meanXTrainMatrix)/sdXTrainMatrix
  
  
  designMatrix=matrix(1,lenInSampleEmb,(m+1)*numLocs+1)
  for(i in 1:lenInSampleEmb){
    designMatrix[i,2:((m+1)*numLocs+1)]=as.vector(altInSampleX[i,,])
  }
  
  ###############################
  ####### Inbetween #############
  ###############################
  altInbetweenXRaw=array(NA,c(curLenInbetween,m+1,numLocs))
  for(i in 1:curLenInbetween){
    altInbetweenXRaw[i,,]=rawDataInput[seq(curInBetweenIndexes[i]-(m*tauEmb),curInBetweenIndexes[i],by=tauEmb),]
  }
  
  
  #######Scale in-sample x and y
  altInbetweenX=(altInbetweenXRaw-meanXTrainMatrix)/sdXTrainMatrix
  
  inbetweenDesignMatrix=matrix(1,curLenInbetween,(m+1)*numLocs+1)
  for(i in 1:curLenInbetween){
    inbetweenDesignMatrix[i,2:((m+1)*numLocs+1)]=as.vector(altInbetweenX[i,,])
  }
  
  
  ####### Out Sample
  altOutSampleXRaw=array(NA,c(curTestLen,m+1,numLocs))
  for(i in 1:curTestLen){
    altOutSampleXRaw[i,,]=rawDataInput[seq(curXTestIndex[i]-(m*tauEmb),curXTestIndex[i],by=tauEmb),]
  }
  
  
  #######Scale in-sample x and y
  altOutSampleX=(altOutSampleXRaw-meanXTrainMatrix)/sdXTrainMatrix
  
  designMatrixOutSample=matrix(1,curTestLen,(m+1)*numLocs+1)
  for(i in 1:curTestLen){
    designMatrixOutSample[i,2:((m+1)*numLocs+1)]=as.vector(altOutSampleX[i,,])
  }
  
  
  
  ######additive scale matrix
  addScaleMat=matrix(altYMean,numLocs,curTestLen)
  
  
  return(list(curInSampMean=altYMean,curInSampSD=altScaleFactor,curYInSamp=altYInSamp,curInSampX=altInSampleX,curOutSampX=altOutSampleX,lenInSampleEmb=lenInSampleEmb,designMatrix=designMatrix,designMatrixOutSample=designMatrixOutSample,curTestLen=curTestLen,addScaleMat=addScaleMat,inbetweenDesignMatrix=inbetweenDesignMatrix ))
  
  
}




################################################################
####  individualLocScaleCreateEmbedRNNData ##################### 
###################### scales all of the locations of the input
###################### indivually instead over all the locs
################################################################
individualLocScaleCreateEmbedRNNData=function(curTrainLen,m,tauEmb,yTrain,rawDataInput,numLocs,curXTestIndex,curTestLen,curInBetweenIndexes,curLenInbetween){
  lenInSampleEmb=curTrainLen-(m*tauEmb)
  
  altInSampleXRaw=array(NA,c(lenInSampleEmb,m+1,numLocs))
  
  for(i in 1:lenInSampleEmb){
    #####using laplacian eigenmaps
    altInSampleXRaw[i,,]=rawDataInput[seq(i,(m*tauEmb+i),by=tauEmb),]
  }
  
  #######Scale in-sample x and y
  altYInSampRaw=yTrain[(m*tauEmb+1):curTrainLen,]
  altYMean=mean(altYInSampRaw)
  altScaleFactor=sd(altYInSampRaw)
  
  altYInSamp=(altYInSampRaw-altYMean)/altScaleFactor
  
  
  
  meanXTrainMatrix=matrix(NA,numLocs,m+1)
  sdXTrainMatrix=matrix(NA,numLocs,m+1)
  altInSampleX=array(NA,c(lenInSampleEmb,m+1,numLocs))
  
  
  for(i in 1:numLocs){
    for(j in 1:(m+1)){
      
      meanXTrainMatrix[i,j]=mean(altInSampleXRaw[,j,i])
      sdXTrainMatrix[i,j]=sd(altInSampleXRaw[,j,i])
      altInSampleX[,j,i]=(altInSampleXRaw[,j,i]-meanXTrainMatrix[i,j])/sdXTrainMatrix[i,j]
      
    }
  }
  
  
  designMatrix=matrix(1,lenInSampleEmb,(m+1)*numLocs+1)
  for(i in 1:lenInSampleEmb){
    designMatrix[i,2:((m+1)*numLocs+1)]=as.vector(altInSampleX[i,,])
  }
  
  
  ###############################
  ####### Inbetween #############
  ###############################
  altInbetweenXRaw=array(NA,c(curLenInbetween,m+1,numLocs))
  for(i in 1:curLenInbetween){
    altInbetweenXRaw[i,,]=rawDataInput[seq(curInBetweenIndexes[i]-(m*tauEmb),curInBetweenIndexes[i],by=tauEmb),]
  }
  
  
  #######Scale in-sample x and y
  # altInbetweenX=(altInbetweenXRaw-meanXTrainMatrix)/sdXTrainMatrix
  
  altInbetweenX=array(NA,c(curLenInbetween,m+1,numLocs))
  for(i in 1:numLocs){
    for(j in 1:(m+1)){
      altInbetweenX[,j,i]=(altInbetweenXRaw[,j,i]-meanXTrainMatrix[i,j])/sdXTrainMatrix[i,j]
      
    }
  }
  
  
  inbetweenDesignMatrix=matrix(1,curLenInbetween,(m+1)*numLocs+1)
  for(i in 1:curLenInbetween){
    inbetweenDesignMatrix[i,2:((m+1)*numLocs+1)]=as.vector(altInbetweenX[i,,])
  }
  
  
  ####### Out Sample
  altOutSampleXRaw=array(NA,c(curTestLen,m+1,numLocs))
  for(i in 1:curTestLen){
    altOutSampleXRaw[i,,]=rawDataInput[seq(curXTestIndex[i]-(m*tauEmb),curXTestIndex[i],by=tauEmb),]
  }
  
  
  #######Scale in-sample x and y
  # altOutSampleX=(altOutSampleXRaw-meanXTrainMatrix)/sdXTrainMatrix
  
  
  altOutSampleX=array(NA,c(curTestLen,m+1,numLocs))
  
  for(i in 1:numLocs){
    for(j in 1:(m+1)){
      altOutSampleX[,j,i]=(altOutSampleXRaw[,j,i]-meanXTrainMatrix[i,j])/sdXTrainMatrix[i,j]
    }
  }
  
  designMatrixOutSample=matrix(1,curTestLen,(m+1)*numLocs+1)
  for(i in 1:curTestLen){
    designMatrixOutSample[i,2:((m+1)*numLocs+1)]=as.vector(altOutSampleX[i,,])
  }
  
  
  
  ######additive scale matrix
  addScaleMat=matrix(altYMean,numLocs,curTestLen)
  
  
  return(list(curInSampMean=altYMean,curInSampSD=altScaleFactor,curYInSamp=altYInSamp,curInSampX=altInSampleX,curOutSampX=altOutSampleX,lenInSampleEmb=lenInSampleEmb,designMatrix=designMatrix,designMatrixOutSample=designMatrixOutSample,curTestLen=curTestLen,addScaleMat=addScaleMat,inbetweenDesignMatrix=inbetweenDesignMatrix ))
  
  
}
################################################################
############### Create Non-Embed Data ##########################
################################################################

createRNNData=function(curTrainLen,yTrain,rawDataInput,numLocs,curXTestIndex,curTestLen){
  lenInSampleEmb=curTrainLen
  
  
  #######Scale in-sample x and y
  altYInSampRaw=yTrain[1:curTrainLen, ]
  altYMean=mean(altYInSampRaw)
  altScaleFactor=sd(altYInSampRaw)
  
  altYInSamp=(altYInSampRaw-altYMean)/altScaleFactor
  
  altInSampleXRaw=rawDataInput[inputInSamplSeq[1:curTrainLen],]
  
  meanXTrainMatrix=mean(altInSampleXRaw)
  sdXTrainMatrix=sd(altInSampleXRaw)
  
  altInSampleX=(altInSampleXRaw-meanXTrainMatrix)/sdXTrainMatrix
  
  designMatrix=cbind(rep(1,lenInSampleEmb),altInSampleX)
  
  ####### Out Sample
  altOutSampleXRaw=rawDataInput[curXTestIndex,]
  altOutSampleX=(altOutSampleXRaw-meanXTrainMatrix)/sdXTrainMatrix
  designMatrixOutSample=cbind(rep(1,curTestLen),altOutSampleX)
  
  
  ######additive scale matrix
  addScaleMat=matrix(altYMean,numLocs,curTestLen)
  
  return(list(curInSampMean=altYMean,curInSampSD=altScaleFactor,curYInSamp=altYInSamp,curInSampX=altInSampleX,curOutSampX=altOutSampleX,lenInSampleEmb=lenInSampleEmb,designMatrix=designMatrix,designMatrixOutSample=designMatrixOutSample,curTestLen=curTestLen,addScaleMat=addScaleMat ))
  
}


################################################################
################################################################
############### Read in SST Data ###############################
################################################################
################################################################


get.SST.data=function(strtyr,bgyrClimateAvg,endyrClimateAvg,dataSetType){
  if(dataSetType=="Feb2015"){
    n_t=1382
  }
  
  if(dataSetType=="Jan2017"){
    n_t=1405
  }
  
  n_x=84
  n_y=46
  
  
  SSTlonlat=as.matrix(read.table("SSTlonlat.dat",sep=","))
  
  ### this data goes from 1900 up to Feb 2015
  # raw.SST=as.matrix(read.table("SST_Data.tsv"))
  
  ### this data goes from 1900 up to Jan 2017
  raw.SST=as.matrix(read.table("SST_1_2017.tsv"))
  
  
  
  
  sst.long=sort(unique(SSTlonlat[,1]))
  sst.lat=sort(unique(SSTlonlat[,2]))
  
  sst.Data=matrix(NA,(n_x*n_y),n_t)
  for(i in 1:n_t){
    sst.Data[,i]=raw.SST[(((i-1)*n_y+1):(i*n_y)),]
  }
  
  SSTlonlatTemp=matrix(NA,nrow(sst.Data),2)
  
  for(i in 1:length(sst.long)){
    SSTlonlatTemp[(1+(i-1)*length(sst.lat)):(i*length(sst.lat)),]=cbind(rep(sst.long[i],length(sst.lat)),sst.lat )
  }
  
  sst.sea=sst.Data[which(sst.Data[,ncol(sst.Data)]>=0),]
  SSTlonlatTemp=SSTlonlatTemp[which(sst.Data[,ncol(sst.Data)]>=0),]
  
  #   year.seq=(which(seq(1900,2014,1)==strtyr)-1)+seq(1,length(seq(1900,2014,1))-(which(seq(1900,2014,1)==strtyr)-1),1)
  #   year.seq=(which(seq(1900,2014,1)==strtyr)-1)+seq(1,length(seq(1900,2014,1))-(which(seq(1900,2014,1)==strtyr)-1),1)-1
  
  
  # sst.sea=sst.sea[,((strtyr-1900)*12+1):(((strtyr-1900)*12+1)+(length(seq(strtyr,lstYear,1))*12))]
  sst.sea=sst.sea[,((strtyr-1900)*12+1):ncol(sst.sea)]
  
  ##################################################################
  ########################## Calculate Anomolies####################
  ##################################################################
  n_t.sea=ncol(sst.sea)
  n_s=nrow(sst.sea)
  sst.anom=matrix(0,n_s,n_t.sea)
  
  sst.anom.land.sea=matrix(NA,nrow(sst.Data),n_t.sea)
  
  moMN=matrix(0,12,n_s)
  # bgyrClimateAvg=1970
  # endyrClimateAvg=1999
  
  # bgyrClimateAvg=1981
  # endyrClimateAvg=2010
  
  bg.time=(bgyrClimateAvg-strtyr)*12
  end.Time=((endyrClimateAvg-strtyr)+1)*12
  
  for(k in 1:12){
    tavg=seq((bg.time+k),end.Time,12)
    moMN[k,]=rowMeans(sst.sea[,tavg])
  }
  
  for(t in 1:n_t.sea){
    tind.x=t%%12
    if(tind.x==0){
      tind.x=12
    }
    sst.anom.land.sea[which(sst.Data[,ncol(sst.Data)]>=0),t]=sst.sea[,t]-moMN[tind.x,]
    sst.anom[,t]=sst.sea[,t]-moMN[tind.x,]
  }
  
  ###########################################################
  ################ Plot SST for May 1998 ####################
  ###########################################################
  
  zlim=c(-3,3)
  month="May"
  year="1998"
  title="Observed SST"
  
  # plot.SST.Field(sst.anom.land.sea[,length(1960:2011)*12+2-18],zlim,month,year,sst.long,sst.lat,n_x,n_y,title)
  
  year="2013"
  zlim=c(-3.5,3.5)
  
  # plot.SST.Field(sst.anom.land.sea[,length(1970:2012)*12+5],zlim,month,year,sst.long,sst.lat,n_x,n_y,title)
  seaIndex=which(sst.Data[,ncol(sst.Data)]>=0)
  
  
  return(list(sst.anom=sst.anom,sst.anom.land.sea=sst.anom.land.sea,sst.long=sst.long,sst.lat=sst.lat,n_y=n_y,n_x=n_x,seaIndex=seaIndex,SSTlonlatTemp=SSTlonlatTemp))
}

#########################################################
############## Function A Single SST Field
############## NOT TO BE USED TO PRODUCE ANALOG FIGURE
#########################################################
####z is what we want to plot 

plot.SST.Field=function(z,zlim,month,year,sst.long,sst.lat,n_x,n_y,title){
  outside.below.color='grey'
  na.color='black'
  outside.above.color='white'
  col=tim.colors(n = 64)
  
  #zlim=c(min(sst.anom[, plot.period]),max(sst.anom[, plot.period]))
  # z=sst.anom.land.sea[, plot.period]
  
  
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA
  
  z[which(z<zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z>zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[1] <- zlim[1] - zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
  
  #col <- c(outside.below.color, col, outside.above.color, na.color) # we construct the new color range by including: na.color and na.outside
  col <- c(outside.below.color, col, na.color)
  
  image.plot(sst.long,sst.lat,z=matrix(z,nrow=n_x,ncol=n_y,by=T),  zlim=zlim, col=col,ylab="Lattitude (deg)",xlab="Longitude (deg)",main=paste(title,"" ,month," ", year),cex.axis=.625,cex.lab=.775, cex.main=.90,xaxs="i") 
  
}


#######################################################
####### Find Month for Plotting  ######################
#######################################################
findMonth=function(curPeriod){
  numMonth=curPeriod%%12
  
  if(numMonth==0){
    curMonth="Dec."
  }
  
  if(numMonth==1){
    curMonth="Jan."
  }
  
  if(numMonth==2){
    curMonth="Feb."
  }
  
  if(numMonth==3){
    curMonth="March"
  }
  if(numMonth==4){
    curMonth="April"
  }
  if(numMonth==5){
    curMonth="May"
  }
  if(numMonth==6){
    curMonth="June"
  }
  if(numMonth==7){
    curMonth="July"
  }
  if(numMonth==8){
    curMonth="Aug."
  }
  if(numMonth==9){
    curMonth="Sept."
  }
  
  if(numMonth==10){
    curMonth="Oct."
  }
  if(numMonth==11){
    curMonth="Nov."
  }
  
  return(curMonth)
}


#######################################################
####### Create month/year label  ######################
#######################################################

yrMonthLabel=function(curPeriod,curYr){
  numMonth=curPeriod%%12
  
  curMontNum=ifelse(numMonth==0,12,numMonth)
  
  tempLastDigit=round((curYr*.01-as.integer(curYr*.01))*100)
  
  if(tempLastDigit<10){
    returnLabel=paste(curMontNum,"/0",tempLastDigit, sep = "")
  }else{
    returnLabel=paste(curMontNum,"/",tempLastDigit, sep = "")
  }
  
  
  return(returnLabel)
}


#######################################################
####### Find Year for Plotting  #######################
#######################################################

findYear=function(strYear,curPerd){
  if((curPerd%%12)==0){
    returnYr=strYear+as.integer(curPerd/12)-1
  }else{
    returnYr=strYear+as.integer(curPerd/12)
  }
  return(returnYr)
}

#######################################################
####### convertFullPlotVec  ###########################
#######################################################
#### To plot on the full domain with missing values
convertFullPlotVec=function(reducedVec,nonNaIndex,rawNumLocs){
  returnVec=rep(NA,rawNumLocs)
  returnVec[nonNaIndex]=reducedVec
  return(returnVec)
}

########################################
##### Scale Function ###################
########################################
scaleFun=function(x,orginalData){
  (x-min(orginalData))/(max(orginalData)-min(orginalData))
}



########################################
##### unscale Function #################
########################################
unscaleFun=function(x,orginalData){
  x*(max(orginalData)-min(orginalData))+min(orginalData)
}

###################################################################
#### Transform Reservior hyperparameters  #########################
###################################################################

resHypTrans=function(curTheta){
  -log(1/(curTheta-0)-1)
}

resHypTransBack=function(curGam){
  0+1/(1+exp(-curGam))
}



########################################################
#### Transform Pi hyperparameters  #####################
########################################################
piTrans=function(curTheta){
  -log(pi_b_a/(curTheta-piMin)-1)
}

piTransBack=function(curGam){
  piMin+pi_b_a/(1+exp(-curGam))
}


###################################################
#### Transform a hyperparameters  #################
###################################################
aParmTrans=function(curTheta){
  -log(aParm_b_a/(curTheta-aMin)-1)
}

aParmTransBack=function(curGam){
  aMin+aParm_b_a/(1+exp(-curGam))
}

###################################################
#### Transform U Parameters #######################
###################################################
uParmTrans=function(curTheta){
  -log(u_b_a/(curTheta-uMin)-1)
}

uParmTransBack=function(curGam){
  uMin+u_b_a/(1+exp(-curGam))
}


#######################################################
#### Transform delta Parameters #######################
#######################################################

deltaTrans=function(curTheta){
  -log(delta_b_a/(curTheta-deltaMin)-1)
}

deltaTransBack=function(curGam){
  deltaMin+delta_b_a/(1+exp(-curGam))
}

###################################################
#### Transform U Parameters #######################
###################################################
rhoTrans=function(curTheta){
  -log(rho_b_a/(curTheta-rhoMin)-1)
}

rhoTransBack=function(curGam){
  rhoMin+rho_b_a/(1+exp(-curGam))
}


#######################################################
#### Spectral Radius  #################################
#######################################################
specRadius=function(specMatrixCalc){
  abs(eigen(specMatrixCalc,only.values=TRUE)$values[1])
}

#######################################################
#### Forward DA Trans  ################################
#######################################################
forwardDATrans=function(alpha,aVal,wMat){
  alpha-log((2*aVal)/(wMat+aVal)-1)
  
}

#########################################################
#### Backwards DA Trans  ################################
#########################################################
backwardDATrans=function(alpha,aVal,wMat){
  -aVal+(2*aVal)/(1+exp(-wMat+alpha))
}

#########################################################
#### sigma MH Alpha Function ############################
#########################################################

calcSigmaMH=function(curSig,curInd,it.Count,p.star){
  
  if(curInd){
    c=curSig/(p.star*(1-p.star))
    returnVal=curSig+(c*(1-p.star))/it.Count
  }else{
    c=curSig/(p.star*(1-p.star))
    returnVal=curSig-(c*(p.star))/it.Count
  }
  
}

################################################################
#### Calculate Fourier Basis Matrix ############################
################################################################
calcFourierAlpha=function(numLocsBasis){
  fourierBasisMat=matrix(NA,numLocsBasis,numLocsBasis)
  fourierBasisMat[,1]=rep(1/sqrt(numLocsBasis),numLocsBasis)
  tempSeq=1:numLocsBasis
  for(i in 1:(numLocsBasis/2-1)){
    
    fourierBasisMat[,(i-1)*2+2]=sqrt(2/numLocsBasis)*cos((2*pi*i*tempSeq)/numLocsBasis   )
    fourierBasisMat[,(i-1)*2+3]=sqrt(2/numLocsBasis)*sin((2*pi*i*tempSeq)/numLocsBasis   )
    
  }
  
  fourierBasisMat[,numLocsBasis]=sqrt(1/numLocsBasis)*cos(pi*tempSeq)
  return(fourierBasisMat)
}

#############################################################
#### Find Fourier Modes with most Energy ####################
############################################################
findTopEnergy=function(fullAlphaMat,numLocsBasis){
  
  returnEng=rep(NA,length(1:(numLocsBasis/2-1)))
  for(i in 1:(numLocsBasis/2-1)){
    returnEng[i]=sum((fullAlphaMat[(i-1)*2+2,]-mean(fullAlphaMat[(i-1)*2+2,]))^2)+sum((fullAlphaMat[(i-1)*2+3,]-mean(fullAlphaMat[(i-1)*2+3,]))^2)
  }
  
  returnEng
}


###################################################################
#### Horse Shoe Prior Functions ###################################
###################################################################

#############################################
##############  Draw Beta ###################
#############################################

draw_beta <- function(mu_n, Lambda_n_inv, sigma){
  rmvn(1,mu_n,sigma * Lambda_n_inv,isChol=F)
}
#############################################
##############  Draw lambda #################
#############################################

draw_lambda <- function(lambda, beta, sigma, tau){
  gamma_l <- 1 / lambda^2
  u1 <- runif(length(lambda), 0, 1 / (1 + gamma_l))
  trunc_limit <- (1 - u1) / u1
  mu2_j <- (beta / (sqrt(sigma) * tau))^2
  rate_lambda <- (mu2_j / 2)
  ub_lambda <- pexp(trunc_limit, rate_lambda)
  u2 <- runif(length(ub_lambda), 0, ub_lambda)
  gamma_l <- qexp(u2, rate_lambda)
  return(1 / sqrt(gamma_l))
}


#############################################
##############  Draw tau ####################
#############################################

draw_tau <- function(lambda, beta, sigma, tau){
  shape_tau <- 0.5 * (length(lambda) + 1)
  gamma_t <- 1 / tau^2
  u1 <- runif(1, 0, 1 / (1 + gamma_t))
  trunc_limit_tau <- (1 - u1) / u1
  mu2_tau <- sum((beta / (sqrt(sigma) * lambda))^2)
  rate_tau <- (mu2_tau / 2)
  ub_tau <- pgamma(trunc_limit_tau, shape=shape_tau, rate=rate_tau)
  u2 <- runif(1, 0, ub_tau)
  gamma_t <- qgamma(p = u2, shape = shape_tau, rate = rate_tau)
  return(1 / sqrt(gamma_t))
}

###########################################################
######### Caculates Coords on a Circle ####################
###########################################################
#### given a center point: (centerX,centerY)
#### radius and angle this function calculates
#### the point on a circle
circleCordFunct=function(centerX,centerY,radius,angle){
  
  xCoord=centerX+radius*cos(angle)
  yCoord=centerY+radius*sin(angle)
  return(c(xCoord,yCoord))
}

###########################################################################
################# C++ function bisquare Function ########################## 
########################################################################### 
cppFunction('NumericVector biSquare(NumericVector gridX, NumericVector gridY, double xVal, double yVal, double radius ) {
            int n = gridX.size();
            NumericVector out(n);
            for(int i = 0; i < n; ++i) {
            float dist=sqrt(pow(gridX[i]-xVal,2)+ pow(gridY[i]-yVal,2));
            if(dist<=radius){
            out[i]=pow(1-pow(dist/radius,2),2);
            }else{
            out[i]=0;
            }
            }
            return out;
            }')


###############################################################################
############### Functions for Laplacian EigenMaps #############################
###############################################################################

mds.edm1 <- function(X) {
  return(as.matrix(dist(X)))
}

graph.knn <- function(D,k) {
  #
  #  Returns an nxn 0-1 matrix KNN.
  #    KNN[i,j]=1 iff j is a neighbor of i.
  #    Note that KNN may be asymmetric.
  #  We assume that D[i,j]>0 except when i=j.
  #
  n <- nrow(D)
  KNN <- matrix(0,nrow=n,ncol=n)
  near <- 2:(k+1)
  for (i in 1:n) {
    v <- D[i,]
    j <- order(v)
    j <- j[near]
    KNN[i,j] <- 1
  }
  return(KNN)
}


##########################################################################
######### either graph.adj or graph.heat can be used #####################
########  to construct adjaceny matrix ###################################
##########################################################################

graph.adj <- function(KNN) {
  #
  #  Uses the output of graph.knn to construct an adjacency matrix.
  #  Vertices i & j are connected iff either j is a neighbor of i
  #    or i is a neighbor of j.
  #
  return(pmax(KNN,t(KNN)))
}


graph.heat <- function(D,sigma) {
  #
  #  Assigns edge weights equal to similarities based on the heat kernel.
  #
  return(exp(-sigma*D^2))
}

graph.laplacian <- function(Gamma) {
  #
  #  Computes the Laplacian matrix of a weighted graph.
  #
  tot <- apply(Gamma,1,sum)
  return(diag(tot)-Gamma)
}


##############################################
##### Lorenz 96 Functions ####################
##############################################

mod = function(i,n) return((i-1)%%n+1)



##########Standard deterministic model
lorenz96 = function(x0, theta, dt, n, M){
  ####################################################
  # Take M steps of the L96 model with step size dt 
  # Uses 1st-order Euler scheme
  #---------------------------------------------------
  # x0: (nx1) initial conditions
  # theta: scalar forcing parameter (often called F)
  # dt: internal step size
  # n: number of spatial locations (n=40 is standard)
  # M: number of internal steps:   M=delta/dt
  ####################################################
  ###### This function is just for one time period 
  ###### mod operator makes 
  
  xx = x0
  dx = rep(0,n)
  
  for (j in 1:M){
    for (i in 1:n){
      dx[i] = ( (xx[mod(i+1,n)] - xx[mod(i-2,n)]) * xx[mod(i-1,n)] - xx[i] + theta) * dt
    }
    xx = xx + dx
  }
  return(xx)
}

##########simulate with error

lorenz40Error = function(x0, theta, dt, n, M,lorenz40Scale,delta_t){
  ####################################################
  # Take M steps of the L96 model with step size dt 
  # Uses 1st-order Euler scheme
  #---------------------------------------------------
  # x0: (nx1) initial conditions
  # theta: scalar forcing parameter (often called F)
  # dt: internal step size
  # n: number of spatial locations (n=40 is standard)
  # M: number of internal steps:   M=delta/dt
  ####################################################
  xx = x0
  dx = rep(0,n)
  
  for (j in 1:M){
    for (i in 1:n){
      dx[i] = ( (xx[mod(i+1,n)] - xx[mod(i-2,n)]) * xx[mod(i-1,n)] - xx[i] + theta+(lorenz40Scale/delta_t)*rnorm(1,0,sqrt(delta_t))) * dt
    }
    xx = xx + dx
  }
  return(xx)
}


###################################################
###################################################
##### Simulates Lotka-Volterra Model ##############
###################################################
###################################################

LotVmod <- function (Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    dx = x*(alpha - beta*y)
    dy = -y*(gamma - delta*x)
    return(list(c(dx, dy)))
  })
}




####################################################################
##### LV INTERACTION Variable Euler Function #######################
####################################################################

##########Standard deterministic model
lvBetaEuler = function(x0, dt, curt, M,lvDelta,lvRho,lvH,lvOmega0,lvGamma,lvSigmaBeta){
  ####################################################
  # Take M steps of the L96 model with step size dt 
  # Uses 1st-order Euler scheme
  
  
  uPrime99=(12*lvH*(.99-(1+lvRho))^2 )/lvDelta^4 +(4*lvH)/lvDelta^2
  uPrime94=(12*lvH*(.94-(1+lvRho))^2 )/lvDelta^4 +(4*lvH)/lvDelta^2
  
  xx = x0
  dx = 0
  
  for (j in 1:M){
    
    (2*pi*exp((2*lvH)/rnorm(1,0,sd=lvSigmaBeta)  ))/sqrt(abs(uPrime99)*uPrime94 )
    
    # dx= (  -((4*lvH*(xx-(1+lvRho))^3 )/lvDelta^4 -(4*lvH*(xx-(1+lvRho)) )/lvDelta^2) +lvGamma*cos(lvOmega0*curt)  ) * dt
    dx= (  -((4*lvH*(xx-(1+lvRho))^3 )/lvDelta^4 -(4*lvH*(xx-(1+lvRho)) )/lvDelta^2) +lvGamma*cos(lvOmega0*curt)+rnorm(1,0,sd=lvSigmaBeta)  ) * dt
    # dx= (  -((4*lvH*(xx-(1+lvRho))^3 )/lvDelta^4 -(4*lvH*(xx-(1+lvRho)) )/lvDelta^2) +lvGamma*cos(lvOmega0*curt) +(2*pi*exp((2*lvH)/rnorm(1,0,sd=lvSigmaBeta)  ))/sqrt(abs(uPrime99)*uPrime94  )) * dt
    xx = xx + dx
  }
  return(xx)
}


#####################################################
##### ST LV  Euler Function ########################
####################################################

stLVEuler = function(x0,y0, dt, n, M,curBeta,lvMu,lvSigmaMain,laggedX0Vec,laggedY0Vec,neighborsMat,lvD){
  ####################################################
  # Uses 1st-order Euler scheme
  #---------------------------------------------------
  # x0: (nx1) initial conditions
  # dt: internal step size
  # n: number of spatial locations 
  # M: number of internal steps:   M=delta/dt
  ####################################################
  
  
  ##### X is species 1
  xx = x0
  dx = rep(0,n)
  
  ##### y is species 2
  yy = y0
  dy = rep(0,n)
  
  
  spec1ReturnMMat=matrix(NA,M,n)
  spec2ReturnMMat=matrix(NA,M,n)
  
  for (j in 1:M){
    for (i in 1:n){
      
      # dx[i] = (lvMu*xx[i]*(1-xx[i]-curBeta*yy[i])+sqrt(lvSigmaMain)*rnorm(1,0,1)*xx[i] ) * dt
      # dy[i] = (lvMu*yy[i]*(1-yy[i]-curBeta*xx[i])+sqrt(lvSigmaMain)*rnorm(1,0,1)*yy[i] ) * dt
      
      # dx[i] = (lvMu*xx[i]*(1-xx[i]-curBeta*yy[i])+sqrt(lvSigmaMain)*rnorm(1,0,1)*xx[i]+lvD*sum(xx[neighborsMat[i,]]-xx[i]) ) * dt
      # dy[i] = (lvMu*yy[i]*(1-yy[i]-curBeta*xx[i])+sqrt(lvSigmaMain)*rnorm(1,0,1)*yy[i]+lvD*sum(yy[neighborsMat[i,]]-yy[i]) ) * dt
      
      # dx[i] = (lvMu*xx[i]*(1-laggedX0Vec[i]-curBeta*yy[i])+sqrt(lvSigmaMain)*rnorm(1,0,1)*xx[i]+lvD*sum(xx[neighborsMat[i,]]-xx[i]) ) * dt
      # dy[i] = (lvMu*yy[i]*(1-laggedY0Vec[i]-curBeta*xx[i])+sqrt(lvSigmaMain)*rnorm(1,0,1)*yy[i]+lvD*sum(yy[neighborsMat[i,]]-yy[i]) ) * dt
      
      # dx[i] = (lvMu*xx[i]*(1-laggedX0Vec[i]-curBeta*yy[i])+sqrt(lvSigmaMain)*rnorm(1,0,1)*xx[i]+lvD*sum(xx[neighborsMat[i,]]-xx[i]) ) * dt
      # dy[i] = (lvMu*yy[i]*(1-laggedY0Vec[i]-curBeta*xx[i])+sqrt(lvSigmaMain)*rnorm(1,0,1)*yy[i]+lvD*sum(yy[neighborsMat[i,]]-yy[i]) ) * dt
      
      dx[i] = (lvMu*xx[i]*(1-laggedX0Vec[i]-curBeta*yy[i])+rtruncnorm(1,0,Inf,0,sqrt(lvSigmaMain))+sqrt(lvSigmaMain)*rnorm(1,0,1)*xx[i]+lvD*sum(xx[neighborsMat[i,]]-xx[i]) ) * dt
      dy[i] = (lvMu*yy[i]*(1-laggedY0Vec[i]-curBeta*xx[i])+rtruncnorm(1,0,Inf,0,sqrt(lvSigmaMain))+sqrt(lvSigmaMain)*rnorm(1,0,1)*yy[i]+lvD*sum(yy[neighborsMat[i,]]-yy[i]) ) * dt
      
      
      # cat("i ", i , " j ", j , " xx ", xx[i], " yy ", yy[i],"\n")
    }
    xx = xx + dx
    yy = yy + dy
    spec1ReturnMMat[j,]=xx
    spec2ReturnMMat[j,]=yy
  }
  return(list(xx,yy,spec1ReturnMMat,spec2ReturnMMat))
}




##############################################
##### gQESNCalcFor    ########################
##############################################
####allows us to do parallel processing

gQESNCalcFor=function(inSampHMat,curInSampY,typeResp,curRidge,outSampHMat){
  glmFitObj=glmnet(inSampHMat,curInSampY ,family=typeResp,alpha=0,nlambda = 20,lambda.min.ratio = .00001)
  outSampFor=predict(glmFitObj,newx=outSampHMat,type='response',s=curRidge)
  return(outSampFor)
}

##############################################
## Funciton For Parallel Processing ##########
##############################################
iblkcol <- function(a, chunks) {
  n <- ncol(a)
  i <- 1
  nextElem <- function() {
    if (chunks <= 0 || n <= 0) stop('StopIteration')
    m <- ceiling(n / chunks)
    r <- seq(i, length=m)
    i <<- i + m
    n <<- n - m
    chunks <<- chunks - 1
    a[,r, drop=FALSE]
  }
  structure(list(nextElem=nextElem), class=c('iblkcol', 'iter'))
}
nextElem.iblkcol <- function(obj) obj$nextElem()

###############################################
##### Calculate EOFS for N X T Data set #######
###############################################
calcEOFS=function(rawDataMatEOF,numBasisFuncts,inSampSeq){
  
  tempDataCalc=rawDataMatEOF[,inSampSeq]
  tempCovMat=((tempDataCalc-rowMeans(tempDataCalc))%*%t((tempDataCalc-rowMeans(tempDataCalc))))/length(inSampSeq)
  tempEigObj=eigen(tempCovMat)
  tempEignMat=tempEigObj$vectors
  
  phiMat=t(tempEignMat[,1:numBasisFuncts])
  returnObj=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%rawDataMatEOF)
  
  leftOverPhiMat=t(tempEignMat[,(numBasisFuncts+1):nrow(rawDataMatEOF)])
  
  return(list(phiMat=phiMat,basisCoefs=returnObj,allEigVals=tempEigObj$values,leftOverPhiMat=leftOverPhiMat))
}

###########################################################################
################# C++ function bisquare Function ########################## 
########################################################################### 

cppFunction('NumericVector biSquare(NumericVector gridX, NumericVector gridY, double xVal, double yVal, double radius ) {
            int n = gridX.size();
            NumericVector out(n);
            for(int i = 0; i < n; ++i) {
            float dist=sqrt(pow(gridX[i]-xVal,2)+ pow(gridY[i]-yVal,2));
            if(dist<=radius){
            out[i]=pow(1-pow(dist/radius,2),2);
            }else{
            out[i]=0;
            }
            }
            return out;
            }')

############################################
##### fourierHiddenFunction#################
############################################

#######Fourier
# eofHiddenFunction=function(curNh,tempDataCalc,numBasisFuncts){
# 
#   fullPhiMat=calcFourierAlpha(curNh)
# 
#   # print(dim(t(fullPhiMat)))
#   # print(dim(tempDataCalc))
# 
#   fullAlphaMat=t(fullPhiMat)%*%tempDataCalc
# 
# 
#   allBasisEng=findTopEnergy(fullAlphaMat,curNh)
# 
#   tempTopEngIndex=(sort(allBasisEng,index.return=T,decreasing=T)$ix)[1:(numBasisFuncts/2)]
# 
#   phiMat=fullPhiMat[,sort(c(tempTopEngIndex*2+1,tempTopEngIndex*2+2))]
#   phiMat=t(phiMat)
# 
#   returnObj=t(phiMat%*%tempDataCalc)
#   scaledObj=scale(returnObj)
# 
#   scaledMeansH2Mat=attr(scaledObj,"scaled:center")
#   scaledSDsH2Mat=attr(scaledObj,"scaled:scale")
# 
# 
#   return(list(scaledInsampH2Mat=scaledObj,scaledMeansH2Mat=scaledMeansH2Mat,scaledSDsH2Mat=scaledSDsH2Mat,phiMat=phiMat))
# }
# 
# 
# #######################################################
# ####### createDimRedESNOutSamp ########################
# #######################################################
# createOutSampH2Mat=function(phiMat, scaledMeans,scaledSDs,inbetweenHMat,outSampH2Mat){
#   # unscaleOutSampObj=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%outSampH2Mat)
#   unscaleOutSampObj=t(phiMat%*%outSampH2Mat)
#   scaleOutSamp=sapply(1:ncol(unscaleOutSampObj),function(funVar) (unscaleOutSampObj[,funVar]-scaledMeans[funVar] )/scaledSDs[funVar] )
# 
#   unscaleInbetween=t(phiMat%*%inbetweenHMat)
#   scaledBetween=sapply(1:ncol(unscaleInbetween),function(funVar) (unscaleInbetween[,funVar]-scaledMeans[funVar] )/scaledSDs[funVar] )
# 
#   return(list(scaleOutSamp=scaleOutSamp,scaledBetween=scaledBetween))
# }










#####EOFs
eofHiddenFunction=function(curNh,tempDataCalc,numBasisFuncts){
  
  tempCovMat=((tempDataCalc-rowMeans(tempDataCalc))%*%t((tempDataCalc-rowMeans(tempDataCalc))))/nrow(tempDataCalc)
  tempEigObj=eigen(tempCovMat)
  tempEignMat=tempEigObj$vectors
  
  phiMat=t(tempEignMat[,1:numBasisFuncts])
  returnObj=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%tempDataCalc)
  
  
  scaledMeansH2Mat=mean(returnObj)
  scaledSDsH2Mat=sd(returnObj)
  
  scaledObj=(returnObj-scaledMeansH2Mat)/scaledSDsH2Mat
  
  
  return(list(scaledInsampH2Mat=scaledObj,scaledMeansH2Mat=scaledMeansH2Mat,scaledSDsH2Mat=scaledSDsH2Mat,phiMat=phiMat))
}



createOutSampH2Mat=function(phiMat, scaledMeans,scaledSDs,inbetweenHMat,outSampH2Mat){
  # unScaledReturnObj=t(solve(phiMat%*%t(phiMat))%*%phiMat%*%outSampH2Mat)
  unScaledReturnObj=t(phiMat%*%outSampH2Mat)
  
  scaleOutSamp=sapply(1:ncol(unScaledReturnObj),function(funVar) (unScaledReturnObj[,funVar]-scaledMeans )/scaledSDs )
  
  unscaleInbetween=t(phiMat%*%inbetweenHMat)
  scaledBetween=sapply(1:ncol(unscaleInbetween),function(funVar) (unscaleInbetween[,funVar]-scaledMeans )/scaledSDs )
  
  return(list(scaleOutSamp=scaleOutSamp,scaledBetween=scaledBetween))
}






#######################################################
##### Rotation Function ###############################
#### For two vectors, y is the target vector ##########
#######################################################
rotation = function(x,y){
  u=x/sqrt(sum(x^2))
  
  v=y-sum(u*y)*u
  v=v/sqrt(sum(v^2))
  
  cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
  
  sint=sqrt(1-cost^2);
  
  diag(length(x)) - u %*% t(u) - v %*% t(v) + 
    cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
}

#############################################
############## vectorRotation ###############
#############################################

vectorRotation=function(x,y){
  rotMat=rotation(x,y)
  x%*% rotMat
}


#############################################
##############  rFunctionBoundedTrans #######
#############################################
rFunctionBoundedTrans=function( curVal, len, minVal){
  -log(len/(curVal-minVal)-1)
}


#################################################
##############  rFunctionBoundedTransBack #######
#################################################
rFunctionBoundedTransBack=function( curVal, len, minVal){
  minVal+len/(1+exp(-curVal))
}


#################################################
######### Create Soil Moisture Embed ############
#################################################

createSMEmbedRNNData=function(curTrainLen,m,tauEmb,yTrain,rawDataInput,numRawInputVars,curXTestIndex,curTestLen,curInBetweenIndexes,curLenInbetween,numLocs){
  lenInSampleEmb=curTrainLen-(m*tauEmb)
  
  altInSampleXRaw=array(NA,c(lenInSampleEmb,m+1,numRawInputVars))
  
  for(i in 1:lenInSampleEmb){
    #####using laplacian eigenmaps
    altInSampleXRaw[i,,]=rawDataInput[seq(i,(m*tauEmb+i),by=tauEmb),]
  }
  
  #######Scale in-sample x and y
  altYInSampRaw=yTrain[(m*tauEmb+1):curTrainLen,]
  altYMean=mean(altYInSampRaw)
  altScaleFactor=sd(altYInSampRaw)
  
  altYInSamp=(altYInSampRaw-altYMean)/altScaleFactor
  
  
  meanXTrainMatrix=mean(rawDataInput[1:curTrainLen,])
  sdXTrainMatrix=sd(rawDataInput[1:curTrainLen,])
  
  
  altInSampleX=(altInSampleXRaw-meanXTrainMatrix)/sdXTrainMatrix
  
  
  designMatrix=matrix(1,lenInSampleEmb,(m+1)*numRawInputVars+1)
  for(i in 1:lenInSampleEmb){
    designMatrix[i,2:((m+1)*numRawInputVars+1)]=as.vector(altInSampleX[i,,])
  }
  
  # print(2:((m+1)*numRawInputVars+1))
  
  ###############################
  ####### Inbetween #############
  ###############################
  altInbetweenXRaw=array(NA,c(curLenInbetween,m+1,numRawInputVars))
  for(i in 1:curLenInbetween){
    altInbetweenXRaw[i,,]=rawDataInput[seq(curInBetweenIndexes[i]-(m*tauEmb),curInBetweenIndexes[i],by=tauEmb),]
  }
  
  
  #######Scale in-sample x and y
  altInbetweenX=(altInbetweenXRaw-meanXTrainMatrix)/sdXTrainMatrix
  
  inbetweenDesignMatrix=matrix(1,curLenInbetween,(m+1)*numRawInputVars+1)
  for(i in 1:curLenInbetween){
    inbetweenDesignMatrix[i,2:((m+1)*numRawInputVars+1)]=as.vector(altInbetweenX[i,,])
  }
  
  
  ####### Out Sample
  altOutSampleXRaw=array(NA,c(curTestLen,m+1,numRawInputVars))
  for(i in 1:curTestLen){
    altOutSampleXRaw[i,,]=rawDataInput[seq(curXTestIndex[i]-(m*tauEmb),curXTestIndex[i],by=tauEmb),]
  }
  
  
  #######Scale in-sample x and y
  altOutSampleX=(altOutSampleXRaw-meanXTrainMatrix)/sdXTrainMatrix
  
  designMatrixOutSample=matrix(1,curTestLen,(m+1)*numRawInputVars+1)
  for(i in 1:curTestLen){
    designMatrixOutSample[i,2:((m+1)*numRawInputVars+1)]=as.vector(altOutSampleX[i,,])
  }
  
  
  
  ######additive scale matrix
  addScaleMat=matrix(altYMean,numLocs,curTestLen)
  return(list(curInSampMean=altYMean,curInSampSD=altScaleFactor,curYInSamp=altYInSamp,curInSampX=altInSampleX,curOutSampX=altOutSampleX,lenInSampleEmb=lenInSampleEmb,designMatrix=designMatrix,designMatrixOutSample=designMatrixOutSample,curTestLen=curTestLen,addScaleMat=addScaleMat,inbetweenDesignMatrix=inbetweenDesignMatrix ))
}

#################################################
##############  createPrecipESNData #############
#################################################
#### Create data for ESN for Precip. applicaiton

createPrecipESNData=function(curM,tauEmb,curTrainLen,curXTestIndex,curTestLen,rawDataInput,yTrain,numRawInputVars,inbetweenIndexes,inbetweenLen,nmfInSampMonthMeans,curYTrainSeq){
  
  lenInSampleEmb=curTrainLen-(curM*tauEmb)
  
  
  altInSampleXRaw=array(NA,c(lenInSampleEmb,curM+1,numRawInputVars))
  
  for(i in 1:lenInSampleEmb){
    curSeq=seq(i,(curM*tauEmb+i),by=tauEmb)
    # cat( "i ", i , " curseq ", curSeq, "\n")
    altInSampleXRaw[i,,]=rawDataInput[curSeq,]
  }
  
  #######Scale in-sample x and y
  insampYSeq=curYTrainSeq[which((curM*tauEmb+tau)<curYTrainSeq)]
  
  
  TempaltYInSamp=matrix(NA,length(insampYSeq),numLocs)
  for(i in 1:length(insampYSeq)){
    
    curIndex=insampYSeq[i]
    # cat( "i ", i , " curIndex ", curIndex, "\n")
    TempaltYInSamp[i,]=rawDataOutput[,curIndex]
  }
  
  
  altScaleFactor=sd(TempaltYInSamp)
  altScaleMean=mean(TempaltYInSamp)
  
  altYInSamp=(TempaltYInSamp-altScaleMean)/altScaleFactor
  
  ######Only scale the EOFs
  altInSampleX=array(NA,c(lenInSampleEmb,curM+1,numRawInputVars))
  meanXTrainMatrix=matrix(NA,curM+1,numRawInputVars)
  sdXTrainMatrix=matrix(NA,curM+1,numRawInputVars)
  for(i in 1:(curM+1)){
    for(j in 1:numRawInputVars){
      
      meanXTrainMatrix[i,j]=mean(altInSampleXRaw[,i,j])
      sdXTrainMatrix[i,j]=sd(altInSampleXRaw[,i,j])
      altInSampleX[,i,j]=(altInSampleXRaw[,i,j]-meanXTrainMatrix[i,j])/sdXTrainMatrix[i,j]
    }
    
  }
  
  
  
  
  designMatrix=matrix(1,lenInSampleEmb,(curM+1)*numRawInputVars+1)
  for(i in 1:lenInSampleEmb){
    designMatrix[i,2:((curM+1)*numRawInputVars+1)]=as.vector(altInSampleX[i,,])
  }
  
  
  ##########inbetween
  altInbetweenXRaw=array(NA,c(inbetweenLen,curM+1,numRawInputVars))
  for(i in 1:inbetweenLen){
    curSeq=seq(inbetweenIndexes[i]-(curM*tauEmb),inbetweenIndexes[i],by=tauEmb)
    # cat( "i ", i , " curseq ", curSeq, "\n")
    altInbetweenXRaw[i,,]=rawDataInput[curSeq,]
  }
  
  
  
  altInbetweenX=array(NA,c(inbetweenLen,curM+1,numRawInputVars))
  for(i in 1:(curM+1)){
    for(j in 1:numRawInputVars){
      altInbetweenX[,i,j]=(altInbetweenXRaw[,i,j]-meanXTrainMatrix[i,j])/sdXTrainMatrix[i,j]
    }
    
  }
  
  designMatrixInbetween=matrix(1,inbetweenLen,(curM+1)*numRawInputVars+1)
  for(i in 1:inbetweenLen){
    designMatrixInbetween[i,2:((curM+1)*numRawInputVars+1)]=as.vector(altInbetweenX[i,,])
  }
  
  
  ####### Out Sample
  altOutSampleXRaw=array(NA,c(curTestLen,curM+1,numRawInputVars))
  for(i in 1:curTestLen){
    curSeq=seq(curXTestIndex[i]-(curM*tauEmb),curXTestIndex[i],by=tauEmb)
    # cat( "i ", i , " curseq ", curSeq, "\n")
    altOutSampleXRaw[i,,]=rawDataInput[curSeq,]
  }
  
  
  altOutSampleX=array(NA,c(curTestLen,curM+1,numRawInputVars))
  for(i in 1:(curM+1)){
    for(j in 1:numRawInputVars){
      altOutSampleX[,i,j]=(altOutSampleXRaw[,i,j]-meanXTrainMatrix[i,j])/sdXTrainMatrix[i,j]
    }
    
  }
  
  
  
  designMatrixOutSample=matrix(1,curTestLen,(curM+1)*numRawInputVars+1)
  for(i in 1:curTestLen){
    designMatrixOutSample[i,2:((curM+1)*numRawInputVars+1)]=as.vector(altOutSampleX[i,,])
  }
  
  ######additive scale matrix
  addScaleMat=matrix(altScaleMean,numLocs,curTestLen)
  
  return(list(lenInSampleEmb=lenInSampleEmb,addScaleMat=addScaleMat,altScaleFactor=altScaleFactor,designMatrix=designMatrix,designMatrixOutSample=designMatrixOutSample,altYInSamp=altYInSamp,designMatrixInbetween=designMatrixInbetween,insampYSeq=insampYSeq))
}


#################################################
##############  sstESNFor #######################
#################################################
#####need a special function for the ESN to get forecasts
#### since the indexing is weird

sstESNFor=function(curNh,lenInSampleEmb,designMatrix,sampVecESN,wWidth,piWESN,piUESN,uWidth,nColsU,curDelta,stratValuesXtemp,quadInd,hMatdim,altYInSamp,ridgeMat,designMatrixOutSample,curLenAllOut,curIndexOutSampPerds,altScaleFactor,addScaleMat,Phi.tr){
  
  reservoirObj=genResCPP(curNh,lenInSampleEmb,designMatrix,sampVecESN,wWidth,piWESN,piUESN,uWidth,nColsU,curDelta,stratValuesXtemp,quadInd,hMatdim)
  
  xTemp=reservoirObj$xTempLast
  hMat=reservoirObj$hMat
  wMatScaled=reservoirObj$wMatESNFun
  uTempESN=reservoirObj$uMatESNFun
  
  vMatESN=t(altYInSamp)%*%t(hMat)%*%solve(hMat%*%t(hMat)+ridgeMat )
  
  uProdMatOutSamp=uTempESN%*%t(designMatrixOutSample)
  
  outSampHMatObj= createHMat( curNh,  curLenAllOut, wMatScaled,  uProdMatOutSamp,xTemp,quadInd,hMatdim)
  
  hMatOutSamp=outSampHMatObj$hMat
  redHMatOutSamp=hMatOutSamp[,curIndexOutSampPerds]
  
  scaledOutSampleFor=altScaleFactor*(vMatESN%*%redHMatOutSamp)+addScaleMat[1,1]
  
  
  tempAllLocsFor=t(scaledOutSampleFor)%*%Phi.tr
  return(tempAllLocsFor)
}


#############################################
################ image.nan ##################
#############################################

image.nan<- function(z,  zlim, col, na.color='white', outside.below.color='gray', outside.above.color='black',...)
{
  zstep <- (zlim[2] - zlim[1]) / length(col); # step in the color palette
  newz.below.outside <- zlim[1] - 2 * zstep # new z for values below zlim
  newz.above.outside <- zlim[2] + zstep # new z for values above zlim
  newz.na <- zlim[2] + 2 * zstep # new z for NA
  
  z[which(z<zlim[1])] <- newz.below.outside # we affect newz.below.outside
  z[which(z>zlim[2])] <- newz.above.outside # we affect newz.above.outside
  z[which(is.na(z>zlim[2]))] <- newz.na # same for newz.na
  
  zlim[1] <- zlim[1] - 2 * zstep # extend lower limit to include below value
  zlim[2] <- zlim[2] + 2 * zstep # extend top limit to include the two new values above and na
  
  col <- c(outside.below.color, col[1], col, outside.above.color, na.color) #correct by including col[1] at bottom of range
  
  return(list(z=z,  zlim=zlim, col=col)) # we finally call image(...)
}


#############################################
### Soil Moisture Genetic Algorithm #########
#############################################
smGAFunction=function(curM,curRidge,layerOneDimRed,curNhTwo,deltaVec,curDeltaTwo){
  
  ######### Forces these values to be discrete
  curM=as.integer(curM)
  layerOneDimRed=as.integer(layerOneDimRed)
  if(!(layerOneDimRed%%2)==0 & dimRedESN=="Fourier"){
    layerOneDimRed=layerOneDimRed+1
    
  }
  
  curNhTwo=as.integer(curNhTwo)
  
  cat("curM ",curM, "\n")
  cat("curRidge ",curRidge, "\n")
  cat("layerOneDimRed ",layerOneDimRed, "\n")
  cat("curNhTwo ",curNhTwo, "\n")
  cat("curDeltaOne ",deltaVec[1], "\n")
  cat("curDeltaTwo ",deltaVec[2], "\n")
  
  curTrainLen=validLen
  curXTestIndex=xValTestIndex
  curYTestLen=validTestLen
  curInbetweenIndexes=validInbetweenIndexes
  curInbetweenLen=validInbetweenLen
  curTestLen=validTestLen
  curYTestAllLocs=yValidALLLocs
  
  
  createESNDataObj=createSMEmbedRNNData(curTrainLen,curM,tauEmb,yTrain,rawDataInput,numRawInputVars,curXTestIndex,curTestLen,curInbetweenIndexes,curInbetweenLen,numLocs )
  
  lenInSampleEmb=createESNDataObj$lenInSampleEmb
  addScaleMat=createESNDataObj$addScaleMat
  altScaleFactor=createESNDataObj$curInSampSD
  designMatrix=createESNDataObj$designMatrix
  designMatrixOutSample=createESNDataObj$designMatrixOutSample
  altYInSamp=createESNDataObj$curYInSamp
  designMatrixInbetween=createESNDataObj$inbetweenDesignMatrix
  
  
  
  ##########first layer
  setParObjFirst=setParsEESN(curRidge ,curNhOne,numLocs,curM)
  sampVecESNOne=setParObjFirst$sampVecESN
  stratValuesXtempOne=setParObjFirst$stratValuesXtemp
  
  
  
  setParObjTwo=setParsEESN(curRidge ,curNhTwo,layerOneDimRed,curMTwo)
  sampVecESNTwo=setParObjTwo$sampVecESN
  stratValuesXtempTwo=setParObjTwo$stratValuesXtemp
  
  if(featureLinkInd){
    
    if(quadInd){
      ridgeMat=diag(curRidge,2*curNhTwo+(numLayers-1)*layerOneDimRed)
    }else{
      ridgeMat=diag(curRidge,curNhTwo+(numLayers-1)*layerOneDimRed)
    }
    
  }else{
    ridgeMat=setParObjTwo$ridgeMat
  }
  
  
  nColsUOne=ncol(designMatrix)
  nColsUTwo=layerOneDimRed+1
  
  
  forMatMultiESN=array(NA,c(ensembleLen,curYTestLen,numLocs))
  for(iESN in 1:ensembleLen){
    forMatMultiESN[iESN,,]=multiLayerDeepESNFor(curNhOne,lenInSampleEmb,designMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUOne,deltaVec,stratValuesXtempOne,curTestLen,dimRedESN,layerOneDimRed,curNhTwo,sampVecESNTwo,uWidthTwo,piWESNTwo,piUESNTwo,wWidthTwo,nColsUTwo,stratValuesXtempTwo,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,designMatrixOutSample,designMatrixInbetween,curInbetweenLen,featureLinkInd,numLayers)$outSampFor
  }
  
  finalEsembelFor=sapply(1:numLocs,function(funVar) colMeans(forMatMultiESN[,,funVar]))
  
  allEnsembleFor=t(phiSM%*%t(finalEsembelFor))
  
  return(mean((allEnsembleFor[validMayIndex,]-curYTestAllLocs[validMayIndex,])^2))
  
}


####################################################
## Soil Moisture Classic ESN Genetic Algorithm #####
####################################################
smEsnGAFunction=function(curM,curRidge,curNh,curDelta){
  
  ######### Forces these values to be discrete
  curM=as.integer(curM)
  curNh=as.integer(curNh)
  
  cat("curM ",curM, "\n")
  cat("curRidge ",curRidge, "\n")
  cat("curNh ",curNh, "\n")
  cat("curDelta ",curDelta, "\n")
  
  ########Variables for validation
  curTrainLen=validLen
  curXTestIndex=xValTestIndex
  curYTestLen=validTestLen
  curInbetweenIndexes=validInbetweenIndexes
  curInbetweenLen=validInbetweenLen
  curTestLen=validTestLen
  curYTestAllLocs=yValidALLLocs
  
  
  createESNDataObj=createSMEmbedRNNData(curTrainLen,curM,tauEmb,yTrain,rawDataInput,numRawInputVars,curXTestIndex,curTestLen,curInbetweenIndexes,curInbetweenLen,numLocs )
  
  lenInSampleEmb=createESNDataObj$lenInSampleEmb
  addScaleMat=createESNDataObj$addScaleMat
  altScaleFactor=createESNDataObj$curInSampSD
  designMatrix=createESNDataObj$designMatrix
  designMatrixOutSample=createESNDataObj$designMatrixOutSample
  altYInSamp=createESNDataObj$curYInSamp
  designMatrixInbetween=createESNDataObj$inbetweenDesignMatrix
  
  
  setParObj=setParsEESN(curRidge ,curNh,numLocs,curM)
  
  
  
  sampVecESN=setParObj$sampVecESN
  stratValuesXtemp=setParObj$stratValuesXtemp
  if(quadInd){
    ridgeMat=diag(curRidge,2*curNh)
  }else{
    ridgeMat=diag(curRidge,curNh)
  }
  nColsU=ncol(designMatrix)
  
  forMatMultiESN=array(NA,c(ensembleLen,curYTestLen,numLocs))
  
  for(iEnsem in 1:ensembleLen){
    # scaledOutSampleFor=calcESNForecastsSpinForward(curNh,lenInSampleEmb,curTestLen,designMatrix,designMatrixOutSample,sampVecESN,wWidth,piWESN,piUESN,uWidth,nColsU,curDelta,stratValuesXtemp,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,quadInd,designMatrixInbetween,curInbetweenLen)
    scaledOutSampleFor=calcESNForecastsSpinForward(curNh,lenInSampleEmb,curTestLen,designMatrix,designMatrixOutSample,sampVecESN,wWidth,piWESN,piUESN,uWidth,nColsU,curDelta,stratValuesXtemp,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,quadInd,designMatrixInbetween,curInbetweenLen)
    forMatMultiESN[iEnsem,,]=t(scaledOutSampleFor)
  }
  
  
  finalEsembelFor=sapply(1:numLocs,function(funVar) colMeans(forMatMultiESN[,,funVar]))
  
  allEnsembleFor=t(phiSM%*%t(finalEsembelFor))
  
  return(mean((allEnsembleFor[validMayIndex,]-curYTestAllLocs[validMayIndex,])^2))
}



#############################################
### Deep Generic Genetic Algorithm ##########
#############################################

deepGAFunction=function(curM,curRidge,layerOneDimRed,curNhTwo,deltaVec,isLogInput,isLogOutput,isExp){
  
  
  ######### Forces these values to be discrete
  curM=as.integer(curM)
  layerOneDimRed=as.integer(layerOneDimRed)
  if(!(layerOneDimRed%%2)==0 & dimRedESN=="Fourier"){
    layerOneDimRed=layerOneDimRed+1
    
  }
  
  curNhTwo=as.integer(curNhTwo)
  
  
  cat("curM ",curM, "\n")
  cat("curRidge ",curRidge, "\n")
  cat("layerOneDimRed ",layerOneDimRed, "\n")
  cat("curNhTwo ",curNhTwo, "\n")
  
  for(i in 1:length(deltaVec)){
    cat("curDelta ", i, " ",deltaVec[i], "\n")
  }
  
  
  ########Variables for validation
  curTrainLen=validLen
  curXTestIndex=xValTestIndex
  curTestLen=validTestLen
  curInbetweenIndexes=validInbetweenIndexes
  curInbetweenLen=validInbetweenLen
  curYTestLen=length(yValTestIndex)
  curYTest=rawData[yValTestIndex,]
  
  gaDataInput=rawDataInput
  if(isLogInput){
    gaDataInput=log(rawDataInput)
  }
  
  gaDataOutput=yTrain
  if(isLogOutput){
    gaDataOutput=log(yTrain)
  }
  
  # createESNDataObj=createEmbedRNNData(curTrainLen,curM,tauEmb,yTrain,log(rawDataInput),numLocs,curXTestIndex,curTestLen,curInbetweenIndexes,curInbetweenLen)
  createESNDataObj=createEmbedRNNData(curTrainLen,curM,tauEmb,gaDataOutput,gaDataInput,numLocs,curXTestIndex,curTestLen,curInbetweenIndexes,curInbetweenLen)
  
  lenInSampleEmb=createESNDataObj$lenInSampleEmb
  addScaleMat=createESNDataObj$addScaleMat
  altScaleFactor=createESNDataObj$curInSampSD
  designMatrix=createESNDataObj$designMatrix
  designMatrixOutSample=createESNDataObj$designMatrixOutSample
  altYInSamp=createESNDataObj$curYInSamp
  designMatrixInbetween=createESNDataObj$inbetweenDesignMatrix
  
  
  ##########first layer
  setParObjFirst=setParsEESN(curRidge ,curNhOne,numLocs,curM)
  sampVecESNOne=setParObjFirst$sampVecESN
  stratValuesXtempOne=setParObjFirst$stratValuesXtemp
  
  
  
  setParObjTwo=setParsEESN(curRidge ,curNhTwo,layerOneDimRed,curMTwo)
  sampVecESNTwo=setParObjTwo$sampVecESN
  stratValuesXtempTwo=setParObjTwo$stratValuesXtemp
  
  if(featureLinkInd){
    
    if(quadInd){
      ridgeMat=diag(curRidge,2*curNhTwo+(numLayers-1)*layerOneDimRed)
    }else{
      ridgeMat=diag(curRidge,curNhTwo+(numLayers-1)*layerOneDimRed)
    }
    
  }else{
    ridgeMat=setParObjTwo$ridgeMat
  }
  
  
  nColsUOne=ncol(designMatrix)
  nColsUTwo=layerOneDimRed+1
  
  forMatMultiESN=array(NA,c(ensembleLen,curYTestLen,numLocs))
  
  for(iESN in 1:ensembleLen){
    tempDeepESNFor=multiLayerDeepESNFor(curNhOne,lenInSampleEmb,designMatrix,sampVecESNOne,wWidthOne,piWESNOne,piUESNOne,uWidthOne,nColsUOne,deltaVec,stratValuesXtempOne,curTestLen,dimRedESN,layerOneDimRed,curNhTwo,sampVecESNTwo,uWidthTwo,piWESNTwo,piUESNTwo,wWidthTwo,nColsUTwo,stratValuesXtempTwo,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,designMatrixOutSample,designMatrixInbetween,curInbetweenLen,featureLinkInd,numLayers)$outSampFor
    if(isExp){
      forMatMultiESN[iESN,,]=exp(tempDeepESNFor)
    }else{
      forMatMultiESN[iESN,,]=tempDeepESNFor
    }
    
  }
  
  finalEsembelFor=sapply(1:numLocs,function(funVar) colMeans(forMatMultiESN[,,funVar]))
  deepESNMSE=mean((finalEsembelFor-curYTest)^2)
  
  cat("deepESNMSE ", deepESNMSE ,"\n")
  
  return(deepESNMSE)
}

#############################################
## Classic ESN Generic Genetic Algorithm ####
#############################################
esnGAFunction=function(curM,curRidge,curNh,curDelta,isLogInput,isLogOutput,isExp){
  
  ######### Forces these values to be discrete
  curM=as.integer(curM)
  curNh=as.integer(curNh)
  
  cat("curM ",curM, "\n")
  cat("curRidge ",curRidge, "\n")
  cat("curNh ",curNh, "\n")
  cat("curDelta ",curDelta, "\n")
  
  ########Variables for validation
  curTrainLen=validLen
  curXTestIndex=xValTestIndex
  curTestLen=validTestLen
  curInbetweenIndexes=validInbetweenIndexes
  curInbetweenLen=validInbetweenLen
  curYTestLen=length(yValTestIndex)
  curYTest=rawData[yValTestIndex,]
  
  
  gaDataInput=rawDataInput
  if(isLogInput){
    gaDataInput=log(rawDataInput)
  }
  
  
  gaDataOutput=yTrain
  if(isLogOutput){
    gaDataOutput=log(yTrain)
  }
  
  
  
  # createESNDataObj=createEmbedRNNData(curTrainLen,curM,tauEmb,yTrain,log(rawDataInput),numLocs,curXTestIndex,curTestLen,curInbetweenIndexes,curInbetweenLen)
  createESNDataObj=createEmbedRNNData(curTrainLen,curM,tauEmb,gaDataOutput,gaDataInput,numLocs,curXTestIndex,curTestLen,curInbetweenIndexes,curInbetweenLen)
  
  lenInSampleEmb=createESNDataObj$lenInSampleEmb
  addScaleMat=createESNDataObj$addScaleMat
  altScaleFactor=createESNDataObj$curInSampSD
  designMatrix=createESNDataObj$designMatrix
  designMatrixOutSample=createESNDataObj$designMatrixOutSample
  altYInSamp=createESNDataObj$curYInSamp
  designMatrixInbetween=createESNDataObj$inbetweenDesignMatrix
  
  
  setParObj=setParsEESN(curRidge ,curNh,numLocs,curM)
  
  
  
  sampVecESN=setParObj$sampVecESN
  stratValuesXtemp=setParObj$stratValuesXtemp
  if(quadInd){
    ridgeMat=diag(curRidge,2*curNh)
  }else{
    ridgeMat=diag(curRidge,curNh)
  }
  nColsU=ncol(designMatrix)
  
  forMatMultiESN=array(NA,c(ensembleLen,curYTestLen,numLocs))
  
  for(iEnsem in 1:ensembleLen){
    # scaledOutSampleFor=calcESNForecastsSpinForward(curNh,lenInSampleEmb,curTestLen,designMatrix,designMatrixOutSample,sampVecESN,wWidth,piWESN,piUESN,uWidth,nColsU,curDelta,stratValuesXtemp,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,quadInd,designMatrixInbetween,curInbetweenLen)
    scaledOutSampleFor=t(calcESNForecastsSpinForward(curNh,lenInSampleEmb,curTestLen,designMatrix,designMatrixOutSample,sampVecESN,wWidth,piWESN,piUESN,uWidth,nColsU,curDelta,stratValuesXtemp,ridgeMat,altYInSamp,altScaleFactor,addScaleMat,quadInd,designMatrixInbetween,curInbetweenLen))
    
    if(isExp){
      forMatMultiESN[iEnsem,,]=exp(scaledOutSampleFor)
    }else{
      forMatMultiESN[iEnsem,,]=scaledOutSampleFor
    }
    
  }
  
  
  finalEsembelFor=sapply(1:numLocs,function(funVar) colMeans(forMatMultiESN[,,funVar]))
  
  ESNMSE=mean((finalEsembelFor-curYTest)^2)
  cat("ESN MSE ", ESNMSE ,"\n")
  return(ESNMSE)
}