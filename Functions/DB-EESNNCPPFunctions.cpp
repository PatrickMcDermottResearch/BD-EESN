#include <RcppArmadillo.h>
#include <Rcpp/Rmath.h>
#include <RcppTN.h>
#include <math.h>
#include <math.h>
#include <RcppEigen.h>
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;  //default from RStudio
// [[Rcpp::depends(RcppTN)]]
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppArmadillo)]]
using Eigen::Map;                 // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers



// [[Rcpp::export]]
MatrixXd tanh(const MatrixXd& x) {
  const int n(x.rows());
  const int p(x.cols());
  MatrixXd expmat(MatrixXd(n, p));
  expmat = (2 * x.array()).array().exp();
  expmat = (expmat.array() - 1) / (expmat.array() + 1);
  return (expmat);
}

// [[Rcpp::export]]
arma::vec gaussianKernel(arma::vec x){
  arma::vec returnValue=exp(-pow(x,2.0));
  return returnValue;
}

// [[Rcpp::export]]
arma::vec bipolarSig(arma::vec x){
  arma::vec returnValue=(1-exp(-x))/(1+exp(-x));
  return returnValue;
}

// [[Rcpp::export]]
double acceptSigmaUTune(double pStarU, double curCount, double sigmaVal){
  double c=sigmaVal/(pStarU*(1-pStarU));
  sigmaVal=sigmaVal+(c*(1-pStarU))/curCount;
  return sigmaVal;
}

// [[Rcpp::export]]
double rejectSigmaUTune(double pStarU, double curCount, double sigmaVal){
  double c=sigmaVal/(pStarU*(1-pStarU));
  sigmaVal=sigmaVal-(c*(pStarU))/curCount;
  return sigmaVal;
}
////function for symetrically bounded parameter spaces
// [[Rcpp::export]]
double boundedTrans(double curVal,double len,double minVal){
  double returnVal= -log(len/(curVal-minVal)-1);
  return returnVal;
}

// [[Rcpp::export]]
double boundedTransBack(double curVal,double len,double minVal){
  double returnVal= minVal+len/(1+exp(-curVal));
  return returnVal;
}

// [[Rcpp::export]]
double forwardPAAlphaTrans( double alphaVal , double aBoundVal,double wMatVal){
  
  double returnVal=alphaVal-log((2*aBoundVal)/(wMatVal+aBoundVal)-1) ;
  return returnVal;
  
  
}

// [[Rcpp::export]]
double backwardPAAlphaTrans( double alphaVal , double aBoundVal,double wMatVal){
  double returnVal= -aBoundVal+(2*aBoundVal)/(1+exp(-wMatVal+alphaVal));
  return returnVal;
}

// forwardDATrans=function(alphaVal,aVal,wMatVal){
//   alphaVal-log((2*aBoundVal)/(wMatVal+aBoundVal)-1)
//   
// }





//////////////////////////////////////////////////////////////
///////////////// Shared function for Create Functions ///////
////////////////////////////////////////////////////////////
// [[Rcpp::export]]
MatrixXd createQuadTerms(MatrixXd xTemp,NumericVector indexMat, int numQuadIntTerms){
  
  
  MatrixXd tempProdMat=xTemp*xTemp.transpose();
  Map<MatrixXd> upperTriVals(tempProdMat.data(), tempProdMat.size(),1);
  
  MatrixXd tempQuad(MatrixXd(numQuadIntTerms, 1));
  for(int i=0; i<numQuadIntTerms; ++i){
    tempQuad.row(i)=upperTriVals.row(indexMat(i)-1);
  }
  
  return tempQuad;
}

///////////////////////////////////////////////////////
///////////////// createHMat FUNCTION /////////////////
///////////////////////////////////////////////////////

// [[Rcpp::export]]
List createHMat(int nh, int trainLen,Eigen::MatrixXd designMat,Eigen::MatrixXd wMatScaled, Eigen::MatrixXd uMat,Eigen::VectorXd startValues,NumericVector indexMat,int numQuadIntTerms,bool quadIntMCMC,NumericVector cppIndexSeq,bool reducedInd ){
  
  
  ///hMat with only the rain season
  Eigen::MatrixXd hMat(MatrixXd(2*nh+1,cppIndexSeq.size())) ;
  
  
  MatrixXd xTemp(MatrixXd(nh, 1));
  xTemp=startValues;
  
  Eigen::MatrixXd hMatFull(MatrixXd(2*nh+1,trainLen).setOnes()) ;
  if(quadIntMCMC){
    hMatFull=(MatrixXd(nh+numQuadIntTerms+1,trainLen).setOnes());
  }
  
  
  Eigen::MatrixXd uProdMat=designMat*uMat.transpose();
  
  
  int rsCount=0;
  
  for(int t=0; t<trainLen; ++t){
    xTemp=tanh(wMatScaled*xTemp+uProdMat.row(t).transpose());
    hMatFull.col(t).segment(1,nh)=xTemp;
    
    
    if(quadIntMCMC){
      hMatFull.col(t).segment(nh+1,numQuadIntTerms)= createQuadTerms(xTemp,  indexMat,  numQuadIntTerms);
    }else{
      hMatFull.col(t).segment(nh+1,nh)= (xTemp.array()).array().pow(2);
    }
    
    
    //put in rs vars
    if(reducedInd && rsCount<cppIndexSeq.size() ){
      
      if(cppIndexSeq(rsCount)==t){
        
        // Rcout << "t " <<t<< std::endl;
        // Rcout << "hMatFull.col(t) " <<hMatFull.col(t)<< std::endl;
        // Rcout << "hMat.col(rsCount)" <<hMat.col(rsCount)<< std::endl;
      
        hMat.col(rsCount)=hMatFull.col(t);
       rsCount=rsCount+1;
      }
      
    }
    
    
    
  }
  
  // Rcout << "cppIndexSeq " <<cppIndexSeq.size()<< std::endl;
  
  return List::create(Named("hMatFull") = hMatFull,Named("hMat") = hMat, Named("xTempLast") = xTemp );
  
}



//////////////////////////////////////////////////////////////
///////////////// FOR W createHMatW FUNCTION ////////////////
////////////////////////////////////////////////////////////


// [[Rcpp::export]]
List createHMatW(int nh, int trainLen,Eigen::MatrixXd uProductMat,Eigen::MatrixXd wMatScaled,Eigen::VectorXd startValues,NumericVector indexMat,int numQuadIntTerms,bool quadIntMCMC,NumericVector cppIndexSeq,bool reducedInd){
  
  ///hMat with only the rain season
  Eigen::MatrixXd hMat(MatrixXd(2*nh+1,cppIndexSeq.size())) ;
  
  
  MatrixXd xTemp(MatrixXd(nh, 1));
  xTemp=startValues;
  
  Eigen::MatrixXd hMatFull(MatrixXd(2*nh+1,trainLen).setOnes()) ;
  if(quadIntMCMC){
    hMatFull=(MatrixXd(nh+numQuadIntTerms+1,trainLen).setOnes());
  }
  
  int rsCount=0;
  
  for(int t=0; t<trainLen; ++t){
    
    xTemp=tanh(wMatScaled*xTemp+uProductMat.row(t).transpose());
    hMatFull.col(t).segment(1,nh)=xTemp;
    
    if(quadIntMCMC){
      hMatFull.col(t).segment(nh+1,numQuadIntTerms)= createQuadTerms(xTemp,  indexMat,  numQuadIntTerms);
    }else{
      
      hMatFull.col(t).segment(nh+1,nh)= (xTemp.array()).array().pow(2);
    }
    
    
    //put in rs vars
    if(reducedInd && rsCount<cppIndexSeq.size() ){
      
      if(cppIndexSeq(rsCount)==t){
        
        // Rcout << "cur index " <<cppIndexSeq(rsCount)<< std::endl;
        hMat.col(rsCount)=hMatFull.col(t);
        rsCount=rsCount+1;
      }
      
    }
    
  }
  
  
  
  
  return List::create(Named("hMatFull") = hMatFull,Named("hMat") = hMat, Named("xTempLast") = xTemp  );
}







//////////////////////////////////////////////////////////////
///////////////// FOR U createHMatW FUNCTION ////////////////
////////////////////////////////////////////////////////////

// [[Rcpp::export]]
List createHMatU(int nh, int trainLen,Eigen::MatrixXd designMat,Eigen::MatrixXd wMatScaled, Eigen::MatrixXd uMat,Eigen::VectorXd startValues,Eigen::MatrixXd uProdMat, int curRow,NumericVector indexMat,int numQuadIntTerms,bool quadIntMCMC,NumericVector cppIndexSeq,bool reducedInd ){
  
  
  ///hMat with only the rain season
  Eigen::MatrixXd hMat(MatrixXd(2*nh+1,cppIndexSeq.size())) ;
  
  MatrixXd xTemp(MatrixXd(nh, 1));
  xTemp=startValues;
  
  Eigen::MatrixXd hMatFull(MatrixXd(2*nh+1,trainLen).setOnes()) ;
  if(quadIntMCMC){
    hMatFull=(MatrixXd(nh+numQuadIntTerms+1,trainLen).setOnes());
  }
  
  
  Eigen::VectorXd curProdVal ;
  MatrixXd curProdVec(MatrixXd(trainLen, 1));
  
  int rsCount=0;
  
  for(int t=0; t<trainLen; ++t){
    curProdVal= uMat.row(curRow)*designMat.row(t).transpose();
    curProdVec.row(t)=curProdVal;
    uProdMat.row(t).col(curRow)=curProdVal;
    xTemp=tanh(wMatScaled*xTemp+uProdMat.row(t).transpose());
    
    hMatFull.col(t).segment(1,nh)=xTemp;
    
    if(quadIntMCMC){
      
      hMatFull.col(t).segment(nh+1,numQuadIntTerms)= createQuadTerms(xTemp,  indexMat,  numQuadIntTerms);
      
    }else{
      
      hMatFull.col(t).segment(nh+1,nh)= (xTemp.array()).array().pow(2);
      
    }
    
    
    //put in rs vars
    if(reducedInd && rsCount<cppIndexSeq.size() ){
      
      if(cppIndexSeq(rsCount)==t){
        
        // Rcout << "cur index " <<cppIndexSeq(rsCount)<< std::endl;
        hMat.col(rsCount)=hMatFull.col(t);
        rsCount=rsCount+1;
      }
      
    }
    
  }
  
  
  
  return List::create(Named("hMatFull") = hMatFull, Named("hMat") = hMat, Named("curProdVec") = curProdVec );
}


// [[Rcpp::export]]
NumericVector convertToNumericVec(Eigen::VectorXd x){
  
  std::vector<double> v2;
  v2.resize(x.size());
  VectorXd::Map(&v2[0], x.size()) = x;
  NumericVector convertedData = wrap(v2);
  return  convertedData;
}

// [[Rcpp::export]]
int mod(double firstVal,double secVal){
  int tempVal=floor(firstVal/secVal);
  int returnVal=firstVal-tempVal*secVal;
  return returnVal;
}


// [[Rcpp::export]]
List multiLorenz40CppAlt(arma::vec x0, arma::vec y0,double theta,  double dt, int n, int M, int numSmScaleLocs, double hx,double hy, double epsilonScale, double largeScaleVar,int numSmScaleLocsTot,bool isError){
  //Chorin and Lu parameterization
  
  
  arma::vec xx = x0;
  arma::vec yy=y0;
  
  arma::vec dx=arma::zeros<arma::vec>(n);
  arma::vec dy=arma::zeros<arma::vec>(numSmScaleLocsTot);
  
  
  for(int j=0; j<M; ++j){
    
    
    
    int smallCounter=0;
    for(int i=0; i<n; ++i){
      
      arma::vec curYVals=yy.rows(i*numSmScaleLocs,i*numSmScaleLocs+numSmScaleLocs-1);
      
      if(isError){
        dx(i) = ( xx(mod(i-1,n))*(xx(mod(i+1,n)) - xx(mod(i-2,n))) - xx(i) + theta +(hx /numSmScaleLocs)*arma::sum(curYVals)  +R::rnorm(0.00,sqrt(largeScaleVar))   ) * dt ;
      }else{
        dx(i) = ( xx(mod(i-1,n))*(xx(mod(i+1,n)) - xx(mod(i-2,n))) - xx(i) + theta +(hx /numSmScaleLocs)*arma::sum(curYVals) ) * dt ;
      }
      
      for(int k=0; k<numSmScaleLocs; ++k){
        
        
        dy(smallCounter)=(  (1/epsilonScale)*(yy(mod(smallCounter+1,numSmScaleLocsTot))*( yy(mod(smallCounter-1,numSmScaleLocsTot))-yy(mod(smallCounter+2,numSmScaleLocsTot)) )-yy(mod(smallCounter,numSmScaleLocsTot)) +hy*xx(i)  ))*dt ;
        smallCounter=smallCounter+1;
      }
      
    }
    
    xx = xx + dx ;
    yy=  yy + dy ;
  }
  return List::create(Named("xx") = xx, Named("yy") = yy );
}





// [[Rcpp::export]]
List multiLorenz40Cpp(arma::vec x0, arma::vec y0,double theta,  double dt, int n, int M, int numSmScaleLocs, double hParm,double cParm, double bParm, double largeScaleVar,int numSmScaleLocsTot){
  //Chorin and Lu parameterization
  
  
  arma::vec xx = x0;
  arma::vec yy=y0;
  
  arma::vec dx=arma::zeros<arma::vec>(n);
  arma::vec dy=arma::zeros<arma::vec>(numSmScaleLocsTot);
  
  
  for(int j=0; j<M; ++j){
    
    
    
    int smallCounter=0;
    for(int i=0; i<n; ++i){
      
      arma::vec curYVals=yy.rows(i*numSmScaleLocs,i*numSmScaleLocs+numSmScaleLocs-1);
      
      dx(i) = ( -xx(mod(i-1,n))*(xx(mod(i-2,n)) - xx(mod(i+1,n))) - xx(i) + theta -((hParm*cParm) /bParm)*arma::sum(curYVals)+R::rnorm(0.00,sqrt(largeScaleVar))    ) * dt ;
      
      for(int k=0; k<numSmScaleLocs; ++k){
        dy(smallCounter)=((-cParm*bParm*yy(mod(smallCounter+1,numSmScaleLocsTot)))*( yy(mod(smallCounter+2,numSmScaleLocsTot))-yy(mod(smallCounter-1,numSmScaleLocsTot)) )-cParm*yy(mod(smallCounter,numSmScaleLocsTot)) +((cParm*hParm)/bParm)*xx(i)  )*dt ;
        smallCounter=smallCounter+1;
      }
      
    }
    
    xx = xx + dx ;
    yy=  yy + dy ;
  }
  return List::create(Named("xx") = xx, Named("yy") = yy );
}



////////////////////////////////////////////////////////////////
///////////////// Multiscale RNN Hidden State Functions ///////
//////////////////////////////////////////////////////////////

// [[Rcpp::export]]
NumericVector csample_num( NumericVector x,
                           int size,
                           NumericVector prob = NumericVector::create()
) {
  NumericVector ret = RcppArmadillo::sample(x, size, FALSE, prob);
  return ret;
}


// [[Rcpp::export]]
NumericMatrix genRandMat(double numRow,double numCol,double width){
  Rcpp::NumericVector draws = Rcpp::runif(numRow*numCol,-width,width);
  
  return NumericMatrix(numRow, numCol, draws.begin());
}


///////////////////////////////////////////////////////
///////////////// ESN createHMat FUNCTION /////////////////
//////////////////////////////////////////////////////


// [[Rcpp::export]]
List esnCreateHMat(int nh, int trainLen,Eigen::MatrixXd wMatScaled, Eigen::MatrixXd uProdMat,Eigen::VectorXd startValues,bool quadInd,int hMatdim){
  
  MatrixXd xTemp(MatrixXd(nh, 1));
  xTemp=startValues;
  
  
  Eigen::MatrixXd hMat(MatrixXd(hMatdim,trainLen).setZero()) ;
  
  for(int t=0; t<trainLen; ++t){
    // Rcout << "t " <<t<< std::endl;
    // Rcout << "r1 " <<uProdMat.col(t).transpose()<< std::endl;
    
    xTemp=tanh(wMatScaled*xTemp+uProdMat.col(t));
    
    // Rcout << "t " <<t<< std::endl;
    // Rcout << "xTemp " <<xTemp<< std::endl;
    
    hMat.col(t).segment(0,nh)=xTemp;
    if(quadInd){
      hMat.col(t).segment(nh,nh)= (xTemp.array()).array().pow(2);
    }
    
    
  }
  
  return List::create(Named("hMat") = hMat, Named("xTempLast") = xTemp);
}


////////////////////////////////////////////////////////////////
///////////////// generate reservoir FUNCTION /////////////////
//////////////////////////////////////////////////////////////


// List genResCPP(double nh,double echoTrainLen, Eigen::MatrixXd designMat, NumericVector sampVecESN,double wWidth, double piWESN,double piUESN,double uWidth,double nColsU, double deltaESN,Eigen::VectorXd startValues,bool quadInd,int hMatdim, Eigen::MatrixXd designMatrixOutSample, double testLen){
//   // List genResCPP(double nh,double echoTrainLen, Eigen::MatrixXd designMat, NumericVector sampVecESN, double numNonZeroW,Eigen::MatrixXd wTempESN,double piUESN,Eigen::MatrixXd uTempESN,double nColsU, double deltaESN,Eigen::VectorXd startValues){
//   
//   //Generate W
//   NumericMatrix wTempESNNumMat=genRandMat(nh,nh,wWidth);
//   Eigen::Map<Eigen::MatrixXd> wTempESN = as<Eigen::Map<Eigen::MatrixXd> >(wTempESNNumMat);
//   
//   //Generate U
//   NumericMatrix uTempESNNumMat=genRandMat(nh,nColsU,uWidth);
//   Eigen::Map<Eigen::MatrixXd> uTempESN = as<Eigen::Map<Eigen::MatrixXd> >(uTempESNNumMat);
//   
//   //Make W Matrix Sparse 
//   for(int i=0; i<nh; ++i){
//     
//     
//     NumericVector numZeroW=Rcpp::rbinom(1,nh,piWESN);
//     double numNonZeroW=nh-numZeroW(0);
//     NumericVector tempIndexW=csample_num(sampVecESN,numNonZeroW);
//     for(int j=0; j<numNonZeroW; ++j){
//       wTempESN(tempIndexW(j),i)=0;
//     }
//     
//   }
//   
//   
//   //Make U Matrix Sparse 
//   for(int i=0; i<nColsU; ++i){
//     
//     NumericVector numZeroU=Rcpp::rbinom(1,nh,piUESN);
//     double numNonZeroU=nh-numZeroU(0);
//     NumericVector tempIndexU=csample_num(sampVecESN,numNonZeroU);
//     for(int j=0; j<numNonZeroU; ++j){
//       uTempESN(tempIndexU(j),i)=0;
//     }
//   }
//   
//   //Scale W Matrix
//   double spectralRadius=wTempESN.eigenvalues().cwiseAbs().maxCoeff();
//   Eigen::MatrixXd  wMatScaled=wTempESN*deltaESN/spectralRadius;
//   
//   //Create H Matrix
//   
//   Eigen::MatrixXd uProdMat=uTempESN*designMat.transpose();
//   
//   List createHMatObj=esnCreateHMat( nh,  echoTrainLen, wMatScaled,  uProdMat, startValues,quadInd,hMatdim);
//   
//   Eigen::VectorXd xTemp=createHMatObj["xTempLast"];
//   Eigen::MatrixXd   hMat=createHMatObj["hMat"];
//   
//   
// 
//   //Out of Sample forecasts
//   Eigen::MatrixXd uProdMatOutSamp=uTempESN*designMatrixOutSample.transpose();
//   List outSampHMatObj= esnCreateHMat( nh,  testLen, wMatScaled,  uProdMatOutSamp,xTemp,quadInd,hMatdim);
// 
//   Eigen::MatrixXd   hMatOutSamp=outSampHMatObj["hMat"];
// 
//   
//   
//   
//   
//   return List::create(Named("hMat") = hMat, Named("xTempLast") = xTemp,Named("wMatESNFun") = wMatScaled,Named("uMatESNFun") = uTempESN,Named("hMatOutSamp") = hMatOutSamp);
// }


// [[Rcpp::export]]
List genResCPP(double nh,double echoTrainLen, Eigen::MatrixXd designMat, NumericVector sampVecESN,double wWidth, double piWESN,double piUESN,double uWidth,double nColsU, double deltaESN,Eigen::VectorXd startValues,bool quadInd,int hMatdim, Eigen::MatrixXd designMatrixOutSample ,double testLen,Eigen::MatrixXd designMatrixInbetween, double inbetweenLen){
  // List genResCPP(double nh,double echoTrainLen, Eigen::MatrixXd designMat, NumericVector sampVecESN, double numNonZeroW,Eigen::MatrixXd wTempESN,double piUESN,Eigen::MatrixXd uTempESN,double nColsU, double deltaESN,Eigen::VectorXd startValues){
  
  //Generate W
  NumericMatrix wTempESNNumMat=genRandMat(nh,nh,wWidth);
  Eigen::Map<Eigen::MatrixXd> wTempESN = as<Eigen::Map<Eigen::MatrixXd> >(wTempESNNumMat);
  
  //Generate U
  NumericMatrix uTempESNNumMat=genRandMat(nh,nColsU,uWidth);
  Eigen::Map<Eigen::MatrixXd> uTempESN = as<Eigen::Map<Eigen::MatrixXd> >(uTempESNNumMat);
  
  //Make W Matrix Sparse 
  for(int i=0; i<nh; ++i){
    
    
    NumericVector numZeroW=Rcpp::rbinom(1,nh,piWESN);
    double numNonZeroW=nh-numZeroW(0);
    NumericVector tempIndexW=csample_num(sampVecESN,numNonZeroW);
    for(int j=0; j<numNonZeroW; ++j){
      wTempESN(tempIndexW(j),i)=0;
    }
    
  }
  
  
  //Make U Matrix Sparse 
  for(int i=0; i<nColsU; ++i){
    
    NumericVector numZeroU=Rcpp::rbinom(1,nh,piUESN);
    double numNonZeroU=nh-numZeroU(0);
    NumericVector tempIndexU=csample_num(sampVecESN,numNonZeroU);
    for(int j=0; j<numNonZeroU; ++j){
      uTempESN(tempIndexU(j),i)=0;
    }
  }
  
  //Scale W Matrix
  double spectralRadius=wTempESN.eigenvalues().cwiseAbs().maxCoeff();
  Eigen::MatrixXd  wMatScaled=wTempESN*deltaESN/spectralRadius;
  
  //Create H Matrix
  
  Eigen::MatrixXd uProdMat=uTempESN*designMat.transpose();
  
  List createHMatObj=esnCreateHMat( nh,  echoTrainLen, wMatScaled,  uProdMat, startValues,quadInd,hMatdim);
  
  Eigen::VectorXd xTemp=createHMatObj["xTempLast"];
  Eigen::MatrixXd   hMat=createHMatObj["hMat"];
  
  
  
  //Inbetween forecasts
  Eigen::MatrixXd uInbetweenProdMat=uTempESN*designMatrixInbetween.transpose();
  List inbetweenHObj= esnCreateHMat( nh,  inbetweenLen, wMatScaled,  uInbetweenProdMat,xTemp,quadInd,hMatdim);
  Eigen::MatrixXd   inbetweenHMat=inbetweenHObj["hMat"];
  Eigen::VectorXd finalHVec=inbetweenHObj["xTempLast"];
  
  
  //Out of Sample forecasts
  Eigen::MatrixXd uProdMatOutSamp=uTempESN*designMatrixOutSample.transpose();
  List outSampHMatObj= esnCreateHMat( nh,  testLen, wMatScaled,  uProdMatOutSamp,finalHVec,quadInd,hMatdim);
  
  Eigen::MatrixXd   hMatOutSamp=outSampHMatObj["hMat"];
  
  
  
  
  return List::create(Named("inbetweenHMat") = inbetweenHMat,Named("hMat") = hMat, Named("xTempLast") = xTemp,Named("wMatESNFun") = wMatScaled,Named("uMatESNFun") = uTempESN,Named("hMatOutSamp") = hMatOutSamp);
}



////////////////////////////////////////////////////////////////
///////////////// DBayesEESNHelper  ///////////////////////////
//////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Eigen::MatrixXd  DBayesEESNHelper(int numVPars, Eigen::MatrixXd curVPars, Eigen::MatrixXd hMat, NumericVector indexMat ){
  
  Eigen::MatrixXd sumMat(MatrixXd(hMat.rows(),1).setZero());
  
  
  for(int i=0; i<indexMat.length(); ++i){
    sumMat=sumMat+hMat.block(0,(indexMat(i)-1)*numVPars, hMat.rows() ,numVPars)*curVPars.block((indexMat(i)-1)*numVPars,0,numVPars,1);
    
    
    // Eigen::MatrixXd temp=   hMat.block(0,(indexMat(i)-1)*numVPars, hMat.rows() ,numVPars);
    // Rcout << "rows1 " <<temp.rows()<< std::endl;
    // Rcout << "cols1 " <<temp.cols()<< std::endl;
    // 
    // Eigen::MatrixXd temp2=curVPars.block((indexMat(i)-1)*numVPars,0,numVPars,1);
    // 
    // Rcout << "rows2 " <<temp2.rows()<< std::endl;
    // Rcout << "cols2 " <<temp2.cols()<< std::endl;
    
    
  }
  
  return sumMat;
  
}



/////////////////////////////////////////////////////////////////////
///////////////// multiDBayesEESNHelper  ///////////////////////////
////////////////////////////////////////////////////////////////////
// [[Rcpp::export]]
Eigen::MatrixXd  multiDBayesEESNHelper(int numVPars, Eigen::MatrixXd curVPars, Eigen::MatrixXd hMat, NumericVector indexMat ){
  
  Eigen::MatrixXd sumMat(MatrixXd(hMat.rows(),curVPars.rows()).setZero());
  
  
  for(int i=0; i<indexMat.length(); ++i){
    sumMat=sumMat+hMat.block(0,(indexMat(i)-1)*numVPars, hMat.rows() ,numVPars)*curVPars.block(0,(indexMat(i)-1)*numVPars,curVPars.rows(),numVPars).transpose();
    
    
    // Eigen::MatrixXd temp=   hMat.block(0,(indexMat(i)-1)*numVPars, hMat.rows() ,numVPars);
    // Rcout << "rows1 " <<temp.rows()<< std::endl;
    // Rcout << "cols1 " <<temp.cols()<< std::endl;
    // 
    // Eigen::MatrixXd temp2=curVPars.block(0,(indexMat(i)-1)*numVPars,curVPars.rows(),numVPars);
    // 
    // Rcout << "rows2 " <<temp2.rows()<< std::endl;
    // Rcout << "cols2 " <<temp2.cols()<< std::endl;
    
    
  }
  
  return sumMat;
  
}


