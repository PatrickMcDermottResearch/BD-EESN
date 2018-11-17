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
using Eigen::SparseMatrix;


// Eigen::VectorXd tanH(Eigen::VectorXd x){
//   
// 
//   std::vector<double> v2;
//   v2.resize(x.size());
//   VectorXd::Map(&v2[0], x.size()) = x;
//   NumericVector convertedData = wrap(v2);
//   NumericVector transformedData=tanh(convertedData);
//   Eigen::Map<Eigen::VectorXd> returnValue(Rcpp::as<Eigen::Map<Eigen::VectorXd> >(transformedData));
//  
//   return returnValue;
// }


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
///////////////// createHMat FUNCTION /////////////////
//////////////////////////////////////////////////////


// [[Rcpp::export]]
List createHMat(int nh, int trainLen,Eigen::MatrixXd wMatScaled, Eigen::MatrixXd uProdMat,Eigen::VectorXd startValues,bool quadInd,int hMatdim){
  
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

// [[Rcpp::export]]
List genResCPP(double nh,double echoTrainLen, Eigen::MatrixXd designMat, NumericVector sampVecESN,double wWidth, double piWESN,double piUESN,double uWidth,double nColsU, double deltaESN,Eigen::VectorXd startValues,bool quadInd,int hMatdim){
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

  List createHMatObj=createHMat( nh,  echoTrainLen, wMatScaled,  uProdMat, startValues,quadInd,hMatdim);
    
    Eigen::VectorXd xTemp=createHMatObj["xTempLast"];
    Eigen::MatrixXd   hMat=createHMatObj["hMat"];

  
  return List::create(Named("hMat") = hMat, Named("xTempLast") = xTemp,Named("wMatESNFun") = wMatScaled,Named("uMatESNFun") = uTempESN);
}





////////////////////////////////////////////////////////////////
///////////////// Ensemble ESN FUNCTION ///////////////////////
//////////////////////////////////////////////////////////////
// [[Rcpp::export]]
MatrixXd calcESNForecasts(double nh,double echoTrainLen,double testLen, Eigen::MatrixXd designMat,Eigen::MatrixXd designMatrixOutSample, NumericVector sampVecESN,double wWidth, double piWESN,double piUESN,double uWidth,double nColsU, double deltaESN,Eigen::VectorXd startValues,Eigen::MatrixXd ridgeMat,Eigen::MatrixXd yInSampMat,double altScaleFactor,Eigen::MatrixXd addScaleMat,bool quadInd){

  int hMatdim=0;
  if(quadInd){
    hMatdim=2*nh;
  }else{
    hMatdim=nh;
  }
  

    //Generate the Reservoir
    List reservoirObj=genResCPP(nh,echoTrainLen,designMat,sampVecESN,wWidth,piWESN,piUESN,uWidth,nColsU,deltaESN,startValues,quadInd,hMatdim);

    Eigen::VectorXd xTemp=reservoirObj["xTempLast"];
    Eigen::MatrixXd   hMat=reservoirObj["hMat"];
    Eigen::MatrixXd   wMatScaled=reservoirObj["wMatESNFun"];
    Eigen::MatrixXd   uTempESN=reservoirObj["uMatESNFun"];

    //Ridge Regression to get vMat
    Eigen::MatrixXd vMatESN=yInSampMat.transpose()*hMat.transpose()*(hMat*hMat.transpose()+ridgeMat).inverse();

    //Out of Sample forecasts
    Eigen::MatrixXd uProdMatOutSamp=uTempESN*designMatrixOutSample.transpose();
    List outSampHMatObj= createHMat( nh,  testLen, wMatScaled,  uProdMatOutSamp,xTemp,quadInd,hMatdim);

    Eigen::MatrixXd   hMatOutSamp=outSampHMatObj["hMat"];

    Eigen::MatrixXd outSampFor=altScaleFactor*(vMatESN*hMatOutSamp)+addScaleMat;


  return outSampFor;
}




// [[Rcpp::export]]
MatrixXd calcESNForecastsSpinForward(double nh,double echoTrainLen,double testLen, Eigen::MatrixXd designMat,Eigen::MatrixXd designMatrixOutSample, NumericVector sampVecESN,double wWidth, double piWESN,double piUESN,double uWidth,double nColsU, double deltaESN,Eigen::VectorXd startValues,Eigen::MatrixXd ridgeMat,Eigen::MatrixXd yInSampMat,double altScaleFactor,Eigen::MatrixXd addScaleMat,bool quadInd,Eigen::MatrixXd inbetweenMat, double inbetweenLen){

  int hMatdim=0;
  if(quadInd){
    hMatdim=2*nh;
  }else{
    hMatdim=nh;
  }


  //Generate the Reservoir
  List reservoirObj=genResCPP(nh,echoTrainLen,designMat,sampVecESN,wWidth,piWESN,piUESN,uWidth,nColsU,deltaESN,startValues,quadInd,hMatdim);

  Eigen::VectorXd xTemp=reservoirObj["xTempLast"];
  Eigen::MatrixXd   hMat=reservoirObj["hMat"];
  Eigen::MatrixXd   wMatScaled=reservoirObj["wMatESNFun"];
  Eigen::MatrixXd   uTempESN=reservoirObj["uMatESNFun"];

  //Ridge Regression to get vMat
  Eigen::MatrixXd vMatESN=yInSampMat.transpose()*hMat.transpose()*(hMat*hMat.transpose()+ridgeMat).inverse();


  //inbetween

    Eigen::MatrixXd uProdMatInbetween=uTempESN*inbetweenMat.transpose();
    List inbetweenObj= createHMat( nh,  inbetweenLen, wMatScaled,  uProdMatInbetween,xTemp,quadInd,hMatdim);

    Eigen::VectorXd xTempFinal=inbetweenObj["xTempLast"];

  //Out of Sample forecasts
  Eigen::MatrixXd uProdMatOutSamp=uTempESN*designMatrixOutSample.transpose();
  List outSampHMatObj= createHMat( nh,  testLen, wMatScaled,  uProdMatOutSamp,xTempFinal,quadInd,hMatdim);

  Eigen::MatrixXd   hMatOutSamp=outSampHMatObj["hMat"];

  Eigen::MatrixXd outSampFor=altScaleFactor*(vMatESN*hMatOutSamp)+addScaleMat;


  return outSampFor;
}


// [[Rcpp::export]]
////Same exact function as calcESNForecastsSpinForward, except it returns multiple objs
List multiObjCalcESNForecastsSpinForward(double nh,double echoTrainLen,double testLen, Eigen::MatrixXd designMat,Eigen::MatrixXd designMatrixOutSample, NumericVector sampVecESN,double wWidth, double piWESN,double piUESN,double uWidth,double nColsU, double deltaESN,Eigen::VectorXd startValues,Eigen::MatrixXd ridgeMat,Eigen::MatrixXd yInSampMat,double altScaleFactor,Eigen::MatrixXd addScaleMat,bool quadInd,Eigen::MatrixXd inbetweenMat, double inbetweenLen){
  
  int hMatdim=0;
  if(quadInd){
    hMatdim=2*nh;
  }else{
    hMatdim=nh;
  }
  
  
  //Generate the Reservoir
  List reservoirObj=genResCPP(nh,echoTrainLen,designMat,sampVecESN,wWidth,piWESN,piUESN,uWidth,nColsU,deltaESN,startValues,quadInd,hMatdim);
  
  Eigen::VectorXd xTemp=reservoirObj["xTempLast"];
  Eigen::MatrixXd   hMat=reservoirObj["hMat"];
  Eigen::MatrixXd   wMatScaled=reservoirObj["wMatESNFun"];
  Eigen::MatrixXd   uTempESN=reservoirObj["uMatESNFun"];
  
  //Ridge Regression to get vMat
  Eigen::MatrixXd vMatESN=yInSampMat.transpose()*hMat.transpose()*(hMat*hMat.transpose()+ridgeMat).inverse();
  
  
  //inbetween
  
  Eigen::MatrixXd uProdMatInbetween=uTempESN*inbetweenMat.transpose();
  List inbetweenObj= createHMat( nh,  inbetweenLen, wMatScaled,  uProdMatInbetween,xTemp,quadInd,hMatdim);
  
  Eigen::VectorXd xTempFinal=inbetweenObj["xTempLast"];
  
  //Out of Sample forecasts
  Eigen::MatrixXd uProdMatOutSamp=uTempESN*designMatrixOutSample.transpose();
  List outSampHMatObj= createHMat( nh,  testLen, wMatScaled,  uProdMatOutSamp,xTempFinal,quadInd,hMatdim);
  
  Eigen::MatrixXd   hMatOutSamp=outSampHMatObj["hMat"];
  
  Eigen::MatrixXd outSampFor=altScaleFactor*(vMatESN*hMatOutSamp)+addScaleMat;
  
  
  return List::create(Named("hMat") = hMat, Named("outSampFor") = outSampFor,Named("hMatOutSamp") = hMatOutSamp);
}