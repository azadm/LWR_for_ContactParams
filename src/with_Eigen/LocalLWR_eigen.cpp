//
//  LocalLWR_eigen.cpp
//  
//
//  Created by Morteza Azad on 24/02/2016.
//
//

#include "LocalLWR_eigen.h"


#include <iostream>
#include <cstdio>
#include <math.h>
#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/local/include/eigen3/Eigen/Core>
#include </usr/local/include/eigen3/Eigen/Cholesky>
#include </usr/local/include/eigen3/Eigen/Eigenvalues>
#include </usr/local/include/eigen3/Eigen/SVD>
//#include <eigen3/Eigen/Dense>

using namespace Eigen;


double localLinearWeightedRegression(MatrixXd X, VectorXd y, VectorXd xTest, MatrixXd Dprecision, double cutoff, int minNumPts, int dimX, int numSamplesY, double *outBeta)
{
    
    double yPred = 888.888;
    
    MatrixXd xDiff(numSamplesY,dimX);
    
    int xId, yId;
    for (yId=0; yId<numSamplesY; yId++) {
        xDiff.row(yId) = X.row(yId) - xTest.transpose();
    }
    
    
    MatrixXd W;
    W = xDiff*Dprecision*xDiff.transpose();
    
    /*
     int numValid = 0; int numWhileIts = 0;
     while(numValid < minNumPts && numWhileIts <= 10){
     numValid = 0;
     printf("%d\n", numValid);
     
     for(yId = 0; yId < numSamplesY; yId++){
     if(numWhileIts < 1){
     W( yId, yId) = exp(-0.5 * W( yId, yId ));
     }
     if(W( yId, yId ) > cutoff){
     numValid++;
     }
     }
     cutoff = cutoff * 0.5;
     numWhileIts++;
     }
     printf("%d\n", numValid);
     
     cutoff = cutoff / 0.5;
     numWhileIts--;
     
     if(numValid < 2){
     //printf("WARN: Less than TWO neighbors found\n");
     return 0.0;
     }
     */
    
    // *** because cutoff = 0, in case cutoff > 0, use code above.
    for(yId = 0; yId < numSamplesY; yId++){
        W( yId, yId ) = exp(-0.5 * W( yId, yId ));
    }
    int numValid = numSamplesY;
    
    
    MatrixXd XSubSet(numValid, dimX+1);
    VectorXd ySubSet(numValid);
    MatrixXd WSubSet(numValid, numValid);
    
    int subSetId=0;
    for(yId = 0; yId < numSamplesY; yId++){
        //if(W( yId, yId ) > cutoff){
        XSubSet( subSetId , 0 ) = 1.0;   //offset
        for(xId = 1; xId < dimX+1; xId++){
            XSubSet( subSetId , xId ) = X( yId , xId-1 );
        }
        ySubSet( subSetId ) = y( yId );
        WSubSet( subSetId ) = sqrt( W( yId, yId ) ) + 1e-6; //reg term
        subSetId++;
        // }
    }
    
    
    
    VectorXd beta;
    MatrixXd XWX;
    VectorXd XWy;
    MatrixXd XSubSetTmp(numValid-2, dimX+1);
    VectorXd ySubSetTmp(numValid-2);
    MatrixXd WSubSetTmp(numValid-2, numValid-2);
    
    
    //remove two outlier
    if(1){
        
        if(numValid < 4){
            //printf("WARN: Cannot remove outlier, Less than Four neighbors found\n");
            fflush(stdout);
            XWX =  XSubSet.transpose() * WSubSet * XSubSet;
            XWy = XSubSet.transpose() * WSubSet * ySubSet;
        } else {
            int minIdY, maxIdY;
            double minY = ySubSet.minCoeff();
            double maxY = ySubSet.maxCoeff();
            
            if (minIdY == maxIdY){
                if(minIdY == 0){
                    maxIdY =1;
                }
                else{
                    maxIdY = 0;
                }
            }
            
            int tmpId = 0;
            for(yId = 0; yId < numValid; yId++){
                if(yId != minIdY && yId != maxIdY){
                    ySubSetTmp( tmpId ) = ySubSet( yId );
                    WSubSetTmp( tmpId, tmpId ) = WSubSet( yId, yId );
                    for(xId = 0; xId < dimX+1; xId++){
                        XSubSetTmp( tmpId , xId ) = XSubSet( yId , xId );
                    }
                    tmpId++;
                }
            }
            XWX =  XSubSetTmp.transpose() * WSubSetTmp * XSubSetTmp;
            XWy = XSubSetTmp.transpose() * WSubSetTmp * ySubSetTmp;
        }
    }
    
    
    
    if(0){
        // Cholesky decomposition of SSQ
        try{
            MatrixXd L = XWX.llt().matrixL();
            beta = L.transpose().inverse() * (L.inverse() * XWy);
        }
        catch(...){
            printf("WARN: Cholesky decomposition failed\n");
            return 0.0;
        }
    }else{
        
        try{
            //Eigenvalues decomposition
            //MatrixXd V, D;
            //EigenSolver<MatrixXd> es(XWX);
            //V = es.eigenvectors(); // XWX = V * D * V.transpose(), XWX.inverse() = V * D.inverse() * V.transpose();
            //D = es.eigenvalues().asDiagonal();
            ////V = XWX.eigensolver().eigenvectors();
            //beta = V * D.inverse() * V.transpose() * XWy;
            
            //SVD decomposition
            MatrixXd U, V;
            VectorXd SingularV;
            JacobiSVD<MatrixXd> svd(XWX, ComputeThinU | ComputeThinV);
            U = svd.matrixU();
            SingularV = svd.singularValues();
            V = svd.matrixV(); // XWX = U * S * V.transpose()
            MatrixXd S;
            S.setZero(U.cols(),V.cols());
            for (int i=0; i<SingularV.rows(); i++) {
                S(i,i) = SingularV(i);
            }
            beta = U * S.inverse() * V.transpose() * XWy;
        }
        catch(...){
            printf("WARN: Matrix decomposition / inverse failed\n");
            return 0.0;
        }
    }
    
    
    yPred = beta( 0 );
    outBeta[0] = yPred;
    
    for(xId = 0; xId < dimX; xId++){
        yPred += xTest( xId ) * beta( xId + 1);
        outBeta[xId + 1] = beta( xId + 1);
    }
    
    return yPred;
}




