// g++ -o test1 test1.cpp



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





int main(int argc, char** argv)
{
    

    // Dimension of the input data. In our case, it should be 10=5*2 (5points underneath the foot and 2 is for position and velocity) for each model. In total, we need 6 models one per each force/torque.
    int dimX = 2;
    // Number of training data. In our case it is 150.
    int numSamplesY = 150;
    
    MatrixXd X(numSamplesY,dimX);
    VectorXd y(numSamplesY);
    VectorXd xTest(dimX);
    MatrixXd Dprecision(dimX,dimX);
    // constant parameters
    double cutoff = 0;
    double sigma = 0.1;
    int minNumPts = 2;
    
    // Position and velocity for the five points. The size should be 10 and we need 6 sets of these.
    xTest(0) = 25;
    xTest(1) = 450;
    
    for (int i=0; i<dimX; ++i) {
        Dprecision(i,i) = sigma;
    }
    
    
    double beta[dimX+1];
    double *outBetaPtr;
    outBetaPtr = &beta[0];
    
    
    // This part will be replaced by the recorded data from previous steps
    for (int i=0; i<numSamplesY; ++i) {
        X(i,0) = i;
        X(i,1) = 2*i*i;
        y(i) = 3*i;
    }
    
    // Calling the learning algorithm
    localLinearWeightedRegression(X, y, xTest, Dprecision, cutoff, minNumPts, dimX, numSamplesY, outBetaPtr);
    
    // output is y = beta0 + beta1*x1 + ... + beta10*x10
    // we need 6 sets of these
    double predict;
    predict = *outBetaPtr + *(outBetaPtr+1) * xTest(0) + *(outBetaPtr+2) * xTest(1);
    printf("predicted value ==> %f\n", predict );
    
    

   return 0; 
}