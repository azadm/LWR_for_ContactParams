// g++ -c LocalLWR.cpp -I/usr/local/include/newmat10
// g++ -o test1 test1.cpp -I/usr/local/include/newmat10 -L/usr/local/include/newmat10 -lnewmat LocalLWR.o



#include <iostream>
#include <cstdio>
#include <math.h>
#include "LocalLWR_eigen.h"





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