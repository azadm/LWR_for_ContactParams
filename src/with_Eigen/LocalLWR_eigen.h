//
//  LocalLWR_eigen.h
//  
//
//  Created by Morteza Azad on 24/02/2016.
//
//

#include </usr/local/include/eigen3/Eigen/Dense>
#include </usr/local/include/eigen3/Eigen/Core>
#include </usr/local/include/eigen3/Eigen/Cholesky>
#include </usr/local/include/eigen3/Eigen/Eigenvalues>
#include </usr/local/include/eigen3/Eigen/SVD>

using namespace Eigen;

double localLinearWeightedRegression(MatrixXd X, VectorXd y, VectorXd xTest, MatrixXd Dprecision, double cutoff, int minNumPts, int dimX, int numSamplesY, double *outBeta);
