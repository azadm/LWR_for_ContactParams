//
//  LocalLWR.h
//  
//
//  Created by Morteza Azad on 17/02/2016.
//
//

#include <newmat.h>
#include <newmatap.h>
#include "math.h"

#if defined use_namespace
using namespace NEWMAT;
#endif

double localLinearWeightedRegression(Matrix X, ColumnVector y,
                                     ColumnVector xTest, DiagonalMatrix Dprecision, double cutoff, int minNumPts, int dimX, int numSamplesY, double *outBeta);

