//
//  LocalLWR.cpp
//  
//
//  Created by Morteza Azad on 17/02/2016.
//
//

#include "LocalLWR.h"


#include <iostream>
#include <cstdio>
#include "math.h"

double localLinearWeightedRegression(Matrix X, ColumnVector y,
                                     ColumnVector xTest, DiagonalMatrix Dprecision, double cutoff, int minNumPts, int dimX, int numSamplesY, double *outBeta)
{
    double yPred = 888.888;
    
    //     printf("cutoff %f, minNumPts %i\n",cutoff, minNumPts);
    //     return yPred;
    
    Matrix xDiff(numSamplesY,dimX);
    //xDiff = bsxfun(@minus, x, xTest(i,:));
    int xId, yId;
    for(yId = 0; yId < numSamplesY; yId++){
        for(xId = 0; xId < dimX; xId++){
            //does not work xDiff.Column( xId ) = X.Column( 1 + xId ) - xTest.element( xId );   //Newmat convention for submatrices 1...numCols
            xDiff.element( yId , xId ) = X.element( yId , xId ) - xTest.element( xId );
        }
    }
    
    //     //DEBUG
    //     printf("XDiff=\n");
    //     for(yId = 0; yId < numSamplesY; yId++){
    //         for(xId = 0; xId < dimX; xId++){
    //             printf("%f ", xDiff.element( yId , xId ) );
    //         }
    //         printf("\n");
    //     }
    //     return yPred;
    
    //printf("Matrix Size %i x %i\n",xDiff.Nrows(), xDiff.Ncols());
    //printf("Matrix Size %i x %i\n",Dprecision.Nrows(), Dprecision.Ncols());
    DiagonalMatrix W;
    W << xDiff*Dprecision*xDiff.t();
    //     //DEBUG
    //     printf("XDiff*Dprec*XDiff.trans=\n");
    //     for(yId = 0; yId < numSamplesY; yId++){
    //         printf("%f ", W.element( yId ) );
    //     }
    //     return yPred;
    
    
    //     int numValid = 0; int numWhileIts = 0;
    //     while(numValid < minNumPts && numWhileIts <= 10){
    //         numValid = 0;
    //
    //         for(yId = 0; yId < numSamplesY; yId++){
    //             if(numWhileIts < 1){
    //                 W.element( yId ) = exp(-0.5 * W.element( yId ));
    //                // printf("%.8f \n", W.element( yId ) );
    //             }
    //
    //             if(W.element( yId ) > cutoff){
    //                 numValid++;
    //             }
    //         }
    //         cutoff = cutoff * 0.5;
    //         numWhileIts++;
    //     }
    //     cutoff = cutoff / 0.5;
    //     numWhileIts--;
    //
    //     if(numValid < 2){
    //         //printf("WARN: Less than TWO neighbors found\n");
    //         return 0.0;
    //     }
    
    
    // *** because cutoff = 0, in case cutoff > 0, use code above.
    for(yId = 0; yId < numSamplesY; yId++){
        W.element( yId ) = exp(-0.5 * W.element( yId ));
        //outSupport[yId] = W.element( yId );
        // printf("%.8f \n", W.element( yId ) );
    }
    int numValid = numSamplesY;
    
    
    
    //     //DEBUG
    //     printf("exp(-.5 W )=\n");
    //     for(yId = 0; yId < numSamplesY; yId++){
    //         printf("%f ", W.element( yId ) );
    //     }
    //     return yPred;
    //
    //     printf("Number of valid points %i, num. while its %i\n",numValid, numWhileIts);
    //     return yPred;
    //
    
    Matrix XSubSet(numValid, dimX+1);           //= [ones(length(valid_ids),1) x(valid_ids,:)];
    ColumnVector ySubSet(numValid);             //= y(valid_ids);
    DiagonalMatrix WSubSet(numValid);
    
    int subSetId=0;
    for(yId = 0; yId < numSamplesY; yId++){
        //if(W.element( yId ) > cutoff){
        XSubSet.element( subSetId , 0 ) = 1.0;   //offset
        for(xId = 1; xId < dimX+1; xId++){
            XSubSet.element( subSetId , xId ) = X.element( yId , xId-1 );
        }
        ySubSet.element( subSetId ) = y.element( yId );
        WSubSet.element( subSetId ) = sqrt( W.element( yId ) ) + 1e-6; //reg term
        subSetId++;
        // }
    }
    
    
    ColumnVector beta;
    SymmetricMatrix XWX;
    ColumnVector XWy;
    Matrix XSubSetTmp(numValid-2, dimX+1);           //= [ones(length(valid_ids),1) x(valid_ids,:)];
    ColumnVector ySubSetTmp(numValid-2);             //= y(valid_ids);
    DiagonalMatrix WSubSetTmp(numValid-2);
    
    
    //remove two outlier
    if(1){
        
        if(numValid < 4){
            //printf("WARN: Cannot remove outlier, Less than Four neighbors found\n");
            fflush(stdout);
            XWX <<  XSubSet.t() * WSubSet * XSubSet;
            XWy << XSubSet.t() * WSubSet * ySubSet;
        } else {
            int minIdY, maxIdY;
            ySubSet.Minimum1(minIdY);minIdY=minIdY-1;//newmat conv. starts with 1
            ySubSet.Maximum1(maxIdY);maxIdY=maxIdY-1;
            //             double minY = ySubSet.Minimum();
            //             double maxY = ySubSet.Maximum();
            //             printf("min Y %f max Y %f \n", minY, maxY );
            
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
                    //if( ySubSet.element( yId ) > minY & ySubSet.element( yId ) < maxY ) {
                    ySubSetTmp.element( tmpId ) = ySubSet.element( yId );
                    WSubSetTmp.element( tmpId ) = WSubSet.element( yId );
                    for(xId = 0; xId < dimX+1; xId++){
                        XSubSetTmp.element( tmpId , xId ) = XSubSet.element( yId , xId );
                    }
                    tmpId++;
                }
            }
            //             printf("numValid=%i numUsed=%i\n",numValid, tmpId);
            XWX <<  XSubSetTmp.t() * WSubSetTmp * XSubSetTmp;
            XWy << XSubSetTmp.t() * WSubSetTmp * ySubSetTmp;
            //                 //DEBUG
            //             printf("WsubsetTMP=\n");
            //             for(xId = 0; xId < numValid-2; xId++){
            //                 printf("%f ", WSubSetTmp.element( xId ) );
            //             }
            //             printf("XSubSetTMP=\n");
            //             for(yId = 0; yId < numValid-2; yId++){
            //                 for(xId = 0; xId < dimX+1; xId++){
            //                     printf("%f ", XSubSetTmp.element( yId , xId ) );
            //                 }
            //                 printf("\n");
            //             }
            //             printf("ySubSetTMP=\n");
            //             for(yId = 0; yId < numValid-2; yId++){
            //                 printf("%f ", ySubSetTmp.element( yId ) );
            //             }
        }
    }
    
    //     //DEBUG
    //     printf("Wsubset=\n");
    //     for(xId = 0; xId < numValid; xId++){
    //         printf("%f ", WSubSet.element( xId ) );
    //     }
    //     printf("XSubSet=\n");
    //     for(yId = 0; yId < numValid; yId++){
    //         for(xId = 0; xId < dimX+1; xId++){
    //             printf("%f ", XSubSet.element( yId , xId ) );
    //         }
    //         printf("\n");
    //     }
    //     printf("ySubSet=\n");
    //     for(yId = 0; yId < numValid; yId++){
    //         printf("%f ", ySubSet.element( yId ) );
    //     }
    
    
    //     return yPred;
    
    
    if(0){
        // Cholesky decomposition of SSQ
        try{
            LowerTriangularMatrix L = Cholesky(XWX);
            beta = L.t().i() * (L.i() * XWy);
        }
        catch(...){
            printf("WARN: Cholesky decomposition failed\n");
            return 0.0;
        }
    }else{
        
        try{
            LinearEquationSolver CLU = XWX; //decomposes XWX
            beta = CLU.i() * XWy;
            
            //Eigenvalues decomposition
            //             Matrix V; DiagonalMatrix D;
            //             EigenValues(XWX,D,V); // XWX = V * D * V.t(), XWX.i() = V * D.i() * V.t();
            //             beta = V * D.i() * V.t() * XWy;
            
            //SVD decomposition
            //             Matrix U, V; DiagonalMatrix D;
            //             SVD(XWX,D,U,V);                              // XWX = U * D * V.t()
            //             beta = U * D.i() * V.t() * XWy;
        }
        catch(...){
            printf("WARN: Matrix decomposition / inverse failed\n");
            return 0.0;
        }
    }
    
    //     //         DEBUG
    //     printf("beta=\n");
    //     for(xId = 0; xId < dimX+1; xId++){
    //         printf("%f ", beta.element( xId ) );
    //     }
    //     return yPred;
    
    
    yPred = beta.element( 0 );
    outBeta[0] = yPred;
    
    for(xId = 0; xId < dimX; xId++){
        yPred += xTest.element( xId ) * beta.element( xId + 1);
        outBeta[xId + 1] = beta.element( xId + 1);
    }
    
    return yPred;
}





