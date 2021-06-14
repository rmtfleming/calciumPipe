#include "mex.h"
#include <cmath>
#include <stdio.h>
#include <string.h>

// mex mexUpdateV.cpp -l blas -largeArrayDims

typedef ptrdiff_t intt;

// Argument 0: V, whose dim is p x r
//          1: X, whose dim is n x p
//          2: U, whose dim is n x r
//          3: inv_Zeta, whose dim is p x r
//          4: lambda
//          5: params

//-------------------------------------------------------------------------

void parse_params(const mxArray *params, intt *max_it, bool *pos, double *normalization);

//=========================================================================

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{

   double *B, *A, *IX, *B_output, *IX_output;
   int nrows = mxGetM(prhs[0]), ncols =  mxGetN(prhs[0]);

   B        = (double*)mxGetPr(prhs[0]);
   IX       = (double*)mxGetPr(prhs[1]);    
   A        = (double*)mxGetPr(prhs[2]);
   //-----------------------------------------------------------------------
   //V_output that is initialized as a copy of V
   
   plhs[0]  = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
   plhs[1]  = mxCreateDoubleMatrix(nrows, ncols, mxREAL);
   
   B_output = (double*)mxGetPr( plhs[0] );
   IX_output = (double*)mxGetPr( plhs[1] );
   
   
   
    for(int i = 0;i<ncols;i++){
        int t = 0; 
        int pos_value = 0;
        bool flag = 1;
        for(int j = 0; (j+pos_value)<nrows;j++){
           if (IX[i*nrows+j+pos_value]!=1){

                if (B[i*nrows+j] < A[i]){
                    B_output[i*nrows+t] = B[i*nrows+(j+pos_value)];
                    IX_output[i*nrows+t] = IX[i*nrows+(j+pos_value)]-1; 
                }else{ 
                    if (B[i*nrows+j] == A[i]){
                        if(flag==1){
                            pos_value = -1;
                            B_output[i*nrows+t] = A[i]; 
                            IX_output[i*nrows+t] = nrows;  
                            flag = 0;
                        }else{
                            B_output[i*nrows+t] = B[i*nrows+(j+pos_value)];
                            IX_output[i*nrows+t] = IX[i*nrows+(j+pos_value)]-1;                             
                        }
                    }else{
                        if(flag ==1){
                            pos_value = -1;
                            B_output[i*nrows+t] = A[i]; 
                            IX_output[i*nrows+t] = nrows;
                            flag = 0;
                        }else{
                            B_output[i*nrows+t] = B[i*nrows+(j+pos_value)];
                            IX_output[i*nrows+t] = IX[i*nrows+(j+pos_value)]-1;
                        }
                    } 
                }
                t++;    
            }else{
                pos_value == 0;
            }           
        } 
        if(flag==1){            
            B_output[i*nrows+t] = A[i]; 
            IX_output[i*nrows+t] = nrows;
        }
 
    }
}
