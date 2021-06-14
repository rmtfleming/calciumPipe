function [Coeffm,Dicts,ReconsDict,ReconsSeq] = infer_spams(data,nBasis,param)
% 
%
% DESCRIPTION:
%   infer_spams is a interface that calls mexTrainDL ( an efficient implementation of the
%   dictionary learning) technique presented in:
%
%     "Online Learning for Matrix Factorization and Sparse Coding"
%     by Julien Mairal, Francis Bach, Jean Ponce and Guillermo Sapiro
%     arXiv:0908.0050
%     
%     "Online Dictionary Learning for Sparse Coding"      
%     by Julien Mairal, Francis Bach, Jean Ponce and Guillermo Sapiro
%     ICML 2009.
%
%     The sparsity patterns addressed on the GUI needs the following
%     parameters:  
%
%       + temporal_sparsity.m
%           param.mode    = 0; 
%           param.modeD   = 0;
%           param.lambda  = x;  % where x>0 enforces sparsity on the temporal coefficients for each frame
%
%       + spatial_sparsity.m
%           param.mode    = 0; 
%           param.modeD   = 1;
%           param.lambda  = 0;
%           params.gamma  = x;  %where x>0 enforces that the dictionary elements has few nonzero pixels
%
%       + spatiotemporal_sparsity.m
%           param.mode    = 0; 
%           param.modeD   = 1;
%           param.lambda  = x1;  % where x1>0 enforces sparsity on the temporal coefficients for each frame
%           params.gamma  = x2;  % where x2>0 enforces that the dictionary elements has few nonzero pixels
%
%
%     Meaning of choosing param.mode and param.modeD:
%
%        1) if param.mode=0
%     min_{D in C} (1/n) sum_{i=1}^n (1/2)||x_i-Dalpha_i||_2^2  s.t. ...
%                                                  ||alpha_i||_1 <= lambda
%        2) if param.modeD=0
%           C={  D in Real^{m x p}  s.t.  forall j,  ||d_j||_2^2 <= 1 }
%        3) if param.modeD=1
%           C={  D in Real^{m x p}  s.t.  forall j,  ||d_j||_2^2 + ... 
%                                                  gamma1||d_j||_1 <= 1 }
% USAGE:
%   [Coeffm,Dicts,ReconsDict,ReconsSeq] = infer_spams(data,nBasis,param)
%
% ARGUMENTS:
%         data     double m (npixels) x n (nframes) matrix   (input signals)
%         nBasis   number of spatial/image patterns to factorize data    
%         param: struct (Important parameters to know, extracted from SPAMS library (mexTrainDL)           
%           param.lambda  (parameter)
%           param.iter (number of iterations).  If a negative number is 
%              provided it will perform the computation during the
%              corresponding number of seconds. For instance param.iter=-5
%              learns the dictionary during 5 seconds.
%            param.mode = 0;
%            param.posAlpha (optional, adds positivity constraints on the
%              coefficients, false by default)
%            param.modeD ( 0 or 1)
%            param.posD (optional, adds positivity constraints on the 
%              dictionary, false by default)
%            param.gamma1 (optional parameter for param.modeD == 1)
%            param.batchsize (optional, size of the minibatch, by default 
%              512)
%            param.clean (optional, true by default. prunes 
%              automatically the dictionary from unused elements).
%            param.numThreads (optional, number of threads for exploiting
%              multi-core / multi-cpus. By default, it takes the value -1,
%              which automatically selects all the available CPUs/cores).
%
%  OUTPUT: 
%         Dicts      : double m (npixels) x p (nBasis) matrix   (vectorized dictionary)
%         Coeffm     : double n (nframes) x p (nBasis) matrix   (temporal activation)
%         ReconsDict : double nrows x ncols x nBasis            (dictionary)
%         ReconsSeq  : double nrows x ncols x nFrames           (reconstructed sequence)
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".


if nargin ~=3
    error('Incorrect number of input parameters')
end

[nrows,ncols,nFrames] = size(data);
param.K = nBasis;

data       =  reshape(data,[nrows*ncols nFrames]);
Dicts      =  mexTrainDL(data,param);
Coeffm     =  full(mexLasso(data,Dicts,param))';

if nargout >= 3,
    for i = 1:nBasis,
        ReconsDict(:,:,i) = reshape(Dicts(:,i),[nrows ncols]);
    end
    ReconsDict = single(ReconsDict);
    if nargout >= 4,
        ReconsSeq = single(reshape(Dicts*Coeffm',[nrows ncols nFrames]));
    end
end