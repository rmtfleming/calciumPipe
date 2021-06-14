function param = initial_parameter_SPAMS()
%  DESCRIPTION:
%   initial parameters used on Dictionary Learning
%
%  USAGE:
%    param = initial_parameter_SPAMS()
%
%  ARGUMENTS:
% 
%  OUTPUT:
%    param       : current parameters
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

% param = [];
% param.mode        = 2; 
% param.lambda      = 0.2;
% param.numThreads  = 1; % number of threads
% param.batchsize   = 256;
% param.clean       = 1;
% param.modeD       = 0;
% param.posD        = 0;
% param.posAlpha    = 1;
% param.gamma1      = 0; 
% param.gamma2      = 0;
% param.lambda2     = 0;
% param.iter        = 1;
% param.pos         = 1; 

param = [];
param.mode        = 2; 
param.lambda      = 0.2;
param.iter        = 1;

