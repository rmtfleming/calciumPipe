function param = initial_parameter_NMF()
%  DESCRIPTION:
%   initial parameters used on NonNegative Matrix Factorization
%
%  USAGE:
%    param =  initial_parameter_NMF()
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
param = [];
param.iter        = 1000;
param.dis         = 0;
param.residual    = 1e-4;
param.tof         = 1e-4;
param.posD        = 1;

param.distance    = 'kl'; %'kl', 'ls'

param.beta        = 0.1;

param.orthogonalA = 0;
param.orthogonalB = 0;
param.orthogonal  = [param.orthogonalA,param.orthogonalB];




