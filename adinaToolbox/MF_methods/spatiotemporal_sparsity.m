function [params] = spatiotemporal_sparsity()
%  DESCRIPTION:
%   initial parameters of Dictionary Learning using only spatio-temporal sparsity
%
%  USAGE:
%
%    param = spatial_sparsity()
%
%  ARGUMENTS:
% 
%  OUTPUT:
%    param       : parameteres needed for the image enhancement
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

if nargin ~= 0
    error('Incorrect number of parameters')
end

params.all           = initial_parameter_SPAMS;
params.Type          = 'SPAMS';
params.subType       = 'spatial_sparsity';
params.all.modeD     = 1;
params.subList       = {'lambda','gamma1','batchsize','numThreads','iter','pos','posD'};