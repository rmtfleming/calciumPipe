function param = initial_parameter_preprocessing()
%  DESCRIPTION:
%   initial parameters of the preprocessing step. The parameters can be
%   modified here, on the main function or set_parameters
%
%  USAGE:
%
%    param = initial_parameter_preprocessing()
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
param = struct();

%% Variables for PreProcessing Step

param.method = 2; % method for f0 estimation ; SEE CODE 
param.halfWinMedian = 10;
param.quantile = 0.20; % which quantile to use for f0
param.stepF0 = 5; % to speed up f0 estimation, f0 is only computed every this many frames
param.tempSigma = 5; % temporal smoothing for compute f0
param.spaSigma = 5; % spatial smoothing for f0
param.wavScale = [5 5 3]; % wavelet in spaceX,spaceY,time ; max is 5
param.wavScaleSpa = 4;
param.wavScaleTime = 3;
param.BorderSize = 16; 


