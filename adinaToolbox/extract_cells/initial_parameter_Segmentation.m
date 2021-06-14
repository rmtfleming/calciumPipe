function [param_out] = initial_parameter_Segmentation()
%  DESCRIPTION:
%   initial parameters of the segmentation
%
%  USAGE:
%
%    param = initial_parameter_Segmentation()
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
% param                = [];
% param.actThr         = 5; 
% param.permutation    = 1;
% param.lastWavelet    = 4;
% param.prevWavelet    = 3;
% param.ThrDenoising   = 4;
% param.AreaThr        = 30;
% param.maxAreaThr     = 250; 
% param.Eccentricity    = 0.9;
% param.MajorALength   = 25;
% param.percentageInt  = 0;
% param.distCentroids  = 6;
% param.SimThr         = 0.005;
% param.lambda         = 0;
% param.sizeFilter     = [ 10 10];
% param.sigmaS         = 3;
% param.thrSmooth      = 0.3;
% param.resetNeurons   = 1;
% param.flag_watershed = 1;
% param.power          = 30;
% param_out.all = param;
% param_out.List = fieldnames(param);
% param_out.all.type          = 'wavelet_denoising';
% param_out.all.Tsimilarity   = 'crosscorrelation';


param                = [];
param.actThr         = 10; 
param.permutation    = 1;
param.AreaThr        = 200;
param.percentageInt  = 0;
param.distCentroids  = 6;
param.SimThr         = 0.005;
param.lambda         = 0;
param.sizeFilter     = [ 10 10];
param.sigmaS         = 3;
param.thrSmooth      = 0.18;
param.flag_watershed = 1;
param.power          = 30;
% param.resetNeurons   = 1;
param.type           = 'image_smoothing';
param.Tsimilarity    = 'crosscorrelation';
param_out.all = param;
