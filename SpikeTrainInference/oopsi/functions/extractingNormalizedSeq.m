function [seq_param] = extractingNormalizedSeq(data,params)
%  DESCRIPTION:
%    Computes different image enhancement methods to extract reliable information as (F-F_0)/F_0  and wavelet transform,
%    and also removes line artifacts from line-to-line motion compensation algorithm.
%
%  USAGE:
%
%    [seq_param] = extractingNormalizedSeq(data,params)
%
%  ARGUMENTS:
%
%    data            : raw data
%    params          : input parameters
%        .method = 1 : original raw data
%        .method = 2 : wavelet transform on original raw data
%        .method = 3 : (F-F_0)/F0 
%        .method = 0 : wavelet transform on (F-F_0)/F0 
%   
%    params.prepro.halfWinMedian :
%    params.prepro.quantile      : which quantile to use for f0 computation
%    params.prepro.stepF0        : to speed up f0 estimation, f0 is only computed every this many frames
%    params.prepro.tempSigma     : temporal smoothing for compute f0
%    params.prepro.spaSigma      : spatial smoothing for f0
%    params.prepro.wavScale      : wavelet in spaceX,spaceY,time ; max is 5
% 
%  OUTPUT:
%    seq_param   : enhanced image sequence.
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

if nargin <2,
    params = initial_parameter_preprocessing;
elseif nargin ==2 && ~isstruct(params),
    param_init = initial_parameter_preprocessing; 
    params = set_parameters_preprocessing(param_init,params);
elseif ~isstruct(params) 
    error('Incorrect input parameters')
end

disp('Starting Preprocessing. Normalizing Block before Cell Detection')
[nrows,ncols,~] = size(data);
bias = 100;
disp('Preprocessing: Removing Lines is processing.')
data_bis = removingLinesMotionArtifacts(data)+bias;
disp('Preprocessing: Removing Lines done.')
disp('Preprocessing: Normalizing Block before Cell Detection')
switch params.method
    case 3
        disp('Preprocessing: delta f')
        h_3(1,1,:) = fspecial('gaussian',[1 5*params.tempSigma],params.tempSigma);
        Wgausstime = imfilter(data_bis, h_3,'corr','symmetric','same');
        [Fbase] = fastExtractingQuantilesSeq(reshape(double(Wgausstime),[],size(data_bis,3)),params.halfWinMedian,params.quantile,params.stepF0);
        Fbase = reshape(Fbase(:,:,end)',[nrows ncols size(data_bis,3)]);
        Fbase_smooth = imfilter(Fbase,fspecial('gaussian',[5*params.spaSigma 5*params.spaSigma],params.spaSigma ),'corr','symmetric','same');
        W = (data_bis-Fbase_smooth)./Fbase_smooth;
    case 2
        disp('Preprocessing: wavelet denoising')
        W = waveletDenoisingDraguhn3DScale(single(data_bis),params.wavScale);
    case 4
        disp('Preprocessing: delta f + wavelet denoising')
        h_3(1,1,:) = fspecial('gaussian',[1 5*params.tempSigma],params.tempSigma);
        Wgausstime = imfilter(data_bis, h_3,'corr','symmetric','same');
        [Fbase] = fastExtractingQuantilesSeq(reshape(double(Wgausstime),[],size(data_bis,3)),params.halfWinMedian,params.quantile,params.stepF0);
        Fbase = reshape(Fbase(:,:,end)',[nrows ncols size(data_bis,3)]);
        Fbase_smooth = imfilter(Fbase,fspecial('gaussian',[5*params.spaSigma  5*params.spaSigma ],params.spaSigma ),'corr','symmetric','same');
        AF = (data_bis-Fbase_smooth)./Fbase_smooth;
        W = waveletDenoisingDraguhn3DScale(AF,params.wavScale);
    case 1
        disp('Preprocessing: none')
        W = data;
    otherwise
        error('extractingNormalizedSeq: not implemented');
end
disp('Preprocessing: Normalizing done.')
seq_param.data_denoise = double(W);
seq_param.data_denoise_pC = seq_param.data_denoise;
