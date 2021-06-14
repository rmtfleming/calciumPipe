function [currentNeurons] = extractSingleCells(Coeffm,ReconsDict,currentNeurons,param)
% DESCRIPTION:
%    Extracts single and unique neurons from currentNeurons and fuse overlapping cells 
%    on the spatial domain.
%
% USAGE:   [cNeurons] = extractSingleCells(Coeffm,ReconsDict[,currentNeurons,param]);
%             currentNeurons, param optional
%
%
% ARGUMENTS:
%         data     double m (npixels) x n (nframes) matrix   (input signals)
%         nBasis   number of spatial/image patterns to factorize data    
%         param: struct (Important parameters to know, extracted from SPAMS library (mexTrainDL)    
%               param.type          = 'wavelet_denoising';  % Choose the extraction algorithm ('wavelet_denoising' (based on Reichnnek et. al. 2012), 'image_smoothing' or 'image_smoothing_powerlaw' )
%               param.Tsimilarity   = 'jointDistance';      % Choose a distance between gaussian 'jointDistance', 'MahalanobisDistance' or 'crosscorrelation'
%               param.actThr         = 10;                  % Only considers the image patterns that have been activated more the actThr frames 
%               param.permutation    = 1;                   % '1' choose a random permutation of ReconsDict, '0' not
%               param.lastWavelet    = 5;                   %  the last wavelet level from ReconsDict(:,:,i) max = 5 min = 1 (only wavelet_denoising)
%               param.prevWavelet    = 4;                   %  the previous last wavelet level from ReconsDict(:,:,i) max = 5 min = 1 (only wavelet_denoising)
%               param.ThrDenoising   = 4;                   % considers only pixels bigger than TheDenoising times the standard deviation of the noise (wavelet_denoising and image_smoothing_powerlaw)
%               param.AreaThr        = 80;                  % consider segments biggers than AreaThr 
%               param.MajorALength   = 25;                  % length of the major axis
%               param.percentageInt  = 0.9;                 % Intensity percentil
%               param.distCentroids  = 16;                  % distance between local maximas
%               param.SimThr         = 0.001;               % similarity measure to fuse two segments
%               param.lambda         = 0.05;                % used for inferring the activation patterns
%               param.sizeFilter     = [ 10 10];            % size of the filter (image_smoothing and image_smoothing_powerlaw)
%               param.sigmaS         = 3;                   % sigma of the gaussian filter
%               param.thrSmooth      = 0.3;                 % considers only pixels bigger than thrSmooth times the maximum intensity of the smoothed image
%               param.flag_watershed = 1;                   % image segmentation based on nearest points using local maximas ('0') or using watershed ('1')
%               param.power          = 30;                  % parameter used for the powerlaw
%
% OUTPUT: 
%         currentNeurons : struct
%            currentNeurons.imMask     % binary mask 
%            currentNeurons.imOrig     % binary mask times original image
%            currentNeurons.Area       % Area of the neuron
%            currentNeurons.maxF       % Maximum Normaliyed intensity
%            currentNeurons.X          % [x,y] coordinates of the mask
%            currentNeurons.obj.mu     % centroid of the cell
%            currentNeurons.obj.sigma  % covariance matrix of the cell (assuing a 2D-Gaussian)  
%            currentNeurons.Centroid   % centroid of the cell
%            currentNeurons.time       % activation pattern from Coeffm
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

if nargin <3,
    param = initial_parameter_Segmentation; param = param.all;
    currentNeurons = [];
elseif nargin ==3,
    param = initial_parameter_Segmentation; param = param.all;
elseif nargin == 4 && ~isstruct(param),
    param_init = initial_parameter_Segmentation; param_init = param_init.all;
    param = set_parameters_Segmentation(param_init,param);
elseif ~isstruct(param)
    error('Incorrect number of input parameters')
end

nBasis = size(ReconsDict,3);


if param.permutation 
    permutBasis = randperm(nBasis);
else
    permutBasis = 1:nBasis;
end

numActCoef = sum(Coeffm>0);
seg_this = extract_cells(param.type,param);

for i = permutBasis,
    if numActCoef(i) > param.actThr 
        im0 = ReconsDict(:,:,i);
        [featNeurons] = inferSingleCells(seg_this,im0,Coeffm(:,i));
        if ~isempty(currentNeurons) && ~isempty(featNeurons)
            currentNeurons = identifyDistinctCells(currentNeurons, featNeurons,param);
        elseif ~isempty(featNeurons)
            currentNeurons = featNeurons;
        end
    end
end

