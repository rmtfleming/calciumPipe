function [neurons] = processingCellDetectionBlock(data_normalized2,nBasis,rMar,cMar,params)
% DESCRIPTION:
%    Extracts single and unique neurons from the block
%
% USAGE:   
%       [neurons] = processingCellDetectionBlock(data_normalized2,nBasis,rMar,cMar,params)  
%
%
% ARGUMENTS:
%         data_normalized2     double m (npixels) x n (nframes) matrix   (input signals)
%         nBasis   number of spatial/image patterns to factorize data    
%         rMar,cMar : specifies the current region of interest of the block
%                     because each block contains a margin to avoid border effects.
%         params: struct (Important parameters to know, extracted from SPAMS library (mexTrainDL)    
%               params.type          = 'wavelet_denoising';  % Choose the extraction algorithm ('wavelet_denoising' (based on Reichnnek et. al. 2012), 'image_smoothing' or 'image_smoothing_powerlaw' )
%               params.Tsimilarity   = 'jointDistance';      % Choose a distance between gaussian 'jointDistance', 'MahalanobisDistance' or 'crosscorrelation'
%               params.actThr         = 10;                  % Only considers the image patterns that have been activated more the actThr frames 
%               params.permutation    = 1;                   % '1' choose a random permutation of ReconsDict, '0' not
%               params.lastWavelet    = 5;                   %  the last wavelet level from ReconsDict(:,:,i) max = 5 min = 1 (only wavelet_denoising)
%               params.prevWavelet    = 4;                   %  the previous last wavelet level from ReconsDict(:,:,i) max = 5 min = 1 (only wavelet_denoising)
%               params.ThrDenoising   = 4;                   % considers only pixels bigger than TheDenoising times the standard deviation of the noise (wavelet_denoising and image_smoothing_powerlaw)
%               params.AreaThr        = 80;                  % consider segments biggers than AreaThr 
%               params.MajorALength   = 25;                  % length of the major axis
%               params.percentageInt  = 0.9;                 % Intensity percentil
%               params.distCentroids  = 16;                  % distance between local maximas
%               params.SimThr         = 0.001;               % similarity measure to fuse two segments
%               params.lambda         = 0.05;                % used for inferring the activation patterns
%               params.sizeFilter     = [ 10 10];            % size of the filter (image_smoothing and image_smoothing_powerlaw)
%               params.sigmaS         = 3;                   % sigma of the gaussian filter
%               params.thrSmooth      = 0.3;                 % considers only pixels bigger than thrSmooth times the maximum intensity of the smoothed image
%               params.flag_watershed = 1;                   % image segmentation based on nearest points using local maximas ('0') or using watershed ('1')
%               params.power          = 30;                  % parameter used for the powerlaw
%
% OUTPUT: 
%         neurons : struct
%            neurons.imMask     % binary mask 
%            nNeurons.imOrig     % binary mask times original image
%            neurons.Area       % Area of the neuron
%            neurons.maxF       % Maximum Normaliyed intensity
%            neurons.X          % [x,y] coordinates of the mask
%            neurons.obj.mu     % centroid of the cell
%            neurons.obj.sigma  % covariance matrix of the cell (assuing a 2D-Gaussian)  
%            neurons.Centroid   % centroid of the cell
%            neurons.time       % activation pattern from Coeffm
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".
disp('Arraging data for Dictionary Learning');
X = preprocessingData(single(data_normalized2),'non','l2-norm');
[Coeff,~,ReconsDict] = inferGroupedCells(X,nBasis,'SPAMS',{'iter',params.main.cellSort.nIter,'lambda',params.main.cellSort.lambda});
ReconsDict(ReconsDict<0) = 0;
disp('Dictionary Learning Done')
disp('Detecting Neurons Block')

[tmp] = extractSingleCells(double(Coeff),double(ReconsDict),[],{'ThrDenoising',params.main.seg.ThrDenoising ,'SimThr',params.main.seg.SimThr,'lastWavelet',params.main.seg.lastWavelet,'prevWavelet',params.main.seg.prevWavelet ,'AreaThr',params.main.seg.AreaThr,'maxAreaThr',params.main.seg.maxAreaThr ,'MajorALength',params.main.seg.MajorALength});  
t = 1;
disp('Neurons extracted');

neurons = cell(-1);

for z  =1:length(tmp),
    if ((sum(sum((tmp(z).imOrig(rMar,cMar))>0)))/(sum(sum((tmp(z).imOrig)>0)))>0.50),
        neurons{t} = tmp(z);
        t = t+1;
    end
end
disp(['Number of neurons found ' num2str(length(tmp)) ' but only ' num2str(t-1) ' found inside the ROI'])


