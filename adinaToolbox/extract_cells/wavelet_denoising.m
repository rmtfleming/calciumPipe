function [featNeurons,segmentation] = wavelet_denoising(dict,time,param)
% DESCRIPTION:
%    Extracts single and unique neurons from base function
%
% USAGE:  [featNeurons,segmentation] = wavelet_denoising(dict,time,param)

% ARGUMENTS:
%         data:      current basis function  
%         time:      temporal activation of the current bsais function
%         param: struct (Important parameters to know, extracted from SPAMS library (mexTrainDL)    
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
%         featNeurons : struct
%            featNeurons.imMask     % binary mask 
%            featNeurons.imOrig     % binary mask times original image
%            featNeurons.Area       % Area of the neuron
%            featNeurons.maxF       % Maximum Normaliyed intensity
%            featNeurons.X          % [x,y] coordinates of the mask
%            featNeurons.obj.mu     % centroid of the cell
%            featNeurons.obj.sigma  % covariance matrix of the cell (assuing a 2D-Gaussian)  
%            featNeurons.Centroid   % centroid of the cell
%            featNeurons.time       % activation pattern from Coeffm
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

std_noise = noiseLevelExtractionTruncated(dict);
[W] = waveletDenoisingDraguhn(dict);
nLevels = size(W,3);

%W = repmat(dict,[1 1 6]);
[im0] = W(:,:,nLevels-param.lastWavelet+1);
[im1] = W(:,:,nLevels-param.prevWavelet+1);
thr = param.ThrDenoising;

imMask =  im0>thr*std_noise;

im0Mask = dict.*imMask;
[imLabel0,numLabel] = bwlabel(im0Mask);

[nrows,ncols] = size(dict);
[X,Y] = meshgrid(1:ncols,1:nrows);

S_label_props = regionprops(imLabel0,'MajorAxisLength','EquivDiameter');
for i = 1:numLabel,
    if S_label_props(i).MajorAxisLength > param.MajorALength,
        % Use im1 (W_4) result
        imMask_1 =  (im1>thr(1)*std_noise).*(imLabel0 == i);
        im0Mask_1 = dict.*imMask_1;
        
        imMask(imLabel0==i) = imMask_1(imLabel0==i);
        im0Mask(imLabel0==i) = im0Mask_1(imLabel0 == i);
    end
end

[imLabel0,numLabel] = bwlabel(imMask);


if sum(im0Mask(:)~=0)
    %% Extracting local maximas
    imLocalMax = imregionalmax(im0Mask);
    [imLabelMax,numMax] = bwlabel(imLocalMax);


    S_props = regionprops(imLabelMax>0,'Centroid');
    Label0tolocalMax = zeros(numLabel,1);
    t = 1;
    for i = 1:numMax,
        row = round(S_props(i).Centroid(2));
        col = round(S_props(i).Centroid(1));
        val =  imLabel0(row,col);

        if val ~=0,
            localMaxtoLabel0(t) = val;
            Label0tolocalMax(val) = Label0tolocalMax(val)  +1;
            t = t+1;
        else
            
        end
    end
    numMax = t-1; % length(localMaxtoLabel0);

    if sum(Label0tolocalMax==0)<0
        idx = find(Label0tolocalMax==0);
        Label0tolocalMax(idx) = [];
        %there is a region without local maxima, it can not happen
        keyboard()
    end


    %% Deleting local maxima that does not exceed the 90th percentile of the pixel intensities

    listDelMax = [];
    for i = 1:numMax,
        intPixRegion = im0Mask(imLabel0 == localMaxtoLabel0(i));
        numPixels = length(intPixRegion);
        row = round(S_props(i).Centroid(2));
        col = round(S_props(i).Centroid(1));
        intLocalMax =  im0Mask(row,col);
        S_props(i).localMax = intLocalMax;
        perPixInt = length(find(intPixRegion<intLocalMax))/numPixels;
        S_props(i).perPixInt = perPixInt;
        if perPixInt < param.percentageInt
            %delete Local Max
            listDelMax = [listDelMax i];
        end
    end



    if ~isempty(listDelMax)
        S_props(listDelMax) = [];
        for i = 1:length(listDelMax),
            Label0tolocalMax(localMaxtoLabel0(listDelMax(i))) = Label0tolocalMax(localMaxtoLabel0(listDelMax(i))) -1;
        end
        localMaxtoLabel0(listDelMax) = [];
    end

    %% Deleting close local maxima below 16 pixels
    numMax = length(S_props);
    distLocal = Inf*ones(numMax);
    for i = 1:numMax,
        for j = i+1:numMax,
            if i == j 
                distLocal(i,i) = Inf;
            else
                distLocal(i,j) = sqrt(sum((S_props(i).Centroid-S_props(j).Centroid).^2));
            end
        end
    end

    [r,c,~] = find(distLocal<param.distCentroids);
    listDelMax = [];
    if ~isempty(r),
        for i = 1:length(r),
               if S_props(r(i)).localMax >=  S_props(c(i)).localMax
                   listDelMax = [listDelMax c(i)];
               else
                   listDelMax = [listDelMax r(i)];
               end   
        end
    end

    if ~isempty(listDelMax)
        S_props(listDelMax) = [];
        for i = 1:length(listDelMax),
            try
                Label0tolocalMax(localMaxtoLabel0(listDelMax(i))) = Label0tolocalMax(localMaxtoLabel0(listDelMax(i))) -1;
            catch
                keyboard()
            end
        end
        localMaxtoLabel0(listDelMax) = [];
    end
    

    %% Deleting regions without local maxima
    
    if sum(Label0tolocalMax==0)>0
        idx = find(Label0tolocalMax==0);
        for i = 1:length(idx),
            im0Mask(imLabel0 == idx(i)) = 0;
            imMask(imLabel0 == idx(i)) = 0;
        end

    end

    %% Cluster the points belonging to the extracted centroids
    
    if ~param.flag_watershed
        %% ClasN segmentation
        mu = zeros(length(S_props),2);
        for i = 1:length(S_props),
            mu(i,:) = S_props(i).Centroid;
        end
        
        [I J] = find(imMask>0);
        distN = zeros(length(I),length(S_props));
        for i = 1:length(S_props),
            distN(:,i) = sqrt((I-mu(i,2)).^2+(J-mu(i,1)).^2);
        end
        if isempty(S_props)
            featNeurons = [];
            segmentation = [];
            return
        end
        
        [~,pos] = min(distN,[],2);
        
        %plot classify pixels;
        segmentation = zeros(size(imMask));
        for i = 1:length(I),
            segmentation(I(i),J(i)) = pos(i);
        end
        
        % extraction of features
        featNeurons = [];
        t = 1;
        for i = 1:length(S_props),
            imMask_tmp = (segmentation == i).*(imLabel0 == localMaxtoLabel0(i));
            S_tmp = regionprops(imMask_tmp ,'Centroid','Area','Eccentricity');
            if length(S_tmp)==1
                if S_tmp.Area> param.AreaThr,
                      [I,J] = find(imMask_tmp>0);
                    tmpImg = imMask_tmp.*dict;
                   
                    featNeurons(t).obj = gridGMM([J,I],tmpImg(imMask_tmp>0),1);
                    tmpProb = reshape(pdf(featNeurons(t).obj,[X(:) Y(:)]),[nrows ncols]);
                    imMask_tmp = imMask_tmp.*(tmpProb>0.001);
                           
                                        
                    featNeurons(t).imMask = imMask_tmp ;
                    featNeurons(t).im = (segmentation == i+1).*im0Mask.*imMask_tmp ;
                    featNeurons(t).imOrig = imMask_tmp.*dict;
                    featNeurons(t).Area = S_tmp.Area;
                    featNeurons(t).maxF = max(dict(:).*(segmentation(:) == i+1).*imMask_tmp(:));
                    [I,J] = find(featNeurons(t).imMask>0);
                    featNeurons(t).X = [J,I];
                    featNeurons(t).X_v =  tmpImg(imMask_tmp>0);
                    featNeurons(t).Centroid = featNeurons(t).obj.mu;
                    featNeurons(t).time = time;
                end
            end
        end
    

    else
        %% Watershed segmentation
        segmentation = watershedSegmentationWAPL(im0Mask,0);
        
        regionsWatershed = unique(segmentation);
        regionsWatershed(regionsWatershed==0) = []; %0's are the region edges
        regionsWatershed(regionsWatershed==1) = []; %1 is the background
        if isempty(regionsWatershed)
            featNeurons = [];
            segmentation = [];
            return
        end
        segmentation(segmentation==1) = 0;
        
        %%%
        for i = 1:length(regionsWatershed)
            tmp = (segmentation == regionsWatershed(i));
            corresp = zeros(1,length(S_props));
            for j = 1:length(S_props)
                corresp(j) = tmp(round(S_props(j).Centroid(2)), round(S_props(j).Centroid(1)));
            end
            if ~any(corresp)
                segmentation(segmentation==regionsWatershed(i)) = 0;
            end
        end
        
        regionsWatershed = unique(segmentation);
        regionsWatershed(regionsWatershed==0) = [];
        %%%
        
        featNeurons = [];
        t = 1;
        for i = 1:length(regionsWatershed),
            imMask_tmp = (segmentation == regionsWatershed(i));%.*(imLabel0 == localMaxtoLabel0(i));
            S_tmp = regionprops(imMask_tmp ,'Centroid','Area','Eccentricity');
            if length(S_tmp)==1
                if (S_tmp.Area> param.AreaThr) && (S_tmp.Area< param.maxAreaThr) && (S_tmp.Eccentricity<param.Eccentricity),
                    
                    [I,J] = find(imMask_tmp>0);
                    tmpImg = imMask_tmp.*dict;
                   
                    %featNeurons(t).obj = gridGMM([J,I],tmpImg(imMask_tmp>0),1);
                    featNeurons(t).obj = gmdistribution.fit([J,I],1);
                    %tmpProb = reshape(pdf(featNeurons(t).obj,[X(:) Y(:)]),[nrows ncols]);
                    imMask_tmp = imMask_tmp;%.*(tmpProb>0.0005);
                           
                                        
                    featNeurons(t).imMask = imMask_tmp ;
                    featNeurons(t).im = imMask_tmp.*dict; %(segmentation == i+1).*im0Mask.*imMask_tmp ;
                    featNeurons(t).imOrig = imMask_tmp.*dict;
                    featNeurons(t).Area = S_tmp.Area;
                    featNeurons(t).maxF = max(dict(:).*(segmentation(:) == i+1).*imMask_tmp(:));
                    [I,J] = find(featNeurons(t).imMask>0);
                    featNeurons(t).X = [J,I];
                    featNeurons(t).X_v =  tmpImg(imMask_tmp>0);
                    featNeurons(t).Centroid = featNeurons(t).obj.mu;
                    featNeurons(t).time = time;
                    t=t+1;
                else
                    segmentation(segmentation==regionsWatershed(i)) = 0;
                end
            end
        end
    
    end

else
    featNeurons = [];
    segmentation = [];
end
