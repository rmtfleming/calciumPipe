function [roi] = extractingROI(dirBase,fileExt,stepChannel)
%  
%  DESCRIPTION:
%    Main script for extracting neural activity. This process load all the
%    data needed into the RAM.
%
%  USAGE:
%
%    function [roi] = extractingROI(dirBase,fileExt,stepChannel)
%
%  ARGUMENTS:
%
%    dirBase     : full path of the base directory where images resides
%    fileExt     : wildcard for matching imaging files within dirBase
%    stepChannel : how many channels you have in your imaging data
%                  (functional/green channel is the first channel)
%  
%  OUTPUT:
%    roi:
%       roi(i).Img : cell appearance (weighted mask) of the i-th cell
%       roi(i).indices  : pixel location of the i-th cell
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".


%% Extracting the information of the sequence

files = dir([dirBase fileExt]);
frames = 0;
for i = 1:length(files),
    fname = files(i).name;
    info = imfinfo([dirBase fname]); frames = frames+length(info)/stepChannel;
end
params.dirBase  = dirBase;
params.files    = files;
params.nrows    = info(1).Height;
params.ncols    = info(1).Width;
params.nFrames  = frames;
params.step     = stepChannel;

%% Main parameters

params.BorderSize = 16; % how much you want to cut from the borders
params.szTemporalBlock = 10000; % how many frames max per block
params.halfWindowTime = 400; % size of window for median filter (half) in frames
params.nSpatialBlocks = 8; % blocks (note: total blocks is this number ^2 ; note that rem(frameSize,nSpatialBlocks) must =0)
params.SpatialConfidence = 20; % size of overlap ; should be size of cell in pixels

%% Variables for PreProcessing Step

params.prepro.method = 2; % method for f0 estimation ; SEE CODE 
params.prepro.halfWinMedian = params.halfWindowTime;
params.prepro.quantile = 0.20; % which quantile to use for f0
params.prepro.stepF0 = 100; % to speed up f0 estimation, f0 is only computed every this many frames
params.prepro.tempSigma = 5; % temporal smoothing for compute f0
params.prepro.spaSigma = 5; % spatial smoothing for f0
params.prepro.wavScale = [ 4 4 3]; % wavelet in spaceX,spaceY,time ; max is 5


%% Variables for Segmentation

params.cellSort.nIter = -2; % dictionary learning iterations (-2 computes the DL only for 2 seconds, while setting it to 200 computes DL for 200 iterations)
params.cellSort.lambda = 0.2;
params.seg.lastWavelet    = params.prepro.wavScale(1);
params.seg.prevWavelet    = params.prepro.wavScale(1)-1;
params.seg.ThrDenoising   = 3;
params.seg.AreaThr        = 30;
params.seg.maxAreaThr     = 250;
params.seg.MajorALength   = 25;
params.seg.SimThr         = 0.4;


%% Reading the files of the directory (Green/Red Channel)

params.idFiles = zeros(1,params.nFrames);
params.idSlice = zeros(1,params.nFrames);

t = 1;
for i = 1:length(files),
    fname = files(i).name;
    info = imfinfo([dirBase fname]); 
    for j = 1:stepChannel:length(info),
        params.idFiles(t) = i;
        params.idSlice(t) = j;
        t = t+1;
    end
end

%% Dividing in temporal blocks
if params.nFrames > params.szTemporalBlock,
    %split sequences
    disp('splitting')
    params.stepFrames = params.szTemporalBlock;
    params.nTemporalBlocks = ceil(params.nFrames/params.stepFrames);
    % 10 seconds/leg of the quantile filter
    for i = 1: params.nTemporalBlocks ,
        rSeq{i}  = max(1,(i-1)*params.stepFrames+1):min(params.nFrames,i*params.stepFrames);  %Actual Temporal Block
        rWSeq{i} = max(1,(i-1)*params.stepFrames+1-params.halfWindowTime):min(params.nFrames,i*params.stepFrames+params.halfWindowTime); % Temporal Block + confident window (paddarray)
        rMar{i}  = max(rSeq{i}(1)-rWSeq{i}(1))+(1:length(rSeq{i}));    % Position of the acutal temporal block in rWSeq
    end
else 
    %no split
    disp('nosplitting')
    params.nTemporalBlocks  = 1;
    rSeq{1}  = 1:params.nFrames;
    rWSeq{1} = 1:params.nFrames;
    rMar{1}  = 1:params.nFrames;
end

params.BorderSize = 16;

rows = 1:params.nrows;
cols = 1:params.ncols;
rows = rows(params.BorderSize+1:end-params.BorderSize);
cols = cols(params.BorderSize+1:end-params.BorderSize);

params.ncolsROI = length(cols);
params.nrowsROI = length(rows);
params.stepSpatial = params.ncolsROI/params.nSpatialBlocks;

for iSpatial = 1:params.nSpatialBlocks,
    rROI{iSpatial} = (iSpatial-1)*params.stepSpatial + (1:params.stepSpatial);
    cROI{iSpatial} = rROI{iSpatial};
    rWin{iSpatial} = max((iSpatial-1)*params.stepSpatial+1-params.SpatialConfidence,1):min(params.ncolsROI,iSpatial*params.stepSpatial+params.SpatialConfidence);
    cWin{iSpatial} = rWin{iSpatial};
    rMarS{iSpatial} = max(rROI{iSpatial}(1)-rWin{iSpatial}(1))+(1:params.stepSpatial);    
    cMarS{iSpatial} = rMarS{iSpatial};
end

paramsOut = [];
curBlock = 1;
for i = 1:params.nTemporalBlocks,
    for j = 1:params.nSpatialBlocks,
            for z = 1:params.nSpatialBlocks,
                paramsOut(curBlock).main          = params;
                paramsOut(curBlock).time.rSeq     = rSeq{i};
                paramsOut(curBlock).time.rWSeq    = rWSeq{i};
                paramsOut(curBlock).time.rMar     = rMar{i};
                paramsOut(curBlock).spatial.rROI  = rROI{j};
                paramsOut(curBlock).spatial.cROI  = cROI{z};
                paramsOut(curBlock).spatial.rWin  = rWin{j};
                paramsOut(curBlock).spatial.cWin  = cWin{z};
                paramsOut(curBlock).spatial.rMarS = rMarS{j};
                paramsOut(curBlock).spatial.cMarS = cMarS{z};
                curBlock = curBlock +1;
            end
    end
end

%% Processing Blocks

[InferCells] = automaticROIextraction(paramsOut);
[roi] = mappingBlockstoCells(InferCells);
