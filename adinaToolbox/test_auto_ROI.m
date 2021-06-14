clear all, clc
close all,

dirBase = '/mnt/siham_hachi/data/fast-oopsi/2lane/stacks/';
fileExt = 'full_small_stack.tif';
stepChannel=1;
matFile = '/mnt/siham_hachi/data/fast-oopsi/2lane/mat/full_small_stack.mat';

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
params.seg.SimThr         = 0.005;


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
params.stepSpatial = params.stepSpatial - mod(params.stepSpatial,1);

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

%%

load(matFile);

for i=1:Im.T
    img=reshape(Im.DataMat(:,i),Im.h,Im.w);
    data(:,:,i)=img(:,:);
end

nBlocks = length(paramsOut);
nrows   = paramsOut(1).main.nrows;
ncols   = paramsOut(1).main.ncols;
nFrames = paramsOut(1).main.nFrames;

for id = 1:nBlocks,
    data_block = double(data(paramsOut(id).spatial.rWin,paramsOut(id).spatial.cWin,paramsOut(id).time.rWSeq));
    disp('Entering Data Normalization')
    seq_param_tmp = extractingNormalizedSeq(data_block,paramsOut(id).main.prepro);
    neurons = processingCellDetectionBlock(seq_param_tmp.data_denoise(:,:,paramsOut(id).time.rMar),round(length(paramsOut(id).time.rSeq)*0.05),paramsOut(id).spatial.rMarS,paramsOut(id).spatial.cMarS,paramsOut(id));
    InferCells(id).neurons = neurons;
    InferCells(id).params = paramsOut(id);
end

%%
% %% infer spams
% 
%     if strcmp(meth,'spams')
% 
%         param_spams.mode        = 2;
%         param_spams.lambda      = 0.2;
% %         param_spams.numThreads  = 1;                                  % number of threads
% %         param_spams.batchsize   = 256;
%     %     param_spams.clean       = 1;
% %         param_spams.modeD       = 0;
%     %     param_spams.posD        = 0;
%     %     param_spams.posAlpha    = 1;
% %         param_spams.gamma1      = 0;
%     %     param_spams.gamma2      = 0;
%     %     param_spams.lambda2     = 0;
%         param_spams.iter        = 1;
%     %     param_spams.pos         = 1;
%         [Coeff,Dict,ReconsDict] = infer_spams(X,1,param_spams);
% 
%     else if strcmp(meth,'nmf')
%             param_nmf = [];
%             param_nmf.iter        = 1000;
%             param_nmf.dis         = 0;
%             param_nmf.residual    = 1e-4;
%             param_nmf.tof         = 1e-4;
%             param_nmf.posD        = 1;
%             param_nmf.distance    = 'kl'; %'kl', 'ls'
%             param_nmf.beta        = 0.1;
%             param_nmf.orthogonalA = 0;
%             param_nmf.orthogonalB = 0;
%             param_nmf.submethodType = 'nmfnnls';
%             param_nmf.orthogonal  = [param_nmf.orthogonalA,param_nmf.orthogonalB];
%             [Coeffm,Dicts,ReconsDict,ReconsSeq] = infer_nmf(data,1,param_nmf);
%         end
%     end
% 
% 
%     %% exract single cells
% 
% %     param.Tsimilarity   = 'jointDistance';                  % Choose a distance between gaussian 'jointDistance', 'MahalanobisDistance' or 'crosscorrelation'
%     param.actThr         = 5;                               % Only considers the image patterns that have been activated more the actThr frames
%     param.permutation    = 1;                               % '1' choose a random permutation of ReconsDict, '0' not
%     % param.lastWavelet    = 5;                             %  the last wavelet level from ReconsDict(:,:,i) max = 5 min = 1 (only wavelet_denoising)
%     % param.prevWavelet    = 4;                             %  the previous last wavelet level from ReconsDict(:,:,i) max = 5 min = 1 (only wavelet_denoising)
%     % param.ThrDenoising   = 4;                             % considers only pixels bigger than TheDenoising times the standard deviation of the noise (wavelet_denoising and image_smoothing_powerlaw)
%     param.AreaThr        = 30;                              % consider segments biggers than AreaThr
%     % param.MajorALength   = 5;                             % length of the major axis
%     param.percentageInt  = 0;                               % Intensity percentil
%     param.distCentroids  = 1;                              % distance between local maximas
%     param.SimThr         = 0.005;                           % similarity measure to fuse two segments
%     param.lambda         = 0;%0.05;                         % used for inferring the activation patterns
%     param.sizeFilter     = [10 10];                         % size of the filter (image_smoothing and image_smoothing_powerlaw)
%     param.sigmaS         = 3;                               % sigma of the gaussian filter
%     param.thrSmooth      = 0.2; %0.09                       % considers only pixels bigger than thrSmooth times the maximum intensity of the smoothed image
%     param.flag_watershed = 0;                               % image segmentation based on nearest points using local maximas ('0') or using watershed ('1')
% %     param.power          = 100;                           % parameter used for the powerlaw
%     param.type = 'image_smoothing';
% 
%     cNeurons = [];
%     if ~isempty(data)
%         if ~isempty(ReconsDict)
%             disp('Extracting single cells ... '),
%             if ~isempty(cNeurons)
%                  % ask if they want to update the current extracted cells                 
%                  [cNeurons] = extractSingleCells(Coeff,ReconsDict,cNeurons,param);
%             else                
%                  [cNeurons] = extractSingleCells(Coeff,ReconsDict,[],param);
%             end
% 
%             if ~isempty(cNeurons)                
%                 disp('Extracted Single Cells is done.')
%                 disp('Computing Temporal Evolution for the extracted single Cells')
%                 [Coeff_single,Dicts_single,ReconsDict_single] = inferCoeffSingleCells(X,cNeurons,param);
% 
%                 %% filter neurons that are not used for reconstructing
%                 [cNeurons,Coeff_single,Dicts_single,ReconsDict_single] = proneNeurons(cNeurons,Coeff_single,Dicts_single,ReconsDict_single);
% 
%                 disp('Temporal Evolution is computed')    
%             else
%                 disp('No found single cells.')
%             end
% 
%         else
%              msgbox('Please run Matrix Factorization algorithm.','Error','error','modal');
%         end
%     else
%         msgbox('Please load a sequence before setting the number of basis functions.','Error','error','modal');
%     end
% 
% %     if ~isfield(Im,'cNeurons')
% %         Im.cNeurons = cNeurons;
% %     end
% %     cNeurons_all{i}=cNeurons;
% % end
% [idFig] = displaySingleCells(data,cNeurons);
% % displayMultipleDetectedSignals(data,Coeff_single,cNeurons,8,2);
