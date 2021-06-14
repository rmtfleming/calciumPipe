clear all, clc
close all,

% following the GUI
addpath('/home/bpf00038/Desktop/Work/spams-matlab/build');

% load('/mnt/siham_hachi/data/leuven_data/Siham_13052015/z_time26/data_maxp.mat');
data_maxp = data2D;
data1 = data_maxp;
data = data1;

% choose method 'spams' or 'nmf'
meth = 'spams';

%% setting parameters

% Parameters for PreProcessing Step
    param_preproc.method            = 2;                            % method for f0 estimation ; SEE CODE 
    param_preproc.halfWinMedian     = 400;
    param_preproc.quantile          = 0.20;                         % which quantile to use for f0
    param_preproc.stepF0            = 5;                            % to speed up f0 estimation, f0 is only computed every this many frames
    param_preproc.tempSigma         = 5;                            % temporal smoothing for compute f0
    param_preproc.spaSigma          = 5;                            % spatial smoothing for f0
    param_preproc.wavScale          = [4 4 3];                      % wavelet in spaceX,spaceY,time ; max is 5
    % param_preproc.wavScaleSpa       = 4;
    % param_preproc.wavScaleTime      = 3;
    % param_preproc.BorderSize        = 16;

% runEnhancement
    % data = preprocessingData(data,'mean','squared l2-norm');
    data_tmp = extractingNormalizedSeq(double(data),param_preproc);
    data = data_tmp.data_denoise;

% preProcessing
    X = data;
    % X = preprocessingData(data,'mean','squared l2-norm');

%% infer spams or nmf

    if strcmp(meth,'spams')

        param_spams.mode        = 2;
        param_spams.lambda      = 0.2;
        param_spams.iter        = 1;
%         param_spams.numThreads  = 1;                                  % number of threads
%         param_spams.batchsize   = 256;
%         param_spams.clean       = 1;
%         param_spams.modeD       = 0;
%         param_spams.posD        = 0;
%         param_spams.posAlpha    = 1;
%         param_spams.gamma1      = 0;
%         param_spams.gamma2      = 0;
%         param_spams.lambda2     = 0;        
%         param_spams.pos         = 1;
        [Coeff,Dict,ReconsDict] = infer_spams(X,1,param_spams);

    else if strcmp(meth,'nmf')
            param_nmf = [];
            param_nmf.iter        = 1000;
            param_nmf.dis         = 0;
            param_nmf.residual    = 1e-4;
            param_nmf.tof         = 1e-4;
            param_nmf.posD        = 1;
            param_nmf.distance    = 'kl'; %'kl', 'ls'
            param_nmf.beta        = 0.1;
            param_nmf.orthogonalA = 0;
            param_nmf.orthogonalB = 0;
            param_nmf.submethodType = 'nmfnnls';
            param_nmf.orthogonal  = [param_nmf.orthogonalA,param_nmf.orthogonalB];
            [Coeffm,Dicts,ReconsDict,ReconsSeq] = infer_nmf(data,1,param_nmf);
        end
    end


    %% exract single cells

    param.Tsimilarity   = 'crosscorrelation';                  % Choose a distance between gaussian 'jointDistance', 'MahalanobisDistance' or 'crosscorrelation'
    param.actThr         = 10;                               % Only considers the image patterns that have been activated more the actThr frames
    param.permutation    = 1;                               % '1' choose a random permutation of ReconsDict, '0' not
    % param.lastWavelet    = 5;                             %  the last wavelet level from ReconsDict(:,:,i) max = 5 min = 1 (only wavelet_denoising)
    % param.prevWavelet    = 4;                             %  the previous last wavelet level from ReconsDict(:,:,i) max = 5 min = 1 (only wavelet_denoising)
    % param.ThrDenoising   = 4;                             % considers only pixels bigger than TheDenoising times the standard deviation of the noise (wavelet_denoising and image_smoothing_powerlaw)
    param.AreaThr        = 30;                              % consider segments biggers than AreaThr
    % param.MajorALength   = 5;                             % length of the major axis
    param.percentageInt  = 0;                               % Intensity percentil
    param.distCentroids  = 16;                              % distance between local maximas
    param.SimThr         = 0.2;                           % similarity measure to fuse two segments
    param.lambda         = 0;%0.05;                         % used for inferring the activation patterns
    param.sizeFilter     = [10 10];                         % size of the filter (image_smoothing and image_smoothing_powerlaw)
    param.sigmaS         = 3;                               % sigma of the gaussian filter
    param.thrSmooth      = 0.2;                       % considers only pixels bigger than thrSmooth times the maximum intensity of the smoothed image
    param.flag_watershed = 0;                               % image segmentation based on nearest points using local maximas ('0') or using watershed ('1')
    param.power          = 30;                           % parameter used for the powerlaw
    param.type = 'image_smoothing';

    cNeurons = [];
    if ~isempty(data)
        if ~isempty(ReconsDict)
            disp('Extracting single cells ... '),
            if ~isempty(cNeurons)
                 % ask if they want to update the current extracted cells                 
                 [cNeurons] = extractSingleCells(Coeff,ReconsDict,cNeurons,param);
            else                
                 [cNeurons] = extractSingleCells(Coeff,ReconsDict,[],param);
            end

            if ~isempty(cNeurons)                
                disp('Extracted Single Cells is done.')
                disp('Computing Temporal Evolution for the extracted single Cells')
                [Coeff_single,Dicts_single,ReconsDict_single] = inferCoeffSingleCells(X,cNeurons,param);

                %% filter neurons that are not used for reconstructing
                [cNeurons,Coeff_single,Dicts_single,ReconsDict_single] = proneNeurons(cNeurons,Coeff_single,Dicts_single,ReconsDict_single);

                disp('Temporal Evolution is computed')    
            else
                disp('No found single cells.')
            end

        else
             msgbox('Please run Matrix Factorization algorithm.','Error','error','modal');
        end
    else
        msgbox('Please load a sequence before setting the number of basis functions.','Error','error','modal');
    end

%     if ~isfield(Im,'cNeurons')
%         Im.cNeurons = cNeurons;
%     end
%     cNeurons_all{i}=cNeurons;
% end
[idFig] = displaySingleCells(double(data1),cNeurons);
% displayMultipleDetectedSignals(data,Coeff_single,cNeurons,8,2);
