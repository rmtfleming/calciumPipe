% This script extracts the activity of a neuronal population in calcium imaging data
% it includes segmentation, fluorescence traces extraction and spike train
% inference.
%
% Created by: S. Hachi, 20-10-2015


%% parameters
crop_data       = 0;    % To analyse a smaller field of view
data_split4     = 0;    % split data into 4 parts in case memory not enough for entire dataset
calcTH          = 0;    % calculate TH positivity
removeNoisySig  = 0;    % Remove noisy signals
runAci          = 0;    % Run aci for morphological analysis
getTransients   = 0;    % find calcium transients in fluorescence traces
saveData        = 0;    % Save segmentation results in MAT files
casa            = 0;    % save results in excel to use in CaSiAn tool

% Load data in .tif or .mat
splitFname   = strsplit(fname,'.');
imname       = splitFname{1};
dataName     = [path fname];
datadir      = [path 'mat' filesep];
resultsdir   = [path 'results' filesep]; 
CaSAdir      = [path 'CaSA' filesep];

% create folder to store data in .mat format
if ~exist ('datadir', 'dir')
    mkdir(datadir);
end
if ~exist([datadir imname '.mat'], 'file')
    Im = load_calcium_data(datadir, imname, dataName, fname);
else
    load([datadir imname '.mat'])
end

%% crop data

if crop_data == 1
    [Im, data3D] = crop_calcium_data(Im, imname);
    save([datadir imname '_cropped'],'Im','-v7.3')
    return
else
    data3D = reshape3D(Im);
end

% Clear Im to save RAM
h = Im.h;
w = Im.w;
meanFrame = Im.MeanFrame;
clear Im

%% divide data into 4 parts

if data_split4 ==1
    split_data_into4(data3D, [path imname], 8);
    return
end


%% cell detection


if ~exist([datadir imname '_cNeurons.mat'],'file')
    
    [cNeurons, param, Coeff_s, param_preproc] = detectNeurons(data3D(:,:,1:20));
    % Display cell detection results
    [idFig] = displaySingleCells(reshape(meanFrame, h, w),cNeurons);
    
else
    load([datadir imname '_cNeurons.mat']);
end

%% Calculate TH psitivity
if calcTH
    immunoFileInfo = dir([path '*immuno*.tif']);

    imDAPI = imread([path immunoFileInfo.name],'Index',3);
    imTH = imread([path immunoFileInfo.name],'Index',1);

    [th_pos] = th_positivity(cNeurons, imDAPI, imTH);

    for i =1:length(cNeurons)
        cNeurons(i).th_pos = th_pos(i);
    end
end


%% fluorescence signal extraction

if ~exist([datadir imname '_time_traces.mat'],'file')
    % erode the segments for a more accurate fluorescence trace measurement
    cNeurons1 = cNeurons;
    se = strel('disk',5);
    for i=1:length(cNeurons)
        cNeurons1(i).imMask = imerode(cNeurons1(i).imMask,se);
    end
    time_traces = raw_fluorescence_traces(data3D, cNeurons1);                    % extract fluorescence signals
else
    load([datadir imname '_time_traces.mat']);
end

% clear cNeurons1

%% fluorescence signal post-processing

% parameters for deltaF_F calculation
fps = 2;                    % images per second
halfwindow      = 10;
quantile        = 0.2;
step            = 5;
sigma           = 3;
sizeGauss       = 10;

% Calculation of deltaF/F
[dF_F_all, Fbase] = deltaF_F(time_traces, halfwindow, quantile, step, sigma);        % calcultate relative changes in fluorescence dF-F0/F0
dF_F_all = dF_F_all';
dF_F_all_z1 = dF_F_all;

a = find(isnan(mean(dF_F_all_z1,2))==1);                 % find signals = NaN
dF_F_all_z1(a,:) = [];
cNeurons(a) = [];

% remove noisy signals after varMap
if removeNoisySig
    
    ranges = range(dF_F_all_z1');
    badSigFromRanges = find(ranges < 0.2);
    figure,
    for i =1:length(badSigFromRanges)
        plot(z1(dF_F_all_z1(badSigFromRanges(i),:))+i), hold on
    end

    goodSigIdx = 1:length(cNeurons);
    goodSigIdx(badSigFromRanges) = [];
    cNeurons = cNeurons(goodSigIdx);
    dF_F_all_z1 = dF_F_all_z1(goodSigIdx,:);

    figure,
    for i =1:length(goodSigIdx)
        plot(z1(dF_F_all_z1(i,:))+i), hold on
    end

    activity_ratio = length(goodSigIdx)/length(cNeurons);
end

%%% calcultate ISI
event = cell(1,size(dF_F_all_z1,1));
for i = 1 : size(dF_F_all_z1,1)
    events{i} = EventDetection2(dF_F_all_z1(i,:));
end
[signalClass, ~] = calculate_ISI(dF_F_all_z1, fps);

% sort the signals by number of events, frequency or h_positivity
method = 'isi';               % sort signals by 'freq', 'number_spikes' or 'th'
[dF_F_all_z1, spike_times, firing_rates, cNeurons, th_percent] = ...
    sort_signals(dF_F_all_z1, cNeurons, fps, method, signalClass);

%% plot fluorescence traces

if ~exist ('resultsdir', 'dir')
    mkdir(resultsdir);
end

plot_spikes = 0;
sig_per_fig = 16;                               % set number of signals to show in one axis
pathFig     = [resultsdir imname];
order       = '';

meanFrame = imread('results\AVG_wellG6_3_crop4-1.tif');
meanFrame = double(meanFrame);

plot_fluorescence_tracesAccordingToClasses(cNeurons, meanFrame, fps, dF_F_all_z1, ...
    sig_per_fig, pathFig, sort(signalClass, 'descend'))

if saveData
    save([datadir imname '_cNeurons'],'cNeurons', '-v7.3');
    save([datadir imname '_time_traces'],'time_traces', '-v7.3');
    save([datadir imname '_activity_ratio'],'activity_ratio', '-v7.3');
end

%% create an excel file for CaSA tool

if casa
    if ~exist ('CaSAdir', 'dir')
        mkdir(CaSAdir);
    end

    T = array2table(time_traces');
    time = [1:size(time_traces,2)]';
    T1 = table(time);
    T = [T1 T];
    writetable(T,[CaSAdir imname '.xls']);
end

%% cross correlation

if runAci

    cMapsDir = [path 'cMaps' filesep];
    if ~exist ('cMapsDir', 'dir')
        mkdir(cMapsDir);
    end
    
    if numel(size(data3D)) == 3
        data3D = permute(data3D,[1 2 4 3]);
    end

    cMaps = calculate_corr_maps(data3D, time_traces);    
    for i = 1 : size(cMaps,4)
        imwrite(uint8(cMaps(:,:,:,i)), [cMapsDir 'cMaps' num2str(i) '.tif']);
    end
end

%% get transients

if getTransients
    for i = 1 : size(dF_F_good,1)
        events{i} = EventDetection(dF_F_good(i,:));    
    end
    s.dF_cell       = dF_F_good;
    s.Spikes_cell   = events;
    s.fps           = 10;

    [A,G,rise_time,fall_time,CV] = getTransients(s);

    for i = 1 : size(dF_F_good,1)
        events{i} = EventDetection(dF_F_good(i,:));
        [~, n{i}] = findpeaks(dF_F_good(i,:),'MinPeakHeight', mean(dF_F_good(i,:)), 'MinPeakDistance',5);
    end
    clc
    spks = zeros(76,400);
    onsets = zeros(76,400);
    for i = 1 : size(dF_F_good,1)
        spks(i,n{i})=1;
        onsets(i,events{i})=1;
    end
end

