% This script allows you to choose to run the electroAnalysis in a batch
% mode. It looks for all the tif files contained in a folder and its
% subfolders and run the analysis for all of them

% clear, close all; clc,

tic,

pathh           = uigetdir;
files           = dirrec(pathh,'.tif');

jjj=1;
for iii=1:1%length(files)
    [path,fname,ext] = fileparts(files{iii});
    path = [path filesep];
    fname = [fname ext];
    electrophysAnalysis;

%     close all
%     clearvars -except files pathh
end

toc
%% calculate average percentage of active neurons

activeNeurons   = 0;

if activeNeurons
    avg = [];
    fileIdx = [];
    activity_avg = 0;
    
    files = dirrec(pathh,'.mat');
    jjj=1;
    for i=1:length(files)
        if strfind(files{i},'_ratio.mat') 
            fileIdx(jjj) = i;
            jjj=jjj+1;
            load(files{i})
            activity_avg = activity_avg + activity_ratio;
        end
    end
    activity_avg = activity_avg/length(fileIdx);
end
