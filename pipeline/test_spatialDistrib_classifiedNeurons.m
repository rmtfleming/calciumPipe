load('signalsWT.mat')
load('signalsMUT.mat')
load('sigClass.mat')

for i = 1 : length(signalsWT)
    splitt = strsplit(signalsWT(i).fileName,'_');
    signalsWT(i).fileName = [splitt{1} '_' splitt{2} '_'];
end

for i = 1 : length(signalsMUT)
    splitt = strsplit(signalsMUT(i).fileName,'_');
    signalsMUT(i).fileName = [splitt{1} '_' splitt{2} '_'];
end

p = uigetdir;

files = dirrec(p,'.mat');

% find the time traces files for WT and MUT
[jF, jFN, jFIm] = find_trace_files(files);
files_timetraces    = files(jF);
files_cNeurons      = files(jFN);
files_IM            = files(jFIm);

% Plot signals of class
fps                 = 5;
numSig_perFig       = 10;

T = 400;

Pl.xlims  =[1 T];
Pl.nticks = 4;
Pl        = PlotParams(Pl);
Pl.fs     = 20;

%%

% add class number at the end
for i = 1 : length(sigClass)
    sigClass{i}(:,end+1) = i;
end 

% merge all classes in one array
sigClassAll = cell2mat(sigClass');

% find WT elements
colWT = find(sigClassAll(:,end-1)==1);
sigClassAllWT = sigClassAll(colWT,:);
% sort by file
sigClassAllWT = sortrows(sigClassAllWT,size(sigClassAllWT,2) - 3);

% find MUT elements
colMUT = find(sigClassAll(:,end-1)==2);
sigClassAllMUT = sigClassAll(colMUT,:);
% sort by file
sigClassAllMUT = sortrows(sigClassAllMUT,size(sigClassAllMUT,2) - 3);

%%
[nb, elem] = find_nb_occurences(sigClassAllWT(:,end-3));
pathFig     = 'Z:\SBG\experimentalProgram\Experiments2016\ELM-Exp-46_analysis\analysis\controlplate12\spatialDistrib\';

[nb, elem] = find_nb_occurences(sigClassAllWT(:,end-3));

for i = 2 : length(elem)
    
    col = find(sigClassAllWT(:,end-3) == elem(i));
    neuronsToPlot = sigClassAllWT(col, end-2);
    fileName = [signalsWT(elem(i)).fileName];
    
    fileMatches = strfind(files_IM,fileName);
    x = files_IM(~cellfun(@isempty,fileMatches));
    load(x{1})
    meanFrame = Im.MeanFrame;
    clear Im
    
    fileMatches = strfind(files_cNeurons,fileName);
    x = files_timetraces(~cellfun(@isempty,fileMatches));
    load(x{1})
    
    fileMatches = strfind(files_cNeurons,fileName);
    x = files_cNeurons(~cellfun(@isempty,fileMatches));
    load(x{1})
    
    cNeurons = cNeurons(neuronsToPlot);
    time_traces = time_traces(neuronsToPlot,:);
    
    pathFigFull = [pathFig fileName];
    plot_fluorescence_tracesAccordingToClasses(cNeurons, meanFrame, fps, time_traces,...
        length(cNeurons), pathFigFull, sigClassAllWT(col, end))
end

%%
[nb, elem] = find_nb_occurences(sigClassAllMUT(:,end-3));

for i = 1 : length(elem)
    
    col = find(sigClassAllMUT(:,end-3) == elem(i));
    neuronsToPlot = sigClassAllMUT(col, end-2);
    fileName = [signalsMUT(elem(i)).fileName];
    
    fileMatches = strfind(files_IM,fileName);
    x = files_IM(~cellfun(@isempty,fileMatches));
    load(x{1})
    meanFrame = Im.MeanFrame;
    clear Im
    
    fileMatches = strfind(files_cNeurons,fileName);
    x = files_timetraces(~cellfun(@isempty,fileMatches));
    load(x{1})
    
    fileMatches = strfind(files_cNeurons,fileName);
    x = files_cNeurons(~cellfun(@isempty,fileMatches));
    load(x{1})
    
    cNeurons = cNeurons(neuronsToPlot);
    time_traces = time_traces(neuronsToPlot,:);
    
    pathFigFull = [pathFig fileName];
    plot_fluorescence_tracesAccordingToClasses(cNeurons, meanFrame, fps, time_traces,...
        70, pathFigFull, sigClassAllMUT(col, end))
end


