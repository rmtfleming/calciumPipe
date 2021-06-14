% Load folder
clear, clc, close all

pathh = uigetdir;
pathWT = [pathh filesep 'K7' filesep];
pathMUT = [pathh filesep 'K7M' filesep];
pathRes = [pathh filesep 'results' filesep];

xlabels = {'WT'; 'LRRK2'};
if ~exist ('pathRes', 'dir')
    mkdir(pathRes);
end

kmeans_clust        = 0;
plot_comparison     = 0;
plot_distrib        = 1;
plot_traces         = 0;
run_tSNE            = 0;
significance_test   = 1;

%% Get data from excel sheets

if ~exist([pathWT 'signalsWT.mat'])
    signalsWT = load_signal_param(pathWT);
    save([pathWT 'signalsWT'], 'signalsWT', '-v7.3');
else
    load([pathWT 'signalsWT.mat'])   
end

%% Get mean and std of interspike intervals and other parameters for WT

paramsCaSA = get_mean_std_params(signalsWT, 10);

sigFileIdWT        = paramsCaSA.sigFileId;
isiWT              = paramsCaSA.isi;
ampWT              = paramsCaSA.amp;
SWidthWT           = paramsCaSA.SWidth;
SAreaWT            = paramsCaSA.SArea;
timeToPeakWT       = paramsCaSA.timeToPeak;
releasingRateWT    = paramsCaSA.releasingRate;
removingRateWT     = paramsCaSA.removingRate;

isiWT.data = isiWT.data*2.5;
isiWT.mean = isiWT.mean*2.5;
isiWT.std= isiWT.std*2.5;

SWidthWT.data = SWidthWT.data*2.5;
SWidthWT.mean = SWidthWT.mean*2.5;
SWidthWT.std= SWidthWT.std*2.5;

% remove signals with width < 0
m = find(SWidthWT.mean < 0);
sigFileIdWT(m,:)           = [];
isiWT.mean(m)              = [];        isiWT.std(m)              = [];
ampWT.mean(m)              = [];        ampWT.std(m)              = [];
SWidthWT.mean(m)           = [];        SWidthWT.std(m)           = [];
SAreaWT.mean(m)            = [];        SAreaWT.std(m)            = [];
timeToPeakWT.mean(m)       = [];        timeToPeakWT.std(m)       = [];
releasingRateWT.mean(m)    = [];        releasingRateWT.std(m)    = [];
removingRateWT.mean(m)     = [];        removingRateWT.std(m)     = [];

for i = 1 : length(signalsWT)
    splitt = strsplit(signalsWT(i).fileName,'_');
    signalsWT(i).fileName = [splitt{1} '_' splitt{2} '_'];
end

%% Get data from excel sheets mutants

if ~exist([pathMUT 'signalsMUT.mat'])
    signalsMUT = load_signal_param(pathMUT);
    save([pathMUT 'signalsMUT'], 'signalsMUT', '-v7.3');
else
    load([pathMUT 'signalsMUT.mat'])    
end

%% Get mean and std of interspike intervals and other parameters for MUT

paramsCaSA = get_mean_std_params(signalsMUT, 10);

sigFileIdMUT         = paramsCaSA.sigFileId;
isiMUT               = paramsCaSA.isi;
ampMUT               = paramsCaSA.amp;
SWidthMUT            = paramsCaSA.SWidth;
SAreaMUT             = paramsCaSA.SArea;
timeToPeakMUT        = paramsCaSA.timeToPeak;
releasingRateMUT     = paramsCaSA.releasingRate;
removingRateMUT      = paramsCaSA.removingRate;

isiMUT.data = isiMUT.data*2.5;
isiMUT.mean = isiMUT.mean*2.5;
isiMUT.std= isiMUT.std*2.5;

SWidthMUT.data = SWidthMUT.data*2.5;
SWidthMUT.mean = SWidthMUT.mean*2.5;
SWidthMUT.std= SWidthMUT.std*2.5;

for i = 1 : length(signalsMUT)
    splitt = strsplit(signalsMUT(i).fileName,'_');
    signalsMUT(i).fileName = [splitt{1} '_' splitt{2} '_'];
end

%% classify signals using kmeans or intervals

paramsWT_MUT = [[[isiWT.mean' isiWT.std']; [isiMUT.mean' isiMUT.std']],...
    [[ampWT.mean' ampWT.std']; [ampMUT.mean' ampMUT.std']], ...
    [[SWidthWT.mean' SWidthWT.std']; [SWidthMUT.mean' SWidthMUT.std']],...
    [[SAreaWT.mean' SAreaWT.std']; [SAreaMUT.mean' SAreaMUT.std']],...
    [[timeToPeakWT.mean' SAreaWT.std']; [timeToPeakMUT.mean' SAreaMUT.std']],...
    [[releasingRateWT.mean' SAreaWT.std']; [releasingRateMUT.mean' SAreaMUT.std']],...
    [[removingRateWT.mean' SAreaWT.std']; [removingRateMUT.mean' SAreaMUT.std']],...
    [sigFileIdWT; sigFileIdMUT]];

paramsWT_MUT2 = paramsWT_MUT(:,1:end-2);

% assign 1 for WT and 2 for mutant - add column at the end
a = size(paramsWT_MUT);
paramsWT_MUT(1:length(isiWT.mean), a(end)+1) = 1;
paramsWT_MUT(length(isiWT.mean)+1:length(isiWT.mean)+length(isiMUT.mean), a(end)+1) = 2;

% separate signals into classes
classMeth = '';             % or 'kmeans' if clustering
[paramsWT_MUT]  = classify_sigs(paramsWT_MUT, isiWT, isiMUT, 3, classMeth);

%% separate classes and plot

sigClass = separate_classes(paramsWT_MUT, 10);

% merge classes
sigClassAll = mergeClasses(sigClass);

% generate figure
figure_classificationAndClasses;

%% Compare parameters between WT and mutants from the same class

if plot_comparison
    
    pathBarPlots = [pathRes 'bar_plots' filesep];
    if ~exist ('pathBarPlots', 'dir')
        mkdir(pathBarPlots);
    end
    
    plot_comparison_bars(sigClass, paramsLabels, xlabels, pathBarPlots)
end

%% Ditribution
parameters = {'Interspike interval', 'amplitude', 'spike width', 'spike area', 'time to peak',...
        'releasing rate', 'removing rate', '\tau'};
    
if plot_distrib
    pathDistrib = [pathRes 'ksDensity' filesep];
    if ~exist ('pathDistrib', 'dir')
        mkdir(pathDistrib);
    end
    
    test_decayTime;
    variables = {'isi', 'amp', 'SWidth', 'SArea', 'timeToPeak', 'releasingRate', 'removingRate', 'tau'};
    
    for i = 1 : length(parameters)
        proba_distribution(eval([variables{i} 'WT.mean']),...
            eval([variables{i} 'MUT.mean']), parameters{i}, 1, pathDistrib)
    end
    close all;
end

%% significance test

if significance_test
    
    pathKLDiv = [pathRes 'KLDivergence' filesep];
    if ~exist ('pathKLDiv', 'dir')
        mkdir(pathKLDiv);
    end
    
    % calculate KL divergence
    
    data1 = paramsWT_MUT2(1:length(isiWT.mean),1:2:end-1);
    data2 = paramsWT_MUT2(length(isiWT.mean)+1:end,1:2:end-1);
    
    % distribution
    pathKLDivPercent = [pathKLDiv 'rand_distribution' filesep];
    if ~exist ('pathKLDivPercent', 'dir')
        mkdir(pathKLDivPercent);
    end
    num_iter = 1000;        % number of randomisations
    dist = KLDivergence_random_distribution(data1, data2, num_iter, parameters, 1, pathKLDivPercent);
    dist = KLDivergence_random_distribution(tauWT.mean', tauMUT.mean', 1000,  {'\tau'}, 1, pathKLDivPercent);
    
    % iterations
    rand_percent = .05;
    pathKLDivPercent = [pathKLDiv 'rand_percentage_' num2str(rand_percent) filesep];
    if ~exist ('pathKLDivPercent', 'dir')
        mkdir(pathKLDivPercent);
    end
    dist = KLDivergence_random_iterations(data1, data2, rand_percent, num_iter, parameters, 0, pathKLDivPercent);

end

%% find traces

if plot_traces
    p = uigetdir;
    files = dirrec(p,'.mat');

    % find the time traces files for WT and MUT
    [jF, jFN] = find_trace_files(files);
    files_timetraces    = files(jF);
    files_cNeurons      = files(jFN);

    % Plot signals of class
    fps                 = 5;
    numSig_perFig       = 10;
    
    T = 400;
    
    Pl.xlims  =[1 T];
    Pl.nticks = 4;
    Pl        = PlotParams(Pl);
    Pl.fs     = 20;
    
    for classToPlot = 3 : 3%length(sigClass)
        
        sigClassTraces = zeros(size(sigClass{1,classToPlot},1),400);
        jj = 1;

        for k = 1 : floor(length(sigClass{classToPlot})/numSig_perFig)+1

            if k*numSig_perFig > length(sigClass{classToPlot})
                num_sig = length(sigClass{classToPlot});
            else
                num_sig = k*numSig_perFig;
            end

            h_signals = figure;

            kk = numSig_perFig*(k-1);
            for i = kk+1 : num_sig
                if sigClass{classToPlot}(i,end) == 1
                    fileName = [signalsWT(sigClass{classToPlot}(i,end-2)).fileName];
                    colorSig = 'b';
                else
                    fileName = [signalsMUT(sigClass{classToPlot}(i,end-2)).fileName];
                    colorSig = 'r';
                end
                matches = strfind(files_timetraces,fileName);
                % time traces to load
                x = files_timetraces(find(~cellfun(@isempty,matches)));
                load(x{1})
                tvec=[];
                T = size(time_traces,2);
                for ii = 2:T+1
                    tvec(1,ii-1)=(1/fps)*ii;
                end
                trace = time_traces(sigClass{classToPlot}(i,end-1),:);
                sigClassTraces(jj,1:length(trace)) = trace;
                jj = jj + 1;
                plot(tvec,z1(detrend(trace)) -1 + (i-kk), colorSig, 'LineWidth', 1.5), hold on
            end

            hold on
            plot([0; 0], [0; -0.4], '-k',  [0; 4.4], [-0.4; -0.4], '-k', 'LineWidth', 2)        
            t1 = text(0, -0.6, '5 sec','FontSize',Pl.fs, 'Color', 'k', 'FontWeight','bold');
            t2 = text(-2, -0.4, '40% \DeltaF/F','FontSize',Pl.fs, 'Color', 'k', 'FontWeight','bold');        
            set(t2,'Rotation',90);
            axis off
        end

        sigClassTracesAll{classToPlot} = sigClassTraces;
    end
%     save('sigClassTracesAll','sigClassTracesAll','-v7.3');
%     for i = 1 : size(sigClassTraces,1)
%         [~, n{i}] = findpeaks(sigClassTraces(i,:),'MinPeakHeight', mean(sigClassTraces(i,:)), 'MinPeakDistance',5);
%     end
end

%% t-SNE

if run_tSNE
    A = [paramsWT_MUT(:,1:2:end-3) paramsWT_MUT(:,end-1)];
    %A=A(:,1:2:end);     % take only the mean
    figure,
    ydata = tsne(A(:,1:end-1), A(:,end), 3, [], 50); 
end

