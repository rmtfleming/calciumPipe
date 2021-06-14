function plot_fluo_traces_class(sigClass, classToPlot, paramsWT_MUT, time_traces, fps, ...
    numSig_perFig, pathFig)

% This function plots the fluorescence traces and their corresponding spike
% trains
% INPUTS:
%     cNeurons      = the mask of each segmented neuron
%     Im            = the calcium data
%     traces         = a matrix nNeurons x ntimePoints of fluorescence traces
%     numSig_perFig = number of plots per figure
%     pathFig       = path to save the figure


load('colord.mat');

T = size(time_traces,2);
for i = 2:T+1
    tvec(1,i-1)=(1/fps)*i;
end


Pl.xlims  =[1 T];
Pl.nticks = 4;
Pl        = PlotParams(Pl);
Pl.fs     = 10;

wid = 0.8;
hei = 0.8;
left = 0.1;
bottom = 0.1;

colord=[0       0       1.0000
        0       0.4000  0
        1.0000  0       0
        0       0.7500  0.7500
        0.7500  0       0.7500
        0.8     0.5     0
        0       0       0.5
        0       0.85    0];

    
numIni = round(size(a,1)/(length(sigClass{classToPlot})+1)); 
colord = a(numIni:numIni:end,:);
        
for k = 1 : floor(length(sigClass{classToPlot})/numSig_perFig)+1
    
    if k*numSig_perFig > length(sigClass{classToPlot})
        num_sig = length(sigClass{classToPlot});
    else
        num_sig = k*numSig_perFig;
    end

    h_signals = figure;
    subplot('position',[left bottom wid hei]); hold on
    axis([-6 (1/fps)*T+7 0 2*(numSig_perFig+1)]);
    set(gca,'yTick',[], 'YColor','w');
    xlabel('Time (sec)','FontSize',Pl.fs);
    ylabel('Fluorescence (a.u.)','FontSize',Pl.fs, 'Color', 'k');
    set(gca,'fontsize',Pl.fs)
    
    kk = numSig_perFig*(k-1);
    for i = kk+1 : num_sig
        
        % plot fluo signal
        figure(h_signals)
        plot(tvec, z1(detrend(time_traces(sigClass{classToPlot}(i,end-1),:))) -1 + 2*(i-kk),'Color',colord(mod(i-1,size(colord,1))+1,:), 'LineWidth', 1);hold on        
        wh=[3 8];
        set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');

    end

    wh=[5 13];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');

    % save signal figures on drive
    if ~strcmp(pathFig, '')
        figure(h_signals)
        figname=[pathFig '_traces' num2str(k)];
        print('-dtiff',figname);
    end
end

