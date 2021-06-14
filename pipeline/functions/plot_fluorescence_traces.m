function plot_fluorescence_traces(cNeurons, meanFrame, fps, F_tab, spike_times, ...
    order, vec_method, numSig_perFig, pathFig, plot_spikes)

% This function plots the fluorescence traces and their corresponding spike
% trains
% INPUTS:
%     cNeurons      = the mask of each segmented neuron
%     Im            = the calcium data
%     F_tab         = a matrix nNeurons x ntimePoints of fluorescence traces
%     spike_times   = time points where spikes occur for each fluorescence
%                       signal
%     order         = order the signals by 'firing' or default
%     vec_method    = vector of sored firing rates or sorted th_positivity
%     numSig_perFig = number of plots per figure
%     pathFig       = path to save the figure
%     plot_spikes   = 1 to plot spikes, 0 otherwise


load('colord.mat');

T = size(F_tab,2);
tvec(1,1) = 0;
for i = 2:T+1
    tvec(1,i-1)=(1/fps)*i; 
end

if strcmp(order,'firing')
    col_name = 'firing rate';
    col_val = roundn(mean(vec_method'),-2);
else if strcmp(order,'th')
        col_name = 'th positivity %';
        col_val = round(vec_method);
    else
        col_name = '';
        col_val = [];
    end
end

Pl.xlims  =[1 T];
Pl.nticks = 4;
Pl        = PlotParams(Pl);
Pl.fs     = 24;

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

    
numIni = round(size(a,1)/(length(cNeurons)+1)); 
colord = a(numIni:numIni:end,:);

h_segm = figure;
maxImg = max(reshape(meanFrame, size(cNeurons(1).im,1), size(cNeurons(1).im,2)),[],3);
dispImg = maxImg - min(maxImg(:));
dispImg = dispImg./max(dispImg(:));
imagesc(dispImg), colormap 'gray', axis off, axis equal
        
for k = 1 : floor(length(cNeurons)/numSig_perFig)+1
    
    if k*numSig_perFig > length(cNeurons)
        num_sig = length(cNeurons);
    else
        num_sig = k*numSig_perFig;
    end

    h_signals = figure;
    subplot('position',[left bottom wid hei]); hold on
    axis([-0.5 (1/fps)*T+1 0 2*(numSig_perFig+1)]);
    set(gca,'yTick',[], 'YColor','w');
    set(gca,'xTick',[], 'XColor','w');    
%     xlabel('Time (sec)','FontSize',Pl.fs);
%     ylabel('\DeltaF/F','FontSize',Pl.fs, 'Color', 'k');
    set(gca,'fontsize',Pl.fs)
    
    kk = numSig_perFig*(k-1);
    if strcmp(col_name,'')
        text((1/fps)*T-2,2*(numSig_perFig+1),'Cell #','FontSize',Pl.fs, 'Color', 'k'); hold on
    else
        text((1/fps)*T-9,2*(numSig_perFig+1),col_name,'FontSize',Pl.fs, 'Color', 'k'); hold on
        text(-4,2*(numSig_perFig+1),'Cell #','FontSize',Pl.fs, 'Color', 'k'); hold on
    end
    for i = (numSig_perFig*(k-1))+1 : num_sig
        
        % plot fluo signal
        figure(h_signals)
        plot(tvec, z1(F_tab(i,:)) -1 + 2*(i-kk),'Color',colord(mod(i-1,size(colord,1))+1,:), 'LineWidth', 1.5);hold on        
        wh=[3 8];
        set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    
        % plot spikes
        if plot_spikes == 1
            spks_ind = find(spike_times(i,:)==1);
            for ii = 1 : length(spks_ind)
                plot([spks_ind(ii)*(1/fps) spks_ind(ii)*(1/fps)], ...
                    [(min(F_tab(i,:))-1.5 + 2*(i-kk)), (min(F_tab(i,:))-1.2+ 2*(i-kk))],'Color', 'b', 'LineWidth',1.2); hold on
            end
        end
        
        % Display cell number and firing rate or th_positivity
        if ~strcmp(col_name,'')
            text((1/fps)*T+1,-0.5+2*(i-kk),num2str(col_val(i)),'FontSize',Pl.fs,'Color', 'k'); hold on      % firing rate or th positivity
        else
            text((1/fps)*T+2,-0.5+2*(i-kk),num2str(i),'FontSize',18, 'Color', 'k'); hold on                                   % cell #
        end
        
        % plot cell contours
        figure(h_segm)
        hold on, contour(cNeurons(i).imMask,1,'Color',colord(mod(i-1,size(colord,1))+1,:),...
            'LineWidth',2);
        axis off, axis square, axis equal
        text(cNeurons(i).obj.mu(2)+2,cNeurons(i).obj.mu(1)+2,num2str(i),'FontSize',18, 'FontWeight', 'Bold', 'Color','y');
    
    end

    % save signal figures on drive
    if pathFig
        figure(h_signals)
        hold on
        plot([-0.4; -0.4], [(2*(numSig_perFig+1))-0.6; 2*(numSig_perFig+1)], '-k', 'LineWidth', 4);
        plot([-0.5; 9.5], [0; 0], '-k', 'LineWidth', 5);
        t1 = text(0.5, 2*(numSig_perFig+1)-0.3, '40% \DeltaF/F','FontSize',Pl.fs, 'Color', 'k', 'FontWeight','bold');
        t2 = text(-0.5, -0.8, '10 sec','FontSize',Pl.fs, 'Color', 'k', 'FontWeight','bold');
%         set(t1,'Rotation',90);
%         plot()
        wh=[7 9];   %width and height
        set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
        figname=[pathFig '_traces' num2str(k)];
        print('-dtiff',figname);
    end
end

% save segmented image on drive
if pathFig
    figure(h_segm)
    wh=[6 6];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[pathFig '_segm'];
    print('-dtiff',figname);
end
