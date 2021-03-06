function plot_fluorescence_traces(cNeurons, Im, F_tab, spike_times, order, vec_method, numSig_perFig, pathFig, plot_spikes)

% This function plots the fluorescence traces and their corresponding spike
% trains

load('colord.mat');
tvec(1,1) = 0;
for i = 2:Im.T+1
    tvec(1,i-1)=(1/Im.fps)*i;
end

if strcmp(order,'firing')
    col_name = 'firing rate';
    col_val = roundn(mean(vec_method'),-2);
else
    col_name = 'th positivity %';
    col_val = round(vec_method);
end

Pl.xlims  =[1 Im.T];
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

    
numIni = round(size(a,1)/(length(cNeurons)+1)); 
colord = a(numIni:numIni:end,:);

h_segm = figure;
maxImg = max(reshape(Im.MeanFrame, Im.h, Im.w),[],3);
dispImg = maxImg - min(maxImg(:));
dispImg = dispImg./max(dispImg(:));
imagesc(adapthisteq(dispImg)), colormap 'gray', axis off, axis equal
        
for k = 1 : floor(length(cNeurons)/numSig_perFig)+1
    
    if k*numSig_perFig > length(cNeurons)
        num_sig = length(cNeurons);
    else
        num_sig = k*numSig_perFig;
    end

    h_signals = figure;
    subplot('position',[left bottom wid hei]); hold on
    axis([-6 (1/Im.fps)*Im.T+7 0 2*(numSig_perFig+1)]);
    set(gca,'yTick',[], 'YColor','w');
    xlabel('Time (sec)','FontSize',Pl.fs);
    ylabel('Fluorescence (a.u.)','FontSize',Pl.fs, 'Color', 'k');
    set(gca,'fontsize',Pl.fs)
    
    kk = numSig_perFig*(k-1);
    text((1/Im.fps)*Im.T-13,2*(numSig_perFig+1),col_name,'FontSize',Pl.fs, 'Color', 'k'); hold on
    text(-6,2*(numSig_perFig+1),'Cell #','FontSize',Pl.fs, 'Color', 'k'); hold on
    
    for i = (numSig_perFig*(k-1))+1 : num_sig
        
        % plot fluo signal
        figure(h_signals)
        plot(tvec, z1(F_tab(i,:)) -1 + 2*(i-kk),'Color',colord(mod(i-1,size(colord,1))+1,:), 'LineWidth', 1);hold on
        spks_ind = find(spike_times(i,:)==1);
        wh=[3 8];
        set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    
        % plot spikes
        if plot_spikes == 1
            for ii = 1 : length(spks_ind)
                plot([spks_ind(ii)*(1/Im.fps) spks_ind(ii)*(1/Im.fps)], ...
                    [(min(F_tab(i,:))-1.5 + 2*(i-kk)), (max(F_tab(i,:))-2+ 2*(i-kk))],'Color', 'b', 'LineWidth',1.2); hold on
            end
        end
        text(-6,-0.5+2*(i-kk),num2str(i),'FontSize',Pl.fs, 'Color', 'k'); hold on                                   % cell #
        text((1/Im.fps)*Im.T+1,-0.5+2*(i-kk),num2str(col_val(i)),'FontSize',Pl.fs,'Color', 'k'); hold on      % firing rate or th positivity
    
        figure(h_segm)
        hold on, contour(cNeurons(i).imMask,1,'Color',colord(mod(i-1,size(colord,1))+1,:),...
            'LineWidth',0.5);
        axis off, axis square, axis equal
        text(cNeurons(i).obj.mu(2)+2,cNeurons(i).obj.mu(1)+2,num2str(i),'FontSize',6,'Color','y','LineWidth', 3);
    
    end

    wh=[5 13];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');

    % save signal figures on drive
    figure(h_signals)
    figname=[pathFig '_traces' num2str(k)];
    print('-dtiff',figname);
end

% save segmented image on drive
figure(h_segm)
wh=[6 6];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[pathFig '_segm'];
print('-dtiff',figname);
