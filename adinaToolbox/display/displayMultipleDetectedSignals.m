function displayMultipleDetectedSignals(data,Coeff_single,cNeurons,framerate,signalsPerPlot)
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
num_signals = size(Coeff_single,2);
nFrames = size(Coeff_single,1);
tlims = [0 nFrames/framerate];
dt = 1/framerate;

load b.mat a;
colord=[0       0       1.0000
        0       0.4000  0
        1.0000  0       0
        0       0.7500  0.7500
        0.7500  0       0.7500
        0.8     0.5     0
        0       0       0.5
        0       0.85    0];

numIni = round(size(a,1)/(min(signalsPerPlot,num_signals)+1)); 
colord = a(numIni:numIni:end,:);

for i = 1:signalsPerPlot:num_signals-signalsPerPlot
%     figure, title(['Detected signals (' num2str(i) ' to ' num2str(i+signalsPerPlot-1) ')']);
    subplot(1,2,1)
    figure,
    %set(gca,'ColorOrder',colord);
    maxImg = max(data,[],3);
    dispImg = maxImg - min(maxImg(:));
    dispImg = dispImg./max(dispImg(:));
    imagesc(dispImg), colormap gray, axis square, axis off;
    %imshow(cat(3,dispImg,dispImg,dispImg),[]), axis square, axis off;
    for j = i:1:i+signalsPerPlot-1
        hold on,plot(cNeurons(j).obj.mu(2),cNeurons(j).obj.mu(1),'.w','MarkerSize',7)
        hold on, contour(cNeurons(j).imMask,1,'Color',colord(mod(j-1,size(colord,1))+1,:),'LineWidth',2);
        text(cNeurons(j).obj.mu(2)+2,cNeurons(j).obj.mu(1)+2,num2str(j),'FontSize',12,'Color','w');
    end    
    axis off, axis square, axis equal
    
    
    ax = figure ; %subplot(1,2,2);
    set(gcf,'DefaultAxesColorOrder',colord)
    complot(Coeff_single(:,(i:i+signalsPerPlot-1))',(i:i+signalsPerPlot-1),1/framerate)
    formataxes
    xlim(tlims)
    xlabel('Time (s)','FontAngle','i')
    ylabel('Cell number','FontAngle','i')
    set(gcf,'Color','w','PaperPositionMode','auto')
    %set(gca,'yticklabel',num2str(fliplr(i:i+signalsPerPlot-1)'))
    
%     axes('Position',get(ax,'Position'),'XAxisLocation','top','Color','none')
%     xt = get(ax,'XTick');
%     xlim(tlims)
%     formataxes
%     set(gca,'XTick',xt,'XTickLabel',num2str(xt'/dt, '%15.0f'))
    
%set(gca,'YTick',[], ...
%        'XTick',xt,'XTickLabel',num2str(xt'/dt, '%15.0f'))
    xlabel('Frame number')
%     axes(ax)
    box on
end

if num_signals <= signalsPerPlot
    i = 1;
else
    i = i+signalsPerPlot;
end
% figure, title(['Detected signals (' num2str(i) ' to ' num2str(i+num_signals-1) ')']);
figure , %subplot(1,2,1)
%set(gca,'ColorOrder',colord);
maxImg = max(data,[],3);
dispImg = maxImg - min(maxImg(:));
dispImg = dispImg./max(dispImg(:));
imagesc(dispImg), colormap gray, axis square, axis off;
for j = i:1:num_signals
    hold on,plot(cNeurons(j).obj.mu(2),cNeurons(j).obj.mu(1),'.w','MarkerSize',7)
    hold on, contour(cNeurons(j).imMask,1,'Color',colord(mod(j-1,size(colord,1))+1,:),'LineWidth',2);
    text(cNeurons(j).obj.mu(2)+2,cNeurons(j).obj.mu(1)+2,num2str(j),'FontSize',12,'Color','w');
end
axis off, axis square, axis equal


ax = figure ; % subplot(1,2,2);
set(gcf,'DefaultAxesColorOrder',colord)
complot(Coeff_single(:,(i:num_signals))',(i:num_signals),1/framerate)
formataxes
xlim(tlims)
xlabel('Time (s)','FontAngle','i','FontSize',18)
ylabel('Cell number','FontAngle','i','FontSize',18)
set(gcf,'Color','w','PaperPositionMode','auto')
%set(gca,'yticklabel',num2str(fliplr(i:num_signals)'))

% axes('Position',get(ax,'Position'),'XAxisLocation','top','Color','none')
% xt = get(ax,'XTick');
% xlim(tlims)
% formataxes
% set(gca,'XTick',xt,'XTickLabel',num2str(xt'/dt, '%15.0f'))
% xlabel('Frame number')
% axes(ax)
box on




function complot(sig, ICuse, dt)

ICuse = 1:length(ICuse);
for i = 1:length(ICuse)
    zsig(i, :) =sig(ICuse(i),:)./max(sig(ICuse(i),:)); % zscore(sig(ICuse(i),:));
end

alpha = mean(max(zsig')-min(zsig'));
if islogical(zsig)
    alpha = 1.5*alpha;
end

zsig2 = zsig;
alpha = 0.5;
for i = 1:size(ICuse,2)
    zsig2(i,:) = zsig(i,:) - alpha*(i-1)*ones(size(zsig(1,:)));
end

tvec = (1:size(zsig,2))*dt;
if islogical(zsig)
    plot(tvec, zsig2','LineWidth',1)
else
    plot(tvec, zsig2','LineWidth',1)
end
axis tight

%set(gca,'YTick',(-size(zsig,1)+1)*alpha:alpha:0);
%set(gca,'YTicklabel',fliplr(ICuse));

% set(gca,'YTick',[ (-size(zsig,1)+1)*alpha:alpha*10:0 0]);
% set(gca,'YTicklabel',fliplr([ICuse(1:10:end) ICuse(end)]));

set(gca,'YTick',fliplr([0:-alpha*10:(-size(zsig,1)+1)*alpha (-size(zsig,1)+1)*alpha]));
set(gca,'YTicklabel',fliplr([ICuse(1:10:end) ICuse(end)]));




function formataxes

set(gca,'FontSize',14,'FontName','Helvetica','LineWidth',1,'TickLength',[1,1]*.02,'tickdir','out')
set(gcf,'Color','w','PaperPositionMode','auto')
