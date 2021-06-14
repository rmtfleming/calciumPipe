function [idFig,dispImg] = displaySingleCellsIdentification(data,cNeurons,idFig)

if nargin == 2,
    idFig = [];
elseif nargin ~=3
    error('Number of input parameters is incorrect')
end

colord=[
        0       0.4000  0
        1.0000  0       0
        0       0.7500  0.7500
        0.7500  0       0.7500
        0.8     0.5     0
        0       0       0.5
        0       0.85    0];
maxImg = max(data,[],3);
dispImg = maxImg - min(maxImg(:));
dispImg = dispImg./max(dispImg(:));
imagesc(dispImg), axis square, axis off;
  
for j = 1:length(cNeurons),
    hold on,plot(cNeurons(j).obj.mu(1),cNeurons(j).obj.mu(2),['.k'],'MarkerSize',7)
    hold on, contour(cNeurons(j).imOrig,1,'Color',colord(mod(j,7)+1,:),'LineWidth',2);
    text(cNeurons(j).obj.mu(1)+2,cNeurons(j).obj.mu(2)+2,num2str(j),'FontSize',12,'Color','black','FontWeight','bold');
end
axis off;
