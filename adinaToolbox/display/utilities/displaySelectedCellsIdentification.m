function [idFig,dispImg] = displaySelectedCellsIdentification(data,cNeurons,idxNeurons,idFig)

if nargin == 3,
    idFig = figure;
elseif nargin ~=4
    error('Number of input parameters is incorrect')
end

figure(idFig);
maxImg = max(data,[],3);
dispImg = maxImg - min(maxImg(:));
dispImg = dispImg./max(dispImg(:));
imagesc(dispImg),
   
for j = idxNeurons,
    hold on,plot(cNeurons(j).obj.mu(2),cNeurons(j).obj.mu(1),'.k','MarkerSize',7)
    hold on, contour(cNeurons(j).imMask,'Color','w');
    text(cNeurons(j).obj.mu(2)+2,cNeurons(j).obj.mu(1)+2,num2str(j),'FontSize',12,'Color','black','FontWeight','bold');
end
axis off
