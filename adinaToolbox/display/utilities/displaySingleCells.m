function [idFig] = displaySingleCells(data,cNeurons,idFig)

if nargin == 2,
    idFig = figure;
elseif nargin ~=3
    error('Number of input parameters is incorrect')
end

figure(idFig),
maxImg = max(data,[],3);
dispImg = maxImg - min(maxImg(:));
dispImg = dispImg./max(dispImg(:));
imagesc(dispImg),
   
for j = 1:length(cNeurons),
    hold on,plot(cNeurons(j).obj.mu(2),cNeurons(j).obj.mu(1),'.k','MarkerSize',7)
    hold on, contour(cNeurons(j).imMask,'Color','w','LineWidth',2);
end
    