function [idFig,dispImg] = displaySelectedCellsIdentification(data,cNeurons,idxNeurons,idFig)
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
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
