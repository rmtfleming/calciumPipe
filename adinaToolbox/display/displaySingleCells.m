function [idFig] = displaySingleCells(data,cNeurons,idFig)
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
if nargin == 2,
    idFig = figure;
elseif nargin ~=3
    error('Number of input parameters is incorrect')
end

figure(idFig);
maxImg = max(data,[],3);
dispImg = maxImg - min(maxImg(:));
dispImg = dispImg./max(dispImg(:));
imagesc(dispImg), axis off, axis square, axis equal, %colormap 'winter'
   
for j = 1:length(cNeurons),
%     hold on,plot(cNeurons(j).obj.mu(2),cNeurons(j).obj.mu(1),'.k','MarkerSize',7)
    hold on, contour(cNeurons(j).imMask,1,'Color','r','LineWidth',2); colormap 'gray'
    drawnow,
%     text(cNeurons(j).obj.mu(2)+2,cNeurons(j).obj.mu(1)+2,num2str(j),'FontSize',8,'Color','y');
%     rectangle('Position', [275-10 150-15 20 30], 'EdgeColor', 'b');
end
end

function circle(x,y,r)
    %x and y are the coordinates of the center of the circle
    %r is the radius of the circle
    %0.01 is the angle step, bigger values will draw the circle faster but
    %you might notice imperfections (not very smooth)
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    plot(x+xp,y+yp,'LineWidth',2);
end