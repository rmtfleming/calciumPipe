function plotSequence(X,Coeff,Dicts,param)
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;

if nargin <4;
    initF = 1;
    endF  = size(X,3);
    step = 1;
elseif nargin == 4
    lengthSeq = size(X,3);
    initF = min(max(param.initF,1),lengthSeq);
    endF  = min(max(param.endF,initF),lengthSeq);
    step  = min(max(param.step,1),endF-initF+1);
end


h1 = figure;


[nrows,ncols,nFrames] = size(X); 
[nFrames, nBasis] = size(Coeff);

if size(Dicts,1) ~= nrows*ncols,
    Dicts = reshape(Dicts,[nrows*ncols nBasis]);
end

reconsImage = reshape(Dicts*Coeff',[nrows ncols nFrames]);


maxX = max(X(:));
minX = min(X(:));
maxXR = max(reconsImage(:));
minXR = min(reconsImage(:));

rangePlot = (-1:nBasis) +0.5;
for i = initF:step:endF,
    figure(h1), subplot(1,3,1), imagesc(Coeff), hold on,  plot(rangePlot,i*ones(1,length(rangePlot)), 'r'); hold off; title('Temporal Evolution of Cell Coefficients'), drawnow
    figure(h1), subplot(1,3,2), imagesc(X(:,:,i)),caxis([minX maxX]),axis square, axis off, colorbar, title('Raw Data')
    figure(h1), subplot(1,3,3), imagesc(reconsImage(:,:,i)),caxis([minXR maxXR]),colorbar, axis square, axis off, title('Reconstructed Data'), 
    pause(0.1)
end


