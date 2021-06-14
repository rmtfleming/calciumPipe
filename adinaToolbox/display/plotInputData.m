function plotInputData(X,param)
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
if nargin <2;
    initF = 1;
    endF  = size(X,3);
    step = 1;
elseif nargin == 2
    lengthSeq = size(X,3);
    initF = min(max(param.initF,1),lengthSeq);
    endF  = min(max(param.endF,initF),lengthSeq);
    step  = min(max(param.step,1),endF-initF+1);
end

[nrows,ncols,nFrames] = size(X); 
maxX = max(X(:));
minX = min(X(:));
for i = initF:step:endF,
    imagesc(X(:,:,i)),caxis([minX maxX]),axis square, axis off, colorbar,drawnow
end