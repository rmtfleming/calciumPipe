function [W,imOut] = waveletDenoisingDraguhn2(im0)
% DESCRIPTION:
%    compute wavelet transform
%
% USAGE:   
%        [W,imOut] = waveletDenoisingDraguhn2(im0)
%
%
% ARGUMENTS:
%         im0 : raw data (nrows x ncols x frames)
%
% OUTPUT: 
%         W: wavelet transform or each order up to 5
%         imOut: smoothed image up to scale order
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".


h = @(i)([1 zeros(1,i-1) 4 zeros(1,i-1) 6 zeros(1,i-1) 4 zeros(1,i-1) 1]/16);

imOut(:,:,6) = im0;
imPrev = im0;
for i = 1:4,
    h_cur = conv2(h(i),h(i)','full');
    [nrows_h,ncols_h]  = size(h_cur);
    %imCur = conv2(padarray(imPrev,[round((nrows_h-1)/2),round((ncols_h-1)/2)],'symmetric'),h_cur,'valid');
    imCur = imfilter(imPrev,h_cur,'symmetric','corr','same'); %(imfilter(imPrev,h_cur,'symmetric','corr','same')+imfilter(imPrev,h_cur,'symmetric','conv','same'))/2;
    W(:,:,5-i) = imPrev - imCur;
    imOut(:,:,5-i) = imCur;
    imPrev = imCur;
end
