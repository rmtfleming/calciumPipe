function [W,imOut] = waveletDenoisingDraguhnScale(im0,scale)
% DESCRIPTION:
%    compute wavelet transform up to scale-th order.
%
% USAGE:   
%        [W,imOut] = waveletDenoisingDraguhnScale(im0,scale)
%
%
% ARGUMENTS:
%         data  : raw data (nrows x ncols x frames)
%         scale : order of the wavelet scale
%
% OUTPUT: 
%         W: wavelet transform or each order up to scale.
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

h{1} = [1 4 6 4 1]/16;
h{2} = [ 1 0 4 0 6 0 4 0 1]/16;
h{3} = [ 1 0 0 0 4 0 0 0 6 0 0 0 4 0 0 0 1]/16;
h{4} = [ 1 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 1]/16;
h{5} = [ 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]/16;

imOut(:,:,6) = im0;
imPrev = im0;
for i = 1:scale,
    h_cur = conv2(h{i},h{i}','full');
    imCur = imfilter(imPrev,h_cur,'symmetric','conv','same'); %(imfilter(imPrev,h_cur,'symmetric','corr','same')+imfilter(imPrev,h_cur,'symmetric','conv','same'))/2;
    W(:,:,scale-i+1) = imPrev - imCur;
    imOut(:,:,scale-i+1) = imCur;
    imPrev = imCur;
end
