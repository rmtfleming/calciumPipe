function [W,imOut] = waveletDenoisingDraguhn3DScale(Vol,scale)
% DESCRIPTION:
%    compute wavelet transform
%
% USAGE:   
%        [W,imOut] = waveletDenoisingDraguhn3D(Vol)
%
%
% ARGUMENTS:
%         Vol : raw data (nrows x ncols x frames)
%         scale : maximum order of the wavelet transform
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

% wavelets mask
h{1} = [1 4 6 4 1]/16;
h{2} = [ 1 0 4 0 6 0 4 0 1]/16;
h{3} = [ 1 0 0 0 4 0 0 0 6 0 0 0 4 0 0 0 1]/16;
h{4} = [ 1 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 1]/16;
h{5} = [ 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]/16;

if length(scale)~= 1 && length(scale) ~= 3
    error('Scale should be an scalar or a vector of three elements');
else
    if length(scale) == 1;
        scale = repmat(scale,[1 3]);
    end
end

if nargout ==2
    W = zeros(size(Vol,1),size(Vol,2),size(Vol,3),scale(3)+1);
    imOut(:,:,:,scale(3)+1) = Vol;
end
imPrev = Vol;
W = single(zeros(size(Vol,1),size(Vol,2),size(Vol,3)));
for i = 1:max(scale),
    imCur = imPrev;
    for z = 1:3,
        switch z
            case 1
                if i <=scale(1),
                    imCur = imfilter(imCur,h{i},'symmetric','conv','same'); %(imfilter(imPrev,h_cur,'symmetric','corr','same')+imfilter(imPrev,h_cur,'symmetric','conv','same'))/2;
                end
            case 2
                if i <=scale(2),
                    imCur = imfilter(imCur,h{i}','symmetric','conv','same');
                end
            case 3
                if i <=scale(3),
                    h_3 = [];
                    h_3(1,1,:) = h{i};
                    imCur = imfilter(imCur,h_3,'symmetric','conv','same');
                end
        end
        %imCur = imfilter(imPrev,h_cur,'symmetric','conv','same'); %(imfilter(imPrev,h_cur,'symmetric','corr','same')+imfilter(imPrev,h_cur,'symmetric','conv','same'))/2;
    end
%     imCur = imfilter(imPrev,h_cur,'symmetric','conv','same'); %(imfilter(imPrev,h_cur,'symmetric','corr','same')+imfilter(imPrev,h_cur,'symmetric','conv','same'))/2;
    W(:,:,:) = imPrev - imCur;
    if nargout ==2
        imOut(:,:,:,scale+1-i) = imCur;
    end
    imPrev = imCur;
end
