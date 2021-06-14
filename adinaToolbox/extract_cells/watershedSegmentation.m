function segmented_dict = watershedSegmentation(dict)
% DESCRIPTION:
%    compute watershed of an image
%
% USAGE:   
%       segmented_dict = watershedSegmentation
%
%
% ARGUMENTS:
%         dict : input image
%
% OUTPUT: 
%         segmented_dict: disjoint regions of the current frame
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

segmented_dict = zeros(size(dict));

for i = 1:size(dict,3)
    image = -dict(:,:,i);
    image(image==0) = -Inf;
    image_close = imclose(image,ones(4));
    L = watershed(image_close);
    segmented_dict(:,:,i) = L;
end
