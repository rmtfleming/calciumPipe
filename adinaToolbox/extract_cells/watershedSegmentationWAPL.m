function [image_inv,segmented_dict] = watershedSegmentationWAPL(dict,flag_closing)
% DESCRIPTION:
%    compute watershed of an image
%
% USAGE:   
%       segmented_dict = watershedSegmentationWAPL(dict,flag_closing)
%
%
% ARGUMENTS:
%         dict : input image
%         flag_closing : smooth operation before watershed
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
    image = dict(:,:,i);
    if flag_closing
        image_morph = imopen(image,ones(3));
    else
        image_morph = image;
    end
    image_inv = -image_morph;
    image_inv(image_inv==0) = -Inf;
    L = watershed(image_inv);
    segmented_dict(:,:,i) = L;
end
