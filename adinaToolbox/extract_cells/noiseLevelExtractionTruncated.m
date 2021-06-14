function [std_noise] = noiseLevelExtractionTruncated(im0)
% DESCRIPTION:
%     extracts the noise level of a image
%
% USAGE:   
%       [std_noise] = noiseLevelExtractionTruncated(im0)
%
% ARGUMENTS:
%        im0 : input image
%       
% OUTPUT: 
%       std_noise: noise level
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
std_all = sqrt(mean(im0(:).^2));
idx = find(im0>(-3*std_all) &  im0<(3*std_all));
std_noise = sqrt(mean(im0(idx(:)).^2)); 