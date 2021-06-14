function [A_sorted] = fastExtractingQuantilesSeq(A,halfWindow,quartiles,step)
%  DESCRIPTION:
%   fast extraction of the quantile on a sliding window for computing F_0
%
%  USAGE:
%
%    [A_sorted] = fastExtractingQuantilesSeq(A,halfWindow,quartiles,step);
%
%  ARGUMENTS:
%
%    A               : vectorized raw data (npixels x nFrames)   
%    halfWindow      : half of the size of the sliding windows
%    quantile        : which quantile to use for f0 computation
%    ste             : to speed up f0 estimation, f0 is only computed every this many frames
% 
%  OUTPUT:
%    A_sorted        : approx. i-th quantile on the vectorized matrix A
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

if size(A,3)>1,
    error('A should be a matrix (npixels x nFrames)')
end

 % should be a factor of nrows and ncols
[nrows,nFrames] = size(A);

if step>nFrames,
    error('please choose an step lower than the length of the sequence')
end

A = A';
cPoints = [1:step:nFrames nFrames];
for i = 1:length(cPoints),
    Quan(i,:,:) = quantile(A(max(1,cPoints(i)-halfWindow):min(nFrames,cPoints(i)+halfWindow),:),quartiles)';
end

A_sorted = interp1(cPoints',Quan,(1:nFrames)');
