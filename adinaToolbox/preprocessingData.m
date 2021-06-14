function [preData] = preprocessingData(data,substract,normalize)
%  DESCRIPTION:
%    function that normalizes each frame in data that subtract an operator
%    defined by "subtract" and divide it by the operator "normalize"
%
%  USAGE:
%
%    [preData] = preprocessingData(data,substract,normalize)
%
%  ARGUMENTS:
%
%    data     : 3D matrix
%    substract: operator to subtract each frame (none,mean,min,...)
%    normalize: operator to normalize each frame (none, l2-norm,...)
%  
%  OUTPUT:
%    predata: pre-processed data
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".




switch nargin
    case 1
        substract = 'non';
        normalize = 'non';
        [nrows,ncols,nFrames] = size(data);
    case 2
        normalize ='non';
        [nrows,ncols,nFrames] = size(data);
    case 3
        [nrows,ncols,nFrames] = size(data);
    otherwise
        error('Number of input parameters are incorrect')
end
preData = reshape(data,[nrows*ncols nFrames]);

switch substract
    case 'mean'
        preData = bsxfun(@minus,preData,mean(preData)); %preData - repmat(mean(preData),[nrows*ncols 1]);
    case 'min'
        preData = bsxfun(@minus,preData,min(preData)); 
    case 'median'
        preData = bsxfun(@minus,preData,median(preData));
    case 'non'
        preData = preData;
    otherwise
        error('Incorrect substract method for the data')
end
disp(['Subtracting each frame by ' substract '.'])

switch normalize
    case 'l2-norm'
        preData = bsxfun(@rdivide,preData,sqrt(sum(preData.^2)));
    case 'squared l2-norm'
        preData = bsxfun(@rdivide,preData,sum(preData.^2));
    case 'l1-norm'
        preData = bsxfun(@rdivide,preData,sum(abs(preData)));      
    case 'non'
        preData = preData;
    otherwise
        error('Incorrect normalization method for the data')
        
end
disp(['Normalizing each frame by ' normalize '.'])

preData = double(reshape(preData,[nrows ncols nFrames]));


