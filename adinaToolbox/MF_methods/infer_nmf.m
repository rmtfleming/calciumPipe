function [Coeffm,Dicts,ReconsDict,ReconsSeq] = infer_nmf(data,nBasis,param)
% 
% DESCRIPTION:
%     This function calls one of the submethods of NMF 
%
% USAGE:
%   [Coeffm, Dicts [,ReconsDict,ReconSeq]] = infer_nmf(X,nBasis,param);
%
% ARGUMENTS:
%   data     double m (npixels) x n (nframes) matrix   (input signals)
%   nBasis   number of spatial/image patterns to factorize data    
%   param.submethodType: STRING
%       * 'nmfnnls'        % For more information, type help nmfnnls
%       * 'nmfrule'        % For more information, type help nmfrule
%       * 'sparsenmfnnls'  % For more information, type help sparsenmfnnls
%       * 'seminmfnnls'    % For more information, type help seminmfnnls
%       * 'seminmfrule'    % For more information, type help seminmfrule
%       * 'convexnmfrule'  % For more information, type help convexnmfrule
%       * 'orthnmfrule'    % For more information, type help orthnmfrule
%       * 'wnmfrule'       % For more information, type help wnmfrule
%
%   param: struct
%       + NMF  : type help method, where method is one of the submethods mentioned above  
%
% OUTPUT: 
%         Dicts      : double m (npixels) x p (nBasis) matrix   (vectorized dictionary)
%         Coeffm     : double n (nframes) x p (nBasis) matrix   (temporal activation)
%         ReconsDict : double nrows x ncols x nBasis            (dictionary)
%         ReconsSeq  : double nrows x ncols x nFrames           (reconstructed sequence)
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

if nargin ~=3
    error('Incorrect number of input parameters')
end

[nrows,ncols,nFrames] = size(data);
K = nBasis;

data           =  reshape(data,[nrows*ncols nFrames]);
options        =  param;
[Dicts,Coeffm] = feval(param.submethodType,data,K,options);
Coeffm = single(Coeffm');

if nargout >= 3,
    for i = 1:nBasis,
        ReconsDict(:,:,i) = reshape(Dicts(:,i),[nrows ncols]);
    end
    ReconsDict = single(ReconsDict);
    if nargout >= 4,
        ReconsSeq = single(reshape(Dicts*Coeffm',[nrows ncols nFrames]));
    end
end