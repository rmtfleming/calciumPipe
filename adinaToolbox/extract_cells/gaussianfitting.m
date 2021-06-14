function obj = gaussianfitting(X)
% DESCRIPTION:
%     compute an approximation of a 2D-Gaussian model
%
% USAGE:   
%       obj = gaussianfitting(X)
%
% ARGUMENTS:
%      X: spatial indices
%       
% OUTPUT: 
%      obj : mean and covariance of X.
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;

obj.mu = mean(X,1);
obj.Sigma = cov(X(:,1),X(:,2));