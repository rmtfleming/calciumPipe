function [cost] = jointDistance(muC,sigmaC,muP,sigmaP)
% DESCRIPTION:
%     compute the joint distance on a location
%
% USAGE:   
%       [cost] = jointDistance(muC,sigmaC,muP,sigmaP)
%
% ARGUMENTS:
%        muC: centroid of the current Cell
%        sigmaC: covariance of the current cell
%        muP: centroid of a new Cell
%        sigmaP: covariance of new cell
%       
% OUTPUT: 
%        cost : distance
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
    cost = 1/(sqrt(2*pi*det(sigmaC+sigmaP)))*exp(-0.5*(muP-muC)'*inv(sigmaC+sigmaP)*(muP-muC));
end
