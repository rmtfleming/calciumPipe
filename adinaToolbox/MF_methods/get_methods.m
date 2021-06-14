function [listTool] = get_methods(pathDir)
%  DESCRIPTION:
%    check is the directory of a matrix factorization is in the current
%    path
%
%  USAGE:
%
%    [listTool] = get_methods(pathDir)
%
%  ARGUMENTS:
% 
%    pathDir       :  current path of the package
%
%  OUTPUT:
%   this        : list of matrix factorization methods
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".
if nargin == 0 
    pathDir = [];
elseif nargin ~=1
    error('get_methods function. Incorrect number of input parameters')
end

listTool = {};
t = 1;

if isdir([pathDir '/SPAMS/'])
    listTool{t} = 'SPAMS';
    t= t+1;
end

if isdir([pathDir '/NMF/'])
    listTool{t} = 'NMF';
    t = t+1;
end

if t ==1
    listTool{t} = 'No methods';
end

