function [params] = get_params(method,submethod)
%  DESCRIPTION:
%   return the list of parameters for a spefic matrix factorization method
%   and type of inference
%
%  USAGE:
%
%    [params] = get_params(method,submethod)
%
%  ARGUMENTS:
%   method      : matrix factorization method
%   submethod   : type of inference or objective to minimize
%
%  OUTPUT:
%    param       : current set of parameters
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".
if nargin ~=2,
    error('Error at get_params. Number of input parameters is incorrect.')
end

 switch lower(method)
     case 'spams'
         params = initial_parameter_SPAMS;
         params.methodType = 'SPAMS';
         params.submethodType = submethod;
         params.ListSub = feval(submethod);
     case 'nmf'
         params.all = initial_parameter_NMF;
         params.all.methodType = 'NMF';
         params.all.submethodType = submethod;
         params.Type = 'NMF';
         params.subType = submethod;
         switch lower(submethod)
             case 'nmfnnls',
                 params.subList = {'iter','dis','residual','tof'};
             case 'nmfrule',
                 params.subList = {'iter','dis','residual','tof','distance'};
             case 'sparsenmfnnls',
                 params.subList = {'iter','dis','residual','tof','beta'};
             case 'seminmfnnls',
                 params.subList = {'iter','dis','residual','tof'};
             case 'seminmfrule',
                 params.subList = {'iter','dis','residual','tof'};
             case 'convexnmfrule',
                 params.subList = {'iter','dis','residual','tof'};
             case 'orthnmfrule',
                 params.subList = {'iter','dis','residual','tof','orthogonalA','orthogonalB'};
             case 'wnmfrule',
                 params.subList = {'iter','dis','residual','tof','distance'};
             otherwise
                 error('Incorrect NMF submethod')
         end
         
     otherwise
         params = [];
 end