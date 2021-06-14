function [Coeffm,Dicts,ReconsDict,ReconsSeq] = inferGroupedCells(data,nBasis,method,param)
% DESCRIPTION:
%      This function calls one of the modalities: SPAMS or NMF, to
%      decompose the image sequence into two low rank matrices.
%
% USAGE:   [Coeffm, Dicts [,ReconsDict,ReconSeq]] = inferGroupedCells(X,nBasis[,param]);
%          param is optional
%
%
% ARGUMENTS:
%         data     double m (npixels) x n (nframes) matrix   (input signals)
%         nBasis   number of spatial/image patterns to factorize data    
%         method: STRING
%           +SPAMS:
%                  method = 'SPAMS' that calls infer_spams using
%                  initial_parameters_SPAM and modified by param
%
%           +NMF: 
%                  method: STRING can be:
%                      * 'nmfnnls'        % For more information, type help nmfnnls
%                      * 'nmfrule'        % For more information, type help nmfrule
%                      * 'sparsenmfnnls'  % For more information, type help sparsenmfnnls
%                      * 'seminmfnnls'    % For more information, type help seminmfnnls
%                      * 'seminmfrule'    % For more information, type help seminmfrule
%                      * 'convexnmfrule'  % For more information, type help convexnmfrule
%                      * 'orthnmfrule'    % For more information, type help orthnmfrule
%                      * 'wnmfrule'       % For more information, type help wnmfrule
%         param: struct
%               + SPAMS: type help infer_spams
%               + NMF  : type help method, where method is one of the submethods mentioned above  
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


if nargin <3,
    method = 'SPAMS';
    param = feval(['initial_parameter_' method]);
    param.type = 'infer_spams';
elseif nargin == 3,
    if strcmp(method,'SPAMS')
        param = feval(['initial_parameter_' method]);        
        param.type = 'infer_spams';
    elseif regexp(method,'nmf')
       [params] = get_params('NMF',method);
        param = params.all;      
        param.type = 'infer_nmf';
    else
      error('MF_methods: Incorrect Matrix Factorization method.')
    end
elseif nargin == 4
    if strcmpi(method,'SPAMS')    
        param_init = feval(['initial_parameter_' method]);
        param = feval(['set_parameters_' method],param_init,param);      
        param.type = 'infer_spams';
    elseif regexp(method,'nmf')        
        [params] = get_params('NMF',method);
        param = feval('set_parameters_NMF',params.all,param);      
        param.type = 'infer_nmf';
    else
        error('MF_methods: Incorrect Matrix Factorization method.')
    end
else
    error('Incorrect number of input parameters')
end

[Coeffm,Dicts,ReconsDict] = feval(param.type,data,nBasis,param);
[nrows ncols nFrames] = size(data);
if nargout >= 3,
    for i = 1:nBasis,
        ReconsDict(:,:,i) = reshape(Dicts(:,i),[nrows ncols]);
    end
    if nargout >= 4,
        ReconsSeq = reshape(Dicts*Coeffm',[nrows ncols nFrames]);
    end
end
