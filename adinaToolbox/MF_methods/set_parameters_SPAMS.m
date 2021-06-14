function [this] = set_parameters_SPAMS(this,param)
%  DESCRIPTION:
%    parse the current parameters (this) with the new parameters specified
%    in param
%
%  USAGE:
%
%    [this] = set_parameters_SPAMS(this,param)
%
%  ARGUMENTS:
% 
%    this        : current set of parameters
%    param       : set of parameters to modify
%
%  OUTPUT:
%   this        :  modified parameteres needed for the image enhancement
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

if nargin == 2,
    if mod(length(param), 2) ~=0
        error('Parse_input_parameter: Input parameters must be given in pairs (name and value)!');
    end;
else
    error('Incorrect number of input parameters')
end

i = 1;

while (i <= length(param))
    
    if ischar(param{i+1})
        param{i+1} = str2double(param{i+1});
    end;
    
    switch lower(param{i})
        
        case 'mode'
            this.mode             =  param{i+1};
            
        case 'lambda'
            this.lambda           =  max(param{i+1},0);
            
        case 'numthreads'
            this.numThreads       =  param{i+1};
            if this.numThreads ~=1,
                this.numThreads       =  max(param{i+1},1);
            end
            
        case 'bacthsize'
            this.batchsize        = max(param{i+1},1);
            
        case 'clean'
            this.clean            = param{i+1};
            if this.clean ~=0
                this.clean = 1;
            end
            
        case 'moded'
            this.modeD            = param{i+1};
            if this.modeD ~= 0
                this.modeD = 1;
            end
            
        case 'posd'
            this.posD             = param{i+1};
            if this.posD ~= 0
                this.posD = 1;
            end
             
        case 'posalpha'
            this.posAlpha         = param{i+1};
            this.pos              = param{i+1};
            
            if this.pos ~= 0
                this.pos = 1;
                this.posAlpha = 1;
            end
            
        case 'gamma1'
            this.gamma1           = max(param{i+1},0);
        
        case 'gamma2' 
            this.gamma2           = max(param{i+1},0);
            
        case 'lambda2' 
            this.lambda2          = max(param{i+1},0);
            
        case 'iter' 
            this.iter             = param{i+1};
            if this.iter ~=-1
                this.iter         = max(abs(param{i+1}),1);
            end
            
        case 'pos' 
            this.posAlpha         = param{i+1};
            this.pos              = param{i+1};
            
            if this.pos ~= 0
                this.pos = 1;
                this.posAlpha = 1;
            end
  
    end;
    
    i = i+2;
end;