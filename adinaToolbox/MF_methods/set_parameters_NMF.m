function [this] = set_parameters_NMF(this,param)
%  DESCRIPTION:
%    parse the current parameters (this) with the new parameters specified
%    in param
%
%  USAGE:
%
%    [this] = set_parameters_NMF(this,param)
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
       
    switch lower(param{i})
        
        case 'iter'
            if ischar(param{i+1})
                param{i+1} = str2double(param{i+1});
            end;
            this.iter             =  max(param{i+1},0);
                      
        case 'dis'
            if ischar(param{i+1})
                param{i+1} = str2double(param{i+1});
            end;
            this.dis              = param{i+1};
            if this.dis ~=0
                this.dis = 1;
            end
            
        case 'residual'
            if ischar(param{i+1})
                param{i+1} = str2double(param{i+1});
            end;
            this.residual         = max(param{i+1},1e-16);
            
        case 'tof'
            if ischar(param{i+1})
                param{i+1} = str2double(param{i+1});
            end;
            this.tof              =  max(param{i+1},1e-16);
             
        case 'distance'
            if strcmpi(param{i+1},'kl')
                this.distance     = 'kl';
            else
                this.distance     = 'ls';
            end
            
        case 'orthogonala'
            if ischar(param{i+1})
                param{i+1} = str2double(param{i+1});
            end;
            this.orthogonalA            = param{i+1};
            if this.orthogonalA ~=0
                this.orthogonalA = 1;
            end
            this.orthogonal(1) = this.orthogonalA;
        
        case 'orthogonalb'
            if ischar(param{i+1})
                param{i+1} = str2double(param{i+1});
            end;
            this.orthogonalB            = param{i+1};
            if this.orthogonalB ~=0
                this.orthogonalB = 1;
            end
            this.orthogonal(2) = this.orthogonalB;
    end;
    
    i = i+2;
end;