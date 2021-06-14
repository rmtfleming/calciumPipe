function [this] = set_parameters_preprocessing(this,param)
%  DESCRIPTION:
%    parse the current parameters (this) with the new parameters specified
%    in param
%
%  USAGE:
%
%    [this] = set_parameters_preprocessing(this,param)
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
        
        case 'method'
            this.method           = param{i+1}; 
            
        case 'halfwinmedian'
            this.halfWinMedian    =  max(param{i+1},0);
            
        case 'stepf0'
            this.stepF0           =  max(param{i+1},0);
            
        case 'tempsigma'
            this.tempSigma        = max(param{i+1},0);
                       
        case 'spasigma'
            this.spaSigma         = max(param{i+1},0);
                        
        case 'wavscale'
            this.wavScale         = param{i+1};  
            
        case 'wavscalespace'
            this.wavScaleSpace         = param{i+1}; 
            this.wavScale(1:2)         = param{i+1}*ones(1,2);
            
        case 'wavscaletime'
            this.wavScaleTime         = param{i+1};
            this.wavScaleTime(3)         = param{i+1};
        case 'bordersize'
            this.BorderSize         = param{i+1};
    end;
    
    i = i+2;
end;