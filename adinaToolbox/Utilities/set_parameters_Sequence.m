function [this] = set_parameters_Sequence(this,param)
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
        case 'non'
            this{1} = ' ';
        case 'image resize'          
            if ischar(param{i+1})
                param{i+1} = str2double(param{i+1}); 
            end;
            value = min(max(param{i+1},0.0625),2);
            this{2} = num2str(value);            
        case 'inverted intensities'
            if strcmpi(param{i+1},'true') || strcmp(param{i+1},'1')
                this{3} = 'true';
            else
                this{3} = 'false'
            end
        case 'wavelet denoising'
            if ischar(param{i+1})
                param{i+1} = str2double(param{i+1}); 
            end;
            value = min(max(param{i+1},1),5);
            this{4} = num2str(value);   
    end;    
    i = i+2;
end;