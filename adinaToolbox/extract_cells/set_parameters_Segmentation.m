function [this] = set_parameters_Segmentation(this,param)
%  DESCRIPTION:
%    parse the current parameters (this) with the new parameters specified
%    in param
%
%  USAGE:
%
%    [this] = set_parameters_Segmentation(this,param)
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
        
        case 'actthr'
            this.actThr             =  max(param{i+1},0);
            
        case 'permutation'
            this.permutation        =  param{i+1};
            if this.permutation ~= 0
                this.permutation = 1;
            end
        case 'lambda'
            this.lambda             =  max(param{i+1},0);
        case 'lastwavelet'
            this.lastWavelet        =  min(max(param{i+1},1),5);
            
        case 'prevwavelet'
            this.prevWavelet        =  min(max(param{i+1},1),5);
            
        case 'thrdenoising'
            this.ThrDenoising       =  max(param{i+1},0);
            
        case 'areathr'
            this.AreaThr            =  max(param{i+1},0);
            
        case 'majoralength'
            this.MajorALength       =  max(param{i+1},0);
            
        case 'percentageint'
            this.percentageInt      =  min(max(param{i+1},0),1);
            
        case 'distcentroids'
            this.distCentroids      =  max(param{i+1},0);
            
        case 'type'
            this.type               =  param{i+1};
            
        case 'tsimilarity'
            this.Tsimilarity        =  param{i+1};
            
        case 'simthr'
            this.SimThr             =  min(max(param{i+1},0),1);
            
        case 'sizefilter'
            if length(param{i+1}) ==1
               this.sizeFilter         =  round([max(param{i+1},1),max(param{i+1},1)]);
            end
        case 'sigmas'
            this.sigmaS             =  max(param{i+1},0);
            
        case 'thrsmooth'
            this.thrSmooth          =  min(max(param{i+1},0),1);
            
        case 'resetneurons'
            this.resetNeurons        =  param{i+1};
            if this.resetNeurons ~= 0
                this.resetNeurons = 1;
            end
            
        case 'power'
            this.power        =  max(param{i+1},1);
            
        case 'flag_watershed'
            this.flag_watershed       =  param{i+1};
            if this.flag_watershed ~= 0
                this.flag_watgershed = 1;
            end

    end;
    
    i = i+2;
end;