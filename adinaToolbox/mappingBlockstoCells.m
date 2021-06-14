function [roi] = mappingBlockstoCells(InferCells)
% DESCRIPTION:
%   Fuse all the neurons detected independently in each block and map the
%   ROI to the original image coordinates,
%
% USAGE:   [roi] = mappingBlockstoCells(InferCells)%
%
% ARGUMENTS:
%       InferCells     contains the list of cells detected at each block
%
% OUTPUT: 
%       roi:
%       roi(i).Img : cell appearance (weighted mask) of the i-th cell
%       roi(i).indices  : pixel location of the i-th cell%
%   
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

%% Extracting the number of files in the folder

InferCellGlobal.neurons =  struct('obj',[],'imMask',[],'im',[],'imOrig',[],'Area',[],'maxF',[],'X',[],'X_v',[],'Centroid',[],'time',[]);
curNNeurons = 1; %Current number of Neurons + 1
for i = 1:length(InferCells),
    for iNeurons = 1:length(InferCells(i).neurons),
	InferCells(i).params.main.seg.Tsimilarity = 'crosscorrelation';
        if ~isempty(InferCells(i).neurons{1})
            tmpImg = zeros(InferCells(i).params.main.nrows,InferCells(i).params.main.ncols);
            tmpImg(InferCells(i).params.spatial.rWin,InferCells(i).params.spatial.cWin) = InferCells(i).neurons{iNeurons}.imOrig; 
            tmpImg = padarray(tmpImg,[InferCells(i).params.main.BorderSize, InferCells(i).params.main.BorderSize]);	
            tmpNeuron = struct('obj',[],'imMask',tmpImg>0,'im',tmpImg,'imOrig',tmpImg,'Area',sum(tmpImg(:)>0),'maxF',max(tmpImg(:)),'X',[],'X_v',[],'Centroid',[],'time',[]);
            if curNNeurons==1,
               InferCellGlobal.neurons(curNNeurons) = tmpNeuron;
               curNNeurons = curNNeurons + 1;
            else
               InferCellGlobal.neurons = identifyDistinctCells(InferCellGlobal.neurons,tmpNeuron,InferCells(i).params.main.seg);
            end
        end
    end    
end

for i = 1:length(InferCellGlobal.neurons)
    roi(i).Img = InferCellGlobal.neurons(i).imOrig;
    roi(i).indices = find(InferCellGlobal.neurons(i).imOrig>0);
end
