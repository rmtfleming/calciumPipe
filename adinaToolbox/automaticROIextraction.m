function [InferCells] = automaticROIextraction(params)
%  DESCRIPTION:
%    Extracts cell soma of each block and the parameters defined in params. 
%    Each element of params corresponds to a block that the algorithm will process. 
%
%  USAGE:
%
%    [InferCells] = automaticROIextraction(params)
%
%  ARGUMENTS:
%
%    params     :  length(params) indicates the number of blocks to process. 
%                  The elements in params indicates the block to process and the current parameters need to segment cell somas.
%         .main :  parameters used for the cell detection
%         .time :  subset of frames from the original image sequence
%         .space:  define the spatial slice of the original image
%         coordinate
%  
%  OUTPUT:
%    InferCells:   contains a list of detected neurons, one for each block.
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".



nBlocks = length(params);
nrows   = params(1).main.nrows;
ncols   = params(1).main.ncols;
nFrames = params(1).main.nFrames;

[data] = extractingDataGreenChannel(params(1).main.dirBase,params(1).main.files,params(1).main.step,nrows,ncols,nFrames);
data = data(params(1).main.BorderSize+1:end-params(1).main.BorderSize,params(1).main.BorderSize+1:end-params(1).main.BorderSize,:);

for id = 1:nBlocks,
    data_block = double(data(params(id).spatial.rWin,params(id).spatial.cWin,params(id).time.rWSeq));
    disp('Entering Data Normalization')
    seq_param_tmp = extractingNormalizedSeq(data_block,params(id).main.prepro);
    neurons = processingCellDetectionBlock(seq_param_tmp.data_denoise(:,:,params(id).time.rMar),round(length(params(id).time.rSeq)*0.05),params(id).spatial.rMarS,params(id).spatial.cMarS,params(id));
    InferCells(id).neurons = neurons;
    InferCells(id).params = params(id);
end



