function [InferCells] = automaticROIextraction_memory(params)
%  DESCRIPTION:
%    Extracts cell soma of each block and the parameters defined in params. 
%    Each element of params corresponds to a block that the algorithm will process. 
%
%  USAGE:
%
%    [InferCells,file] = automaticROIextraction(params,file)
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
%    InferCells:   contains the detected neurons and also meta information
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

nrowsBlock = length(params.spatial.rWin);
ncolsBlock = length(params.spatial.cWin);
nFrames    = length(params.time.rWSeq);
data = zeros(nrowsBlock,ncolsBlock,nFrames);

for i_data = 1:nFrames,
    tmp = imread([params.main.dirBase params.main.files(params.main.idFiles(params.time.rWSeq(i_data))).name],params.main.idSlice(params.time.rWSeq(i_data)));
    tmp = tmp(params.main.BorderSize+1:end-params.main.BorderSize,params.main.BorderSize+1:end-params.main.BorderSize,:);
    data(:,:,i_data) = double(tmp(params.spatial.rWin,params.spatial.cWin));
end
disp('Entering Data Normalization')
seq_param_tmp = extractingNormalizedSeq(data,params.main.prepro);
neurons = processingCellDetectionBlock(seq_param_tmp.data_denoise(:,:,params.time.rMar),round(length(params.time.rSeq)*0.05),params.spatial.rMarS,params.spatial.cMarS,params);
InferCells.neurons = neurons;
InferCells.params = params;

