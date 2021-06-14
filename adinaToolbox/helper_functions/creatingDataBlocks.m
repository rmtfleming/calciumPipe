function [] = creatingDataBlocks(params)
% DESCRIPTION:
%    create blocks specified in params for extractingROI_memory script
%
% USAGE:   
%       [] = creatingDataBlocks(params)
%
%
% ARGUMENTS
%         params: parameters to divide the raw data into blocks
%
% OUTPUT: 
%
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
    save(['./' params(1).main.files(1).name(1:end-4) '_tmpBlock_' num2str(id,'%04d') '.mat'],'data_block','-v7.3');
end


 

