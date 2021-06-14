function [nucleus_mask] = nucleus_segm(imN)
% 
% This function segments the nuclei 
% 
%   INPUTS: imN:        DAPI image
%   OUTPUTS: cellbw4:  nucleus mask
% 
% Created by S. Hachi: 25/07/2016


% image preprocessing
% imN_preproc = function_ImPreProcessing(imN, 20, 'Nucleus');

% segment nuclei
prm.method = 'thrs';
prm.thrsth = 0.2;
prm.split = 0;
prm.smoothim.method = 'eed';
prm.smoothim.eed.kappa = 0.1;
[cellbw3,wat,imsegmout,prmout] = cellsegm.segmct(double(imN),1,500,'prm',prm);
% cellsegm.show(cellbw3,3);title('Cell segmentation by THRS');axis off;

% improving the results by splitting of cells
n = prmout.minvolvox;
h = [0.5 0.5 1.5];
splitth = 1;
nucleus_mask = cellsegm.splitcells(cellbw3,splitth,n,h);
