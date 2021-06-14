function [nucleus, Seeds_x, Seeds_y, th_int] = calculate_th_positivity(imC, nucleus_mask)
% 
% This function segments the nuclei 
% INPUTS: imC:          TUJ1 image
%         nuclus_mask:  mask of nuclei after segmentation
%         cNeurons:   neurons detected from calcium time-series
% OUTPUTS: th_intensity:  fluorescence intensity for each nucleus
% 
% Created by S. Hachi: 25/07/2016


% image preprocessing
imC_preproc = function_ImPreProcessing(imC, 50, 'Cell');

% get nucleus boundaries
[Boundary, Segment, SeedNumber] = bwboundaries(nucleus_mask,'nohole');
Seeds_x =zeros(1,SeedNumber);
Seeds_y =zeros(1,SeedNumber);

% get the center of the seeds.
for i =1:1:SeedNumber
    [x y] = find(Segment==i);
    Seeds_x(i) =round(mean(x));
    Seeds_y(i) =round(mean(y));
end

se = strel('disk',10);
for i = 1 : length(Seeds_x)
    nucleus(i).boundary = Boundary{i};
    imm = zeros(size(imC,1),size(imC,2));
    imm(Segment==i)=1;
    nucleus(i).mask = imm;
    imm_dilate = imdilate(imm,se);
    imMask = imm_dilate(:).*double(imC_preproc(:));
    nucleus(i).imMask = reshape(imMask,size(imC,1),size(imC,2));
    nnz = nonzeros(imMask);
    th_int(i) = mean(nnz);
end