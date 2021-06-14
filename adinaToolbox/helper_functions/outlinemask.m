function [ maskb ] = outlinemask( im,mask )
%OUTLINEMASK Summary of this function goes here
%   Detailed explanation goes here

maskb=mask<0;
maskb=medfilt2(single(imdilate(imerode(maskb,ones(3)),ones(3))));
maskb= edge(maskb,'log');



figure;
imagesc(single(im));colormap gray;colorbar;hold on;imagesc(maskb,'AlphaData',maskb);hold off;
figure
imagesc(mask);colorbar