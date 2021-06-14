function detrend_ca_movie(Im, path)
% This function uses Emperical Mode Decomposition to detrend calcium movies
% INPUTS:
%     Im: a structure containing 
%         Im.DataMat: a n x m matrix of your dataset 
%             (where n is the number of pixels of one image and m numberof images)
%         Im.MeanFrame: the meanframe of the dataset
%     path: where the result should be saved
% 
% by S. Hachi: 13/05/2016


data = Im.DataMat;
meanframe = Im.MeanFrame;
mask = ones(Im.h,Im.w);

% extract the fluorescence trae from the movie
f = fluoTraceExtraction(data, meanframe, mask);

% EMD
[imf,ort,nbits] = emd(f');
parfor i = 1 : size(data,1)
    data_emd(i,:) = data(i,:) - imf(end,:);     % subtract residue to detrend
end

% reshape matrix from 2D to 3D
parfor i = 1 : size(data,2)
    img = reshape(data_emd(:,i), Im.h, Im.w);
    data_detrend(:,:,i) = img;
end

% save detrended time series
for K=1:Im.T
    imwrite(uint16(data_detrend(:,:,K)), path, 'WriteMode', 'append');
end