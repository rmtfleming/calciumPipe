function f = fluoTraceExtraction(data, meanFrame, mask)

% This function extracts the fluorescence trace from a region of interest (ROI)
% INPUTS:
%     data: calcium imaging data
%     meanframe: mean frame of the data
%     mask: mask corresponding to ROI
%     
% OUTPUTS:
%     f: extracted fluorescent trace

if ndims(data)==3
    w = size(data,1);
    h = size(data,2);
    t = size(data,3);
    dataMat = zeros(w*h,t);
    for j=1:t
        X = data(:,:,j);        
        dataMat(:,j)=X(:);
    end
else
    dataMat = data;
end

a = ones(size(meanFrame));
weighted_ROI = meanFrame.*mask(:);

f = weighted_ROI'*dataMat(:,1:end)/sum(weighted_ROI(:));
% f = detrend(f);
    
end