function [Im, data3D] = crop_calcium_data(Im, imname)

% This function allows the user to crop a region in a calcium dataset by
% selecting two points corresponding to the corners of the region to crop


    h=figure; imagesc(reshape(Im.DataMat(:,1), Im.h, Im.w)); colormap('hot');
    pause(5);
    [x, y] = getpts(h);
    x = round(x);
    y = round(y);
    close all;

    for i=1:Im.T
           image = reshape(Im.DataMat(:,i), Im.h, Im.w);
           data3D(:,:,i)= image(y(1):y(2),x(1):x(2));
           im = data3D(:,:,i);
           dataMat_crop(:,i) = im(:);
    end
    Im.DataMat = dataMat_crop;
    Im.MeanFrame = mean(Im.DataMat,2);
    Im.h = size(data3D,1);
    Im.w = size(data3D,2);
    Im.Np = Im.h*Im.w;
    Im.fname = [imname '_cropped.tif'];