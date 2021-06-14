clear all, clc;

tic

path = 'D:\LCSB calcium imaging data\Test Zyla\image_2\stack\';
FileTif= [path 'Original_stack_8bits.tif'];
outputFileName = [path 'small_stack.tif'];
InfoImage=imfinfo(FileTif);
NumberImages=length(InfoImage);

 % delimit interesting area
image = imread(FileTif,'Index',1);
h=figure; imagesc(image); colormap('hot');
pause(5);
[x, y] = getpts(h);
x = round(x);
y = round(y);
close all;

for i=1:NumberImages
       image = imread(FileTif,'Index',i);
       FinalImage1(:,:,i)= image(y(1):y(1)+599,x(1):x(1)+599);
end

% FileTif='D:\Microscope Cameras files\Test Zyla\image_6\stack\stack\stack2.tif';
% InfoImage=imfinfo(FileTif);
% NumberImages=length(InfoImage);
% 
% for i=1:NumberImages
%        image = imread(FileTif,'Index',i);
%        FinalImage2(:,:,i)= image(y(1):y(2),x(1):x(2));
% end

%%
% FinalImage=cat(3,FinalImage1,FinalImage2);
for K=1:length(FinalImage1(1, 1, :))
   imwrite(FinalImage1(:, :, K), outputFileName, 'WriteMode', 'append');
end

toc