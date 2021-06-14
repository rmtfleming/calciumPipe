clear all, clc;

path = 'D:\LCSB calcium imaging data\Test Neo camera_21-8-14\3 lane uFlu plate_25-8-14\image_5_well M24_25-8-14\stack\';
FileTif1= [path 'stack1_8bits.tif'];
FileTif2=[path 'stack2_8bits.tif'];
FileTif3=[path 'stack_MMStack_2.ome.tif'];
outputFileName = [path 'full_stack.tif'];
delimit_area=0;
numStacks = 2;

% if delimit_area==1
%     image = imread(FileTif,'Index',1);
%     h=figure; imshow(image);
%     pause(7);
%     [x, y] = getpts(h);
%     x = round(x);
%     y = round(y);
%     close
% else
%     x=[1 2048];
%     y=[1 2048];    
% end

for i=1:numStacks
    if i==1
        FileTif=FileTif1;
        InfoImage=imfinfo(FileTif);
        NumberImages=length(InfoImage);
        for j=1:NumberImages
            image = imread(FileTif,'Index',j);
            FinalImage1(:,:,j)= image;%(y(1):y(2),x(1):x(2));
        end
    else if i==2
            FileTif=FileTif2;
            InfoImage=imfinfo(FileTif);
            NumberImages=length(InfoImage);
            for j=1:NumberImages
                image = imread(FileTif,'Index',j);
                FinalImage2(:,:,j)= image; %(y(1):y(2),x(1):x(2));
            end
        else
            FileTif=FileTif3;
            InfoImage=imfinfo(FileTif);
            NumberImages=length(InfoImage);
            for j=1:NumberImages
                image = imread(FileTif,'Index',j);
                FinalImage3(:,:,j)= image; %(y(1):y(2),x(1):x(2));
            end
        end
    end
end

FinalImage=cat(3,FinalImage1, FinalImage2);
%%
img=FinalImage;
for K=1:length(img(1, 1, :))
   imwrite(img(:, :, K), outputFileName, 'WriteMode', 'append');
end