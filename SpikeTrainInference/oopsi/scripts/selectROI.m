% Run the staup script first

% clear all, clc, close all
tic

path         = '/mnt/siham_hachi/data/spikeTrainInference_data/image_5/';
imname       = 'full_small_stack';
tifname     = [path 'stacks/' imname '.tif'];
figdir      = [path 'figs/'];
datadir     = [path 'mat/'];
fname       = 'full_small_stack';
fig_num     = 2;
if fig_num==1
    F_TH        = load([datadir 'F_TH_noisy_image.mat']);
    prettytif   = [path 'stacks/AVG/' 'MAX_full_stack_TH_red20.tif'];
    prettytiff  = imread(prettytif);
else
    prettytiff = [];
    F_TH=[];
end

% set switches of things to do
LoadTif     = 0;
small_stack = 0;
GetROI      = 1;

%% get raw image data

if LoadTif == 1                                     % get whole movie
    MovInf  = imfinfo(tifname);                     % get number of frames
    Im.T  = numel(MovInf);                          % only alternate frames have functional data
    Im.h  = MovInf(1).Height;
    Im.w  = MovInf(1).Width;
    Im.Np = Im.w*Im.h;
    Im.fname = fname;
    PSF = fspecial('gaussian',7,2.);
    
    Im.DataMat = zeros(Im.w*Im.h,Im.T);             % initialize mat to store movie
    for j=1:Im.T
        X = imread(tifname,j);        
        Im.DataMat(:,j)=X(:);
    end
    Im.MeanFrame=mean(Im.DataMat,2);
    save([datadir fname],'Im','-v7.3')
else
    load([datadir fname])
end
toc

%% Create (F-F_0)/F0 data

for i=1:Im.T
    img=reshape(Im.DataMat(:,i),Im.h,Im.w);
    data(:,:,i)=img(:,:);
end
param = initial_parameter_preprocessing();
param.method = 3;                                   % choose (F-F_0)/F0 method
[seq_param] = extractingNormalizedSeq(data,param);
for i=1:Im.T
    x = seq_param.data_denoise(:,:,i);
    deltaF_data(:,i) = x(:);
end

%% create small stack

if small_stack == 1
    InfoImage=imfinfo(tifname);
    NumberImages=length(InfoImage);

    % delimit area of interest
    image = imread(tifname,'Index',1);
    h=figure; imshow(image);
    pause(5);
    [x, y] = getpts(h);
    x = round(x);
    y = round(y);

    %FinalImage=zeros(nImage,mImage,10);%NumberImages);
    for i=1:NumberImages
       image = imread(tifname,'Index',i);
       FinalImage(:,:,i)= image(y(1):y(2),x(1):x(2));
    end

    img=FinalImage;
    outputFileName = [path 'small_stack2.tif'];
    for K=1:length(img(1, 1, :))
       imwrite(img(:, :, K), outputFileName, 'WriteMode', 'append');
    end
end


%% select roi

if GetROI == 1
    numROI=input('How many ROIs? (max 3): ');
    for i=1:numROI
        figure(1); clf,
        imagesc(reshape(Im.MeanFrame,Im.h,Im.w)')
        colormap('hot')
        title('select roi radius, double click when complete')
        set(gca,'BusyAction','queu','Interruptible','on');
        ellipse0=imellipse;
        wait(ellipse0);
        vertices0=getVertices(ellipse0);
        xdat=vertices0(:,1);
        ydat=vertices0(:,2);
        Im.x0 = min(xdat) + .5*(max(xdat)-min(xdat));
        Im.y0 = min(ydat) + .5*(max(ydat)-min(ydat));
        a = max(xdat)-Im.x0;
        b = max(ydat)-Im.y0;
        Im.radius0=mean([a b]);
        [pixmatx pixmaty] = meshgrid(1:Im.h,1:Im.w);

        if i==1
            BW1 = createMask(ellipse0);
            Im.roi1 = (((pixmatx-Im.x0).^2 + (pixmaty-Im.y0).^2 )<= Im.radius0^2);
            Im.roi_edge1 = edge(uint8(Im.roi1));
            roi_tab(:,:,i)=Im.roi1';
        else if i==2
                BW2 = createMask(ellipse0);
                Im.roi2 = (((pixmatx-Im.x0).^2 + (pixmaty-Im.y0).^2 )<= Im.radius0^2);
                Im.roi_edge2 = edge(uint8(Im.roi2));
                roi_tab(:,:,i)=Im.roi2';
            else if i==3
                    BW3 = createMask(ellipse0);
                    Im.roi3 = (((pixmatx-Im.x0).^2 + (pixmaty-Im.y0).^2 )<= Im.radius0^2);
                    Im.roi_edge3 = edge(uint8(Im.roi2));
                    roi_tab(:,:,i)=Im.roi3';
                end
            end
        end
    end

    % ROIs masks and fluorescence matrices
    for i=1:numROI
        roi1=roi_tab(:,:,i);
        % weighted_ROI= Im.MeanFrame.*roi1(:);
        a = ones(size(Im.MeanFrame));
        weighted = a.*roi1(:);
        if i==1
            mask1 = edge(BW1);    
            Im.MeanPretty(~mask1) = NaN;
            % out = imoverlay(prettytiff, mask1', [1 0 0]);
            % F=weighted_ROI'*Im.DataMat(:,1:end)/sum(weighted_ROI(:)); 
            F = weighted'*deltaF_data(:,1:end);
            F = detrend(F); 
            % F = F-min(F)+eps; 
            % F=z1(F(1:end))';
            % mini(1)=min(F); maxi(1)=max(F);
            F_tab(:,i)=F;
        else if i==2
                mask2 = edge(BW2);    
                Im.MeanPretty(~mask2) = NaN;
                F1 = weighted'*deltaF_data(:,1:end);
                F_tab(:,i)=F1;
            else
                mask3=edge(BW3);    
                Im.MeanPretty(~mask3) = NaN;
                F2 = weighted'*deltaF_data(:,1:end);
                F2=detrend(F2);
                F_tab(:,i)=F2;
            end
        end
    end

    %% plot ROI

    V.fast_plot     = 0;
    V.fast_iter_max = 5;
    V.dt            = 0.118;
    V.T             = Im.T;
    V.est_sig       = 1;
    V.est_lam       = 1;
    V.est_gam       = 1;
    V.est_b         = 1;
    V.est_a         = 1;

    % figure,
    plotFigures(fig_num, numROI, F_tab, V, prettytiff, F_TH);
end
