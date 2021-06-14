
clear all, clc, close all

path         = '/mnt/siham_hachi/data/fast-oopsi/3lane/';
imname       = 'full_small_stack';
tifname     = [path 'stacks/' imname '.tif'];
figdir      = [path 'figs/'];
datadir     = [path 'mat/'];
fname       = 'full_stack_TH_noisy_image_cNeurons';
LoadTif     = 0;

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

%% Create (F-F_0)/F0 data

if ~isfield(Im, 'deltaF_data')
    for i=1:Im.T
        img=reshape(Im.DataMat(:,i),Im.h,Im.w);
        data(:,:,i)=img(:,:);
    end
    param = initial_parameter_preprocessing();
    param.method = 3;                                   % choose (F-F_0)/F0 method
    [seq_param] = extractingNormalizedSeq(data,param);
    for i=1:Im.T
        x = seq_param.data_denoise(:,:,i);
        Im.deltaF_data(:,i) = x(:);
    end
    Im.deltaF_data = Im.deltaF_data;
    imWavelet = waveletDenoisingDraguhn3DScale(seq_param.data_denoise(:,:,1),[4 4 3]);
    Im.imWavelet = imWavelet;
    save([datadir fname],'Im','-v7.3')
end

%% get fluo signals from ROIs

numROIs = length(Im.cNeurons);
F_tab = [];

for i = 1:numROIs
    
    roi = Im.cNeurons(i).imMask;
    
    a = ones(size(Im.MeanFrame));
    weighted = a.*roi(:);
    
    F = weighted'*Im.deltaF_data(:,1:end);
    F = detrend(F);
    F_tab(:,i) = F;
    
end


%% plot ROI
numROIs = 8;
V.fast_plot     = 0;
V.fast_iter_max = 5;
V.dt            = 2;%0.118;
V.T             = Im.T;
V.est_sig       = 1;
V.est_lam       = 1;
V.est_gam       = 1;
V.est_b         = 1;
V.est_a         = 1;

tvec(1,1) = 0;
for i = 2:V.T+1
    tvec(1,i-1)=V.dt*i;
end

Pl.xlims=[1 V.T];
Pl.nticks=4;
Pl = PlotParams(Pl);
Pl.fs=12;

wid = 0.8;
hei = 0.8;
left = 0.1;
bottom = 0.1;
for i=1:numROIs
    F_tab(:,i) = z1(F_tab(:,i)');
end

% a=reshape(Im.DataMat(:,1),Im.h,Im.w);
% displaySingleCells(Im.imWavelet,Im.cNeurons);
% figure,
% imagesc(Im.imWavelet), axis off, axis square, axis equal, colormap 'gray'; %hold on
% for j = 1:length(Im.cNeurons),
% %     hold on,plot(Im.cNeurons(j).obj.mu(2),Im.cNeurons(j).obj.mu(1),'.k','MarkerSize',7)
%     hold on, contour(Im.cNeurons(j).imMask,1,'Color','r','LineWidth',1);
%     text(Im.cNeurons(j).obj.mu(2)+2,Im.cNeurons(j).obj.mu(1)+2,num2str(j),'FontSize',12,'Color','w');
% end

figure,
subplot('position',[left bottom wid hei]); hold on
axis([0 V.dt*V.T 0 numROIs]);
set(gca,'ytick',[]);
xlabel('time (sec)','FontSize',Pl.fs);

for i = 1:numROIs
    
    F = F_tab(:,i);
    
    V.F = F;
    [n_best(:,i) Phat V C] = fast_oopsi(V.F,V);    
    n_best(1,i) = 0;
    
    plot(tvec, F + i-1); hold on
    
    [row col] = find(n_best(:,i)<(max(n_best(2:end,i))+min(n_best(2:end,i)))/1.4);
    n_bestt(:,i) = n_best(:,i);
    n_bestt(row,i) = 0;

    [row col] = find(n_bestt(:,i)~=0);
    n_bestt(row,i) = 1.2;
    
end

wh=[4 4];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');

hei = 0.8/numROIs;
figure,
h2 = subplot('position',[left bottom wid hei]); hold on
axis([0 V.dt*V.T 0 1.5]);
set(gca,'ytick',[]);
xlabel('time (sec)','FontSize',Pl.fs);

for i = 1:numROIs
    
%     bar(tvec,z1(n_best(:,i)'));
    bar(tvec,n_bestt(:,i),'EdgeColor','r');
    
    if i ~= numROIs
        bottom = bottom + (0.8/numROIs);
        subplot('position',[left bottom wid hei]); hold on
        axis([0 V.dt*V.T 0 1.5]);
        set(gca,'ytick',[]);
        set(gca,'xtick',[]);
    end    
end

wh=[4 4];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
%     figname=[figdir 'fig7'];

