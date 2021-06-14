clear all, clc,

original= imread('/mnt/siham_hachi/MATLAB codes/fast-oopsi-master/data/2lane/small_stack2_raw.tif','Index',2);

se = strel('disk',100);
tophatFiltered = imtophat(original,se);

background = imopen(original,strel('disk',100));
im=original-background;

im_adapthisteq = adapthisteq(im,'Range','full');

med_im=medfilt2(im_adapthisteq);

figure, imshow(original); figure, imshow(tophatFiltered); figure, imshow(im);
figure, imshow(im_adapthisteq); figure, imshow(med_im);

%%
% FileTif='../../../data/2lane/stack1_8bits.tif';
% InfoImage=imfinfo(FileTif);
% NumberImages=length(InfoImage);
% for j=1:NumberImages
%     img8 = imread(FileTif,'Index',j);
%     image1(:,:,j)=img8;
% end
% 
% FileTif='../../../data/2lane/stack2_8bits.tif';
% InfoImage=imfinfo(FileTif);
% NumberImages=length(InfoImage);
% for j=1:NumberImages
%     img8 = imread(FileTif,'Index',j);
%     image2(:,:,j)=img8;
% end

% FileTif3='/mnt/siham_hachi/MATLAB codes/fast-oopsi-master/data/New folder/stack3.tif';
% InfoImage=imfinfo(FileTif);
% NumberImages=length(InfoImage);
% for j=1:NumberImages
%     image3(:,:,j) = imread(FileTif,'Index',j);
% end

% FinalImage=cat(3,image1, image2);
%%
% clear all, clc,
% image=imread('/mnt/siham_hachi/MATLAB codes/fast-oopsi-master/data/3lane/stack_small_TH.tif','Index',1);
% im_median=medfilt2(image);
% H = fspecial('gaussian',[3 3],2);
% im_gaussian=imfilter(image,H,'replicate');
% figure, imshow(image);
% figure, imshow(im_median);
% figure, imshow(im_gaussian);

%%
cd '/mnt/siham_hachi/code/Imaging/analysis/fast-oopsi/functions/';

% Pl.xlims=[1 Im.T];
% Pl.nticks=4;
% Pl = PlotParams(Pl);
% Pl.fs=12;
V.fast_plot = 0;
V.fast_iter_max=5;
V.dt=0.118;
V.T=Im.T;
V.est_sig=1;
V.est_lam=    1 ;
V.est_gam=    1 ;
V.est_b=    1 ;
V.est_a=      1 ;
tvec(1,1)=0;
for i=2:Im.T
    tvec(1,i)=V.dt*i;
end
% prettytiff = imread(prettytif);
% 
% fig=figure(2); clf,
% height=0.4;
% width=Im.h/Im.w*height;
% left = (1-width)/2;
% bottom=(.5-height)/2+.5;
% 
% %plot stack_TH
% subplot('position',[0.12 0.5 0.3 0.4]); hold all
% prettytiff=imresize(prettytiff,4);
% imagesc(flipud(prettytiff));
% colormap('hot');
% title('mean frame','FontSize',Pl.fs);
% axis off;
% set(gca,'YDir','normal')
% 
% %fast oopsi
% F_TH=load([datadir 'F_TH_noisy_image.mat']);
% V.F=F_TH.F;
% [n_best Phat V C] = fast_oopsi(V.F,V);
% [raw col] = find(n_best<0.7);
% n_best(raw,col)=0;
% z=n_best(1:V.T);
% [raw1 col1] = find(z>0.0001);
% z(raw1,col1) = 0.8;
% 
% %plot fluo trace TH and spikes
% subplot('position',[0.1 0.1 0.35 0.3]); hold all
% plot(tvec,0.4*(Phat.a\F_TH.F)+1, 'r');
% axis([0 V.dt*Im.T 0 3]);
% set(gca,'ytick',[]);
% xlabel('time (sec)','FontSize',Pl.fs);
% ylabel('Spikes    F (a.u.)','FontSize',Pl.fs);
% bar(tvec,z);

if plott== 1
    wid=0.8; %0.35;
    hei=0.3; %0.17;
    left=0.1; %0.55;
    bottom=0.1;
    F_tab=z1(F_tab');    
    
%     figure,
%     subplot('position',[left bottom wid hei]); hold on
%     plot(tvec, F_tab(:,1), 'r'); hold on
%     plot(tvec, F_tab(:,2), 'b'); hold on
%     plot(tvec, F_tab(:,3), 'k'); hold on
%     set(gca,'ytick',[]);
%     xlabel('time (sec)','FontSize',Pl.fs);
%     ylabel('F','FontSize',Pl.fs);
%     title('Fluorescence observations','FontSize',Pl.fs);
    
    hei=0.27; %0.2;
%     bottom=bottom+0.3; %0.25;
    subplot('position',[left bottom wid hei]); hold on
    axis([0 V.dt*Im.T 0 4]);
    set(gca,'ytick',[]);
    xlabel('time (sec)','FontSize',Pl.fs);
    ylabel('Spikes     F (a.u.)','FontSize',Pl.fs);

    j=1; yy=0;
    colorr=[{'k'},{'k'},{'k'}];
    cutoffs=[0.2 0.45 0.5];
    text_tab=[{'soma'},{'dendrite'},{'dendritic'}];
    for i=1:numROI
        F=F_tab(:,i);
        c=char(colorr(i));
        %fast oopsi
        V.F=F;
        [n_best Phat V C] = fast_oopsi(V.F,V);
        n_best2;
        plot(tvec, 0.5*(Phat.a\F)+1, 'r'); hold on
        text(120,1.7,text_tab(i),'color',char(colorr(i)));
        if i==3
            text(118,1.5,'terminal','color',char(colorr(i)));
        end
        [raw col] = find(n_best<cutoffs(i));
        n_best(raw,col)=0;
        z=n_best(1:V.T);
        [raw1 col1] = find(z>0.0001);
        z(raw1,col1) = 0.8;
        bar(tvec,z);
        bottom=bottom+0.32;
        subplot('position',[left bottom wid hei]); hold on
        axis([0 V.dt*Im.T 0 3]);
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        ylabel('Spikes     F (a.u.)','FontSize',Pl.fs);
        %     percent=(raw1*100)/V.T;
        %     x1_z=(percent*V.dt)*10;
        %     for ii=1:length(x1_z)
        %         plot([x1_z(ii) x1_z(ii)], [y1_z y1_z+0.8],'k'); hold on;
        %     end
    end
    
%     axis([0 V.dt*Im.T 0 2]);
%     set(gca,'ytick',[]);
%     set(gca,'xtick',[]);
%     xlabel('time (sec)','FontSize',Pl.fs);
%     ylabel('Spikes    F (a.u.)','FontSize',Pl.fs);

    wh=[4 4];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
    figname=[figdir 'fig7']; 
%     print('-dpng',figname);
%     print('-dtiff',figname);
%     print('-depsc',figname);
%     print('-dpdf',figname);
%     saveas(fig,figname);
end
%%
V.F=F;
tvec(1,1)=0;
for i=2:V.T
tvec(1,i)=V.dt*i;
end
[n_best Phat V C] = fast_oopsi(V.F,V);
n_best(1,1)=0;

[raw col]=find(n_best<(max(n_best)+min(n_best))/1.3);
n_bestt=n_best;
n_bestt(raw,col)=0;
figure,
plot(tvec, 0.5*(Phat.a\F)+1, 'k'); hold on
bar(tvec,n_bestt);

PARMHAT = expfit(n_best);
[raw col]=find(n_best<(PARMHAT*5));
n_bestt=n_best;
n_bestt(raw,col)=0; figure,
plot(tvec, 0.5*(Phat.a\F)+1, 'k'); hold on
bar(tvec,n_bestt);