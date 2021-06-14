% '~/code/Imaging/analysis/fast-oopsi/scripts/';

clear all, clc, close all
path        = '/mnt/siham_hachi/codeInGit/Imaging/analysis/fast-oopsi/';
datadir     = [path 'data/'];
fileName1   = '3lane_TH_positive';
fileName2   = '3lane_TH_negative';
figdir     = [path 'figs/'];
avg_tif   = 'AVG_TH.tif';

load([datadir fileName1]);

Pl.xlims=[1 V.T];
Pl.nticks=4;
Pl.fs=12;
tvec(1,1)=0;
for i=2:V.T
    tvec(1,i)=V.dt*i;
end
V.F=z1(V.F);

fig=figure(1); clf,
    
%plot stack_TH
subplot('position',[0.12 0.5 0.3 0.4]); hold all
avg_tiff=imread([datadir avg_tif]);
avg_tiff=imresize(avg_tiff,3);
imagesc(flipud(avg_tiff));
colormap('hot');
title('mean frame','FontSize',Pl.fs);
axis off;
set(gca,'YDir','normal')

%fast oopsi
[n_best Phat V C] = fast_oopsi(V.F,V);
n_best(1,1)=0;

    
[raw col]=find(n_best<(max(n_best)+min(n_best))/1.5);
n_bestt=n_best;
n_bestt(raw,col)=0;
[raw col] = find(n_bestt~=0);
n_bestt(raw,col)=1;


%plot fluorescence trace TH and spikes
subplot('position',[0.1 0.1 0.35 0.3]); hold all
plot(tvec,0.4*(Phat.a\V.F)+1, 'b');
axis([0 V.dt*V.T 0 3.2]);
set(gca,'ytick',[]);
xlabel('time (sec)','FontSize',Pl.fs);
ylabel('Spikes    F (a.u.)','FontSize',Pl.fs);
%     bar(tvec,z);
bar(tvec,n_best);
bar(tvec,n_bestt, 'EdgeColor','r');

%%
load([datadir fileName2]);
F_tab=V.F_tab;
F_tab=z1(F_tab');

wid=0.35;
hei=0.3;
left=0.55;
bottom=0.1;
c='b';
cutoffs=[0.55 0.4 0.5];

subplot('position',[left bottom wid hei]); hold on
axis([0 V.dt*V.T 0 3.2]);
set(gca,'ytick',[]);
xlabel('time (sec)','FontSize',Pl.fs);
ylabel('Spikes     F (a.u.)','FontSize',Pl.fs);

j=1; yy=0;
colorr=[{'k'},{'k'},{'k'}];
text_tab=[{'soma'},{'neurite'},{'neurite'}];
for i=1:3
    F=F_tab(:,i);
    
    %fast oopsi
    V.F=F;
    [n_best Phat V C] = fast_oopsi(V.F,V);
    n_best(1,1)=0;
    plot(tvec, 0.5*(Phat.a\F)+1, c); hold on

    [raw col]=find(n_best<(max(n_best)+min(n_best))/1.5);
    n_bestt=n_best;
    n_bestt(raw,col)=0;
    [raw col] = find(n_bestt~=0);
    n_bestt(raw,col)=0.8;
    
    bar(tvec,n_best);
    bar(tvec,n_bestt,'EdgeColor','r');

    bottom=bottom+0.32;
    subplot('position',[left bottom wid hei]); hold on
    axis([0 V.dt*V.T 0 3.2]);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    ylabel('Spikes     F (a.u.)','FontSize',Pl.fs);
    
    wh=[5 4];   %width and height
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
end

figname=[figdir 'fig5'];
print('-dtiff',figname);