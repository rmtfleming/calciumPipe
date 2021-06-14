%cd '~/code/Imaging/analysis/fast-oopsi/scripts/';

clear all, clc, close all
path        = '/mnt/siham_hachi/codeInGit/Imaging/analysis/fast-oopsi/';
datadir     = [path 'data/'];
figdir     = [path 'figs/'];
fileName   = '2lane';

load([datadir fileName]);

Pl.xlims=[1 V.T];
Pl.nticks=4;
Pl.fs=12;
tvec(1,1)=0;
for i=2:V.T
    tvec(1,i)=V.dt*i;
end
F_tab=V.F_tab;

cd '../functions/';

fig=figure(2); clf,

wid=0.8;
hei=0.3;
left=0.1;
bottom=0.1;
c='b';
cutoffs=[0.2 0.45 0.5];

subplot('position',[left bottom wid hei]); hold on
axis([0 V.dt*V.T 0 3.2]);
set(gca,'ytick',[]);
xlabel('time (sec)','FontSize',Pl.fs);
ylabel('Spikes     F (a.u.)','FontSize',Pl.fs);

j=1; yy=0;
colorr=[{'k'},{'k'},{'k'}];
text_tab=[{'neurite'},{'neurite'},{'soma'}];
for i=1:3
    F=F_tab(:,i);
    
    %fast oopsi
    V.F=F;
    [n_best Phat V C] = fast_oopsi(V.F,V);
    n_best(1,1)=0;
    plot(tvec, 0.5*(Phat.a\F)+1.2, c); hold on
    text(120,1.7,text_tab(i),'color',char(colorr(i)));
    if i==1
        text(119,1.4,'terminal','color',char(colorr(i)));
    end

    [raw col]=find(n_best<(max(n_best)+min(n_best))/1.5);
    n_bestt=n_best;
    n_bestt(raw,col)=0;

    [raw col] = find(n_bestt~=0);
    n_bestt(raw,col)=1; 
    
    bar(tvec,n_best);
    bar(tvec,n_bestt,'EdgeColor','r');

    bottom=bottom+0.32;
    subplot('position',[left bottom wid hei]); hold on
    axis([0 V.dt*V.T 0 3.2]);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    ylabel('Spikes     F (a.u.)','FontSize',Pl.fs);
end

wh=[4 4];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
figname=[figdir 'fig6'];
print('-dtiff',figname);