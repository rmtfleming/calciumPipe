function plotFigures(fig_num, numROI, F_tab, V, prettytiff, F_TH)

if nargin<5
    prettytiff = [];
    F_TH=[];
end

Pl.xlims=[1 V.T];
Pl.nticks=4;
Pl = PlotParams(Pl);
Pl.fs=12;

tvec(1,1)=0;
for i=2:V.T
    tvec(1,i)=V.dt*i;
end

if fig_num==1
    
    fig=figure(2); clf,
    
    %plot stack_TH
    subplot('position',[0.12 0.5 0.3 0.4]); hold all
    prettytiff=imresize(prettytiff,4);
    imagesc(flipud(prettytiff));
    colormap('hot');
    title('mean frame','FontSize',Pl.fs);
    axis off;
    set(gca,'YDir','normal')
    
    %fast oopsi
    V.F=F_TH.F;
    [n_best Phat V C] = fast_oopsi(V.F,V);
    n_best(1,1)=0;
    
    PARMHAT = expfit(n_best);
    [raw col]=find(n_best<(PARMHAT*5));
    n_bestt=n_best;
    n_bestt(raw,col)=0;
%     
%     [raw col]=find(n_best<(max(n_best)+min(n_best))/1.5);
%     n_bestt=n_best;
%     n_bestt(raw,col)=0;
    [raw col] = find(n_bestt~=0);
    n_bestt(raw,col)=1;
    
%     [raw col] = find(n_best<0.7);
%     n_best(raw,col)=0;
%     z=n_best(1:V.T);
%     [raw1 col1] = find(z>0.0001);
%     z(raw1,col1) = 0.8;
        
    %plot fluo trace TH and spikes
    subplot('position',[0.1 0.1 0.35 0.3]); hold all
    plot(tvec,0.4*(Phat.a\F_TH.F)+1, 'r');
    axis([0 V.dt*V.T 0 3]);
    set(gca,'ytick',[]);
    xlabel('time (sec)','FontSize',Pl.fs);
    ylabel('Spikes    F (a.u.)','FontSize',Pl.fs);
%     bar(tvec,z);
    bar(tvec,n_best);
    bar(tvec,n_bestt, 'EdgeColor','r');
    
    wid=0.35;
    hei=0.17;
    left=0.55;
    bottom=0.1;
    F_tab=z1(F_tab');
    c='k';
    cutoffs=[0.55 0.4 0.5];
else
    wid=0.8;
    hei=0.3;
    left=0.1;
    bottom=0.1;
    F_tab=z1(F_tab');
    c='b';
    cutoffs=[0.2 0.45 0.5];
end       
    
hei=0.27; %0.2;
%     bottom=bottom+0.3; %0.25;
subplot('position',[left bottom wid hei]); hold on
axis([0 V.dt*V.T 0 4]);
set(gca,'ytick',[]);
xlabel('time (sec)','FontSize',Pl.fs);
ylabel('Spikes     F (a.u.)','FontSize',Pl.fs);

j=1; yy=0;
colorr=[{'k'},{'k'},{'k'}];
text_tab=[{'soma'},{'neurite'},{'neurite'}];
for i=1:numROI
    F=F_tab(:,i);
%         c=char(colorr(i));
    %fast oopsi
    V.F=F;
    [n_best Phat V C] = fast_oopsi(V.F,V);
    n_best(1,1)=0;
    plot(tvec, 0.5*(Phat.a\F)+1, c); hold on
    if fig_num~=1
        text(120,1.7,text_tab(i),'color',char(colorr(i)));
        if i==3
            text(118,1.5,'terminal','color',char(colorr(i)));
        end
    end

%     PARMHAT = expfit(n_best);
%     [raw col]=find(n_best<(PARMHAT*5));
%     n_bestt=n_best;
%     n_bestt(raw,col)=0;

%         [raw col] = find(n_best<cutoffs(i));
%         n_best(raw,col)=0;
%         z=n_best(1:V.T);
%         [raw1 col1] = find(z>0.0001);
%         z(raw1,col1) = 0.8;
%         bar(tvec,z);

        [raw col]=find(n_best<(max(n_best)+min(n_best))/1.5);
        n_bestt=n_best;
        n_bestt(raw,col)=0;

    [raw col] = find(n_bestt~=0);
    n_bestt(raw,col)=0.8; 
    
    bar(tvec,n_best);
    bar(tvec,n_bestt,'EdgeColor','r');

    bottom=bottom+0.32;
    subplot('position',[left bottom wid hei]); hold on
    axis([0 V.dt*V.T 0 3]);
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    ylabel('Spikes     F (a.u.)','FontSize',Pl.fs);
end

wh=[4 4];   %width and height
set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh],'Color','w');
%     figname=[figdir 'fig7'];
end