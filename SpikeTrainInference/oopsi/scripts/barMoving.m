load('F.mat');
wh=[2.66 1];

for i=1:Im.T
    x=tvec(1,i);
    figure,
    plot(tvec, F, 'k');hold on,
    plot([x x],[0 1],'r','LineWidth',1);
    axis([0 V.dt*Im.T+25 0 1]);
    set(gcf,'PaperSize',wh,'PaperPosition',[0 0 wh], 'color','w');
    plot([125 135],[0.04 0.04],'k','LineWidth',1);
    text(126,0,'10 s');
    plot([135 135],[0.04 0.24],'k','LineWidth',1);
    text(136,0.15,'20\% $\Delta$F/F','interpreter','latex','FontSize',5);
    axis 'off';
    imname=['/mnt/siham_hachi/data/fast-oopsi/traces_for_video/' num2str(i) '.tif'];
%     imname=['P:\data\fast-oopsi\traces_for_video\' num2str(i) '.tif'];
    print('-dtiff',imname);
    close;
end