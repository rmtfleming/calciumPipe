function fluo_trace_spikes(F, events, C, fps)

% This function plots a fluorescence trace and its underlying spike train
% INPUTS:
%     F       : vector of fluorescence trace 
%     events  : vector of spikes
%     C       : a vector of a color in RGB for the fluorescence trace
%     fps     : frame rate
    

T = size(F,2);
tvec(1,1) = 0;
for i = 2:T+1
    tvec(1,i-1)=(1/fps)*i;
end

figure, 

plot(tvec, z1(F), 'color', C, 'LineWidth', 1.5)
hold on

plot([-0.5; 9.5], [0; 0], '-k', 'LineWidth', 2);
plot([-2; -2], [0.3; 0.5], '-k', 'LineWidth', 2);
t1 = text(-5, 0.2, '20% \DeltaF/F','FontSize',14, 'Color', 'k', 'FontWeight','bold');
t2 = text(-0.5, -0.05, '10 sec','FontSize',14, 'Color', 'k', 'FontWeight','bold');
set(t1,'Rotation',90);

for j = 1:length(events)
    lh=line([events(j)/fps events(j)/fps],[1 1.5]);
    set(lh,'linestyle','-','color',[.49 .49 .49],'linewidth',1);
    lh=line([events(j)/fps events(j)/fps],[0.3 1]);
    set(lh,'linestyle','--','color',[.49 .49 .49],'linewidth',0.5);
end

axis off