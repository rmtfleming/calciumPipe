
% plot 2 spikes

% load('wellG6_3_8bits_time_traces.mat');

halfwindow      = 100;
quantile        = 0.2;
step            = 5;
sigma           = 3;
i               = 16;

[dF_F Fbase] = deltaF_F(time_traces(16,:), halfwindow, quantile, step, sigma);

tvec=[];
T = size(time_traces,2);
for ii = 2:T+1
    tvec(1,ii-1)=(1/5)*ii;
end

dF_F = z1(dF_F);
plot(tvec, dF_F), hold on
% plot([12; 12], [-0.1; 0], '-k',  [12; 15], [-0.1; -0.1], '-k', 'LineWidth', 2)
% t1 = text(12, -0.15, '5 sec','FontSize',12, 'Color', 'k', 'FontWeight','bold');
% t2 = text(10, -0.12, '20% \DeltaF/F','FontSize',12, 'Color', 'k', 'FontWeight','bold');
% set(t2,'Rotation',90); 
% axis off

[~, n] = findpeaks(dF_F,'MinPeakHeight', ...
    mean(dF_F), 'MinPeakDistance',5);

plot(tvec(n), dF_F(n), 'r.');