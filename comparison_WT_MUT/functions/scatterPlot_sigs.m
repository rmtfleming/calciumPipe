function [a, slope] = scatterPlot_sigs(s, colors, cluster_flag, slope_flag)

% This function scatter plots the mean of ISI against std of ISI of each signal
% INPUTS:
%     s             : a matrix of mean and std of ISI for each signal
%     cluster_flag  : flag for showing clusters or not
%     slope_flag    : flag for estimating slope
%
% OUTPUTS:
%     slope     : slope of the line passing through the points

% figure,

if cluster_flag == 1
    idx = s(:,end);
    a = gscatter(s(:,1),s(:,2),idx, colors);
    set(gca,'linewidth',2)
%     set(gca, 'xlim', [floor(min(s(:,1))) floor(max(s(:,1)))+1])
else
    a = gscatter(s(:,1),s(:,2));
end

set(gca,'FontSize', 14,'FontWeight','bold')
xlabel('mean ISI (s)')
ylabel('std ISI (s)')
box off
% grid on

slope = 0;
if slope_flag
    % regression line
    format long
    b1 = s(:,1)\s(:,2);
    yCalc1 = b1*s(:,1);
    hold on
    plot(s(:,1),yCalc1, 'b', 'LineWidth',2)

    slope = (yCalc1(2) - yCalc1(1)) ./ (s(2,1) - s(1,1));
    hold on
    text(min(s(:,1)),max(s(:,2)), ['slope = ' num2str(slope)], 'FontSize', 14, 'FontWeight', 'bold')
end
