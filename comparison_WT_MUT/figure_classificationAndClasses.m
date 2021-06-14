
% scatter plot WT vs MUT
colors = parula(5);
C = [colors(1:3,:); 0.768 0.078 0.752; 0.494 0.494 0.494];

figure,
s               = paramsWT_MUT(:,[1,2,end]);
cluster_flag    = 1;
slope_flag      = 0;
[a2, slope] = scatterPlot_sigs(s, C, cluster_flag, slope_flag);

%% percentage cell line and ISI e.g. 30% of WT neurons have 1 < ISI < 2
% look at TH+ neurons

figure,
plot_pies(sigClass, xlabels, C(1:4,:))

%%
% parameters for subplots
numRows = 2;
numCol  = length(sigClass);

figure,

% scatter classified signals
s = sigClassAll(:,[1,2,end-1]);
cluster_flag    = 1;                % scatter plot using clusters or not
slope_flag      = 0;
ax1 = subplot(numRows,numCol,[1,numRows]);
[a1, ~] = scatterPlot_sigs(s, [0.85 0.325 0.098; 0 0.447 0.741],  cluster_flag, slope_flag);
set(a1,'parent',ax1)
xlabel(ax1, []);
legend({'WT' 'LRRK2'})

% scatter classified signals
s = sigClassAll(:,[1,2,end]);
cluster_flag    = 1;                % scatter plot using clusters or not
slope_flag      = 0;
ax1 = subplot(numRows,numCol,[numRows+1,numCol]);
[a1, ~] = scatterPlot_sigs(s, C, cluster_flag, slope_flag);
set(a1,'parent',ax1) 
xlabel(ax1, 'Class 1', 'Color', 'w', 'BackgroundColor', C(1,:));
ylabel(ax1, []);
for ii = 1 : length(sigClass)
    leg{ii} = ['Class ' num2str(ii)];
end
legend(leg)

%% scatter single classes
% figure,
for ii = 1 : length(sigClass)

    classNum = ii;
    s = sigClass{classNum}(:,[1,2,end]);
    cluster_flag    = 1;                % scatter plot using clusters or not
    slope_flag      = 0;
    ax1 = subplot(numRows,numCol,ii+numCol);
    [a2, slope] = scatterPlot_sigs(s, [0.85 0.325 0.098; 0 0.447 0.741], cluster_flag, slope_flag);
    set(a2,'parent',ax1)
    legend(ax1, 'off');
    if ii ~= 1
        ylabel(ax1, []);
        xlabel(ax1, ['Class ' num2str(ii)], 'Color', 'w', 'BackgroundColor', C(ii,:));
    end

end
