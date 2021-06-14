function plot_pies(sigClass, cellLabels, C)

% find number of WT signals after removing the classes with less than 10
% signals
out = cat(1,sigClass{:});

% WT
for k = 1 : length(cellLabels)
    labels = {};
    X = [];
    for i = 1 : length(sigClass)
        X = [X length(find(sigClass{i}(:,end)==k))/length(find(out(:,end)==1))];
        labels = cat(1,labels,...
            [num2str(round((min(sigClass{i}(:,1)))*100)/100) '-'...
            num2str(round((max(sigClass{i}(:,1)))*100)/100)]);
    end
    ax1 = subplot(1,2,k);
    p = pie(ax1,X);       
    title(ax1,[], 'FontSize', 26);
    title(ax1,cellLabels{k});

    hText = findobj(p,'Type','text');
    for j = 1 : length(sigClass)    
        hText(j).Color = 'k';      
        hText(j).FontSize = 22;
        hText(j).FontWeight = 'bold';
    end
end

legend(labels,'Orientation','horizontal', 'FontSize', 20,...
    'FontWeight', 'bold', 'Position', [0.344 0.138 0.329 0.058])

% colors = parula(length(sigClass));
colormap(C)
