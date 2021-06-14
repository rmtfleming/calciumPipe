function plot_comparison_bars(sigClass, paramsLabels, xlabels, pathBarPlots)

    for i = 1 : length(sigClass)
        clust = sigClass{i};
        % consider only big clusters
        if size(clust,1) > 0
            for j = 1 : length(paramsLabels)
                % plot parameter
                clust_WT = clust(find(clust(:,end)==1),j*2-1);
                clust_MUT = clust(find(clust(:,end)==2),j*2-1);
                figure
                bar(1:2,[mean(clust_WT), mean(clust_MUT)], 0.5, 'FaceColor', [.2 .5 .7])
                hold on, errorbar(1:2,[mean(clust_WT), mean(clust_MUT)], ...
                    [std(clust_WT), std(clust_MUT)],'.', 'LineWidth', 3, 'Color', 'k')
                set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold')
                title(gca, ['Cluster' num2str(i)]);
                ylabel(paramsLabels{j})
                grid on
                figname=[pathBarPlots regexprep(paramsLabels{j},' ','_') num2str(i)];
                print('-dtiff',figname);
            end
            close all;
        end
    end