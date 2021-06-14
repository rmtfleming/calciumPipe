function proba_distribution(WTmean, MUTmean, parameters, saveFlag, pathDistrib)

    [fWT,xiWT] = ksdensity(WTmean,'support','positive');
    [fMUT,xiMUT] = ksdensity(MUTmean,'support','positive');
    figure,
    h1 = plot(xiWT,fWT,'Color', [0.85 0.325 0.098]);
    hold on
    h2 = plot(xiMUT,fMUT, 'Color', [0 0.44 0.74]);
    hold off
    alpha(.8)
    set(h1,'LineWidth',3)
    set(h2,'LineWidth',3)
    set(gca,'FontSize', 30,'FontWeight','bold', 'linewidth',3)

    legend([h1 h2],' WT',' LRRK2');
    xlabel([parameters ' (s)'])
    ylabel('Probability')
    box off

    if saveFlag
        figname=[pathDistrib parameters];
        print('-dtiff',figname);
    end
