function dist = KLDivergence_random_iterations(data1, data2, rand_percent, num_iter, parameters, saveFlag, pathKLDivPercent)

% This function calculates the Kullback-Leibler divergence between two distributions 
% then runs a certain number of itrations where the labels are randomised
% and the KLD is calculated. Then, the distribution of the divergences is
% plotted 

% INPUTS:
%     data1       : n x m matrix where each column reprensents the values
%                   of a certain parameter for one of the populations (e.g. controls)
%     data2       : same as data1 for the population to compare
%     rand_percent: the percentage of data to randomise (e.g. 0.2)
%     num_iter    : number of ransomisations
%     parameters  : a cell array of strings containing the name of each
%                   parameter
%     
% OUTPUTS:
%     dist    : an array of KL divergence at each iteration


data1(:,end+1) = zeros(1,size(data1,1));
data2(:,end+1) = ones(1,size(data2,1));
data = [data1; data2];

data_shuff = data;

% randomise a percentage
dist = [];
for jj = 1 : size(data,2)-1

    [fWT,xiWT] = ksdensity(data_shuff(1:length(find(data_shuff(:,end)==0)),jj));
    [fMUT,xiMUT] = ksdensity(data_shuff(length(find(data_shuff(:,end)==0))+1:end,jj));
    dist(1) = KLDiv(fMUT, fWT);

    % randomly shuffle
    for ii = 2 : num_iter
        % get random idexes of rows to radomise
        idx_to_shuff = randperm(size(data_shuff,1),...
            round(size(data_shuff,1)*rand_percent));

        % randomise the labels
        shuff_percent = data_shuff(idx_to_shuff,end);
        shuff_percent = shuff_percent(randperm(size(shuff_percent,1)));
        data_shuff(idx_to_shuff,end) = shuff_percent;

        data_WT = data_shuff(find(data_shuff(:,end)==0),:);
        data_MUT = data_shuff(find(data_shuff(:,end)==1),:);

        [fWT,xiWT] = ksdensity(data_WT(:,jj));
        [fMUT,xiMUT] = ksdensity(data_MUT(:,jj));

        dist(ii) = KLDiv(fMUT, fWT);
    end

    figure, plot(dist,'LineWidth', 2), hold on, plot(1, dist(1), 'r*', 'Markersize', 12,'LineWidth', 2)
        if find(dist == max(dist)) ~= 1
            hold on, plot(find(dist == max(dist)), max(dist), 'g*', 'Markersize', 12,'LineWidth', 2)
        end
    set(gca,'FontSize', 20,'FontWeight','bold')
    xlabel('Iteration'), xlim([-5 ii])
    ylabel('Distance')
    title(parameters{(jj)})
    grid on
    
    if saveFlag
        figname=[pathKLDivPercent parameters{(jj)} '_iterationPlot'];
    %     print('-dtiff',figname);
    end
end