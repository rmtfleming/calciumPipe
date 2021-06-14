function [paramsWT_MUT] = classify_sigs(paramsWT_MUT, isiWT, isiMUT, t_interval, clustMeth)

% This function classifies the signals using kmeans or ISI intervals
% INPUTS:
%     paramsWT_MUT: matrix containing all the estimated parameters
%     isiWT       : mean and std of ISI for WT
%     isiMUT      : mean and std of ISI for MUT
%     clustMeth   : 'kmeans' or intervals
%     
% OUTPUT:
%     paramsWT_MUT : parameter matrix with the clusters

if strcmp(clustMeth, 'kmeans')
    
    numClusters = 8;
    s = [[isiWT.mean' isiWT.std']; [isiMUT.mean' isiMUT.std']];
    clustIdx = kmeans(s, numClusters);
    a = size(paramsWT_MUT);
    paramsWT_MUT(:, a(end)+1) = clustIdx;
    
else
    % sort vector according to ISI
    paramsWT_MUT = sortrows(paramsWT_MUT);
    paramsWT_MUT(:,end+1)=0;
    
    j = 1;
    for i = paramsWT_MUT(1:1) : t_interval : round(max(paramsWT_MUT(:,1)))
        idx = find(paramsWT_MUT(:,1)>=i & paramsWT_MUT(:,1)<i+t_interval);
        paramsWT_MUT(length(find(paramsWT_MUT(:,end)~=0))+1: ...
            length(find(paramsWT_MUT(:,end)~=0))+length(idx),end) = j;
        j = j+1;
        
    end
end