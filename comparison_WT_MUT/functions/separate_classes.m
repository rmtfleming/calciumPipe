function sigClass = separate_classes(paramsWT_MUT, numSigThr)

% sort the cluster/class column
paramsWT_MUT = sortrows(paramsWT_MUT,size(paramsWT_MUT,2));

% get the number of signals in each cluster
[nSig_per_Clust, ~] = find_nb_occurences(paramsWT_MUT(:,end));

% remove clusters/classes with less than 10 signals
nSig_per_Clust(find(nSig_per_Clust < numSigThr)) = [];

% separate the clusters
for i = 1 : length(nSig_per_Clust)
    sigClass{i} = paramsWT_MUT(sum(nSig_per_Clust(1:i-1))+1 : ... 
        sum(nSig_per_Clust(1:i-1))+nSig_per_Clust(i),1:end-1);
    sigClass{i} = sortrows(sigClass{i},size(sigClass{i},2));
    
    % Remove the clusters with number of WT or MUT signals < 10 
    if (length(find(sigClass{i}(:,end)==1) < numSigThr) || ...
            length(find(sigClass{i}(:,end)==2)) < numSigThr)
        sigClass_rm = i;
    end
end
sigClass{sigClass_rm} = [];
sigClass = sigClass(~cellfun('isempty',sigClass));

% plot percentage of K7 and K7M signals in each cluster
txt = {'K7 ','K7M '}; % strings
for i = 1 : length(sigClass)
    if size(sigClass{i},1) > 0
        X = [length(find(sigClass{i}(:,end)==1)) length(find(sigClass{i}(:,end)==2))];
        ax1 = subplot(floor(length(sigClass)/3)+mod(length(sigClass),3),3,i);
        p = pie(ax1,X);       
        title(ax1,['ISI : ' num2str(round((min(sigClass{i}(:,1)))*100)/100) ' - '...
            num2str(round((max(sigClass{i}(:,1)))*100)/100)], 'FontSize', 20);
        
        hText = findobj(p,'Type','text');
        percentValues = get(hText,'String');
        combinedtxt = strcat(txt,percentValues');
        
        hText(1).String = combinedtxt(1); hText(2).String = combinedtxt(2);        
        hText(1).Position = 0.5*hText(1).Position; hText(2).Position = 0.5*hText(2).Position;        
        hText(1).Color = 'w'; hText(2).Color = 'k';       
        hText(1).FontSize = 16; hText(2).FontSize = 16;
        
        hText(1).FontWeight = 'bold'; hText(2).FontWeight = 'bold';
    end
end
colormap summer