function reassignment = findCorrespondence(dict_GT,dict_data,thr_area)

% inputs:
% 
% dict_GT: matrix of nrows x nrows x nDicts. It contains the GT dictionaries.
% dict_data: matrix of nrows^2 x nPredDicts or nrows x nrows x nPredDicts. It contains the predicted dicts.
% thr_area: threshold of matching between areas
% 
% 
% output:
% 
% reassignment: vector which indicates the relation between data and GT

if size(dict_data,3) == 1
    dict_data = reshape(dict_data,size(dict_GT,1),size(dict_GT,2),[]);
end

corresp_mat = zeros(size(dict_data,3),size(dict_GT,3));

maxGT = max(max(max(dict_GT(:,:,:))));
for i = 1:size(dict_GT,3)
    dict_GT(:,:,i) = dict_GT(:,:,i)>0.05*maxGT;
end

for i = 1:size(dict_data,3)
    imD = imopen(dict_data(:,:,i),ones(4));
    [imLabel,numLabels] = bwlabel(imD);
    if numLabels>1 || numLabels==0
        corresp_mat(i,:) = 0;
    else
        for j = 1:size(dict_GT,3)
            intersection = dict_GT(:,:,j).*imLabel;
            area_intersection = sum(sum(intersection>0));
            union = dict_GT(:,:,j)+imLabel;
            area_union = sum(sum(union>0));
            measure = area_intersection/area_union;
            if measure>=thr_area
                corresp_mat(i,j) = 1;
            end
        end
    end
end

% if a row sums more than 1, we have to put the whole row to 0. That is for
% preventing cases in which some pred. dictionaries has big areas such that
% covers more than one GT_dict cell (and the user has put a low area threshold).
sum_rows = sum(corresp_mat,2);
err = find(sum_rows>1);
for i = 1:length(err)
    corresp_mat(err(i),:) = 0;
end

reassignment = zeros(1,size(dict_GT,3));

for i = 1:size(dict_GT,3)
    pos = find(corresp_mat(:,i)==1);
    if ~isempty(pos)
        reassignment(i) = pos(1); % here we take the first pred. dict if there is more than one predicted dicts for that GT_dict
    end
end