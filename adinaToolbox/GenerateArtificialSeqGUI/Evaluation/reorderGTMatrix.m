function [Ucells_GT_reord] = reorderGTMatrix(dict_GT,dict_data,Ucells_GT,Ucells_data,thr_area)

reassignment = findCorrespondence(dict_GT,dict_data,thr_area);
num_data_dicts = size(dict_data);

reassignment_inv = zeros(1,num_data_dicts(end));
for i = 1:num_data_dicts(end)
    pos = find(reassignment==i);
    if ~isempty(pos)
        reassignment_inv(i) = pos;
    end
end

Ucells_GT_reord = zeros(size(Ucells_data,1),size(Ucells_data,2));
for i = 1:length(reassignment_inv)
    if reassignment_inv(i)~=0
        Ucells_GT_reord(:,i) = Ucells_GT(:,reassignment_inv(i));
    end
end