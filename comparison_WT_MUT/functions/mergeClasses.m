function sigClassAll = mergeClasses(sigClass)

% This function merges all the classes in one array

sigClassAll   = cell2mat(sigClass');
sigClassAll(:,end+1) = 0;
sigClassAll(1:length(sigClass{1}), end) = 1;
classlengths = cellfun(@length, sigClass);

for ii = 2 : length(sigClass)
    sigClassAll(sum(classlengths(1:ii-1))+1 : sum(classlengths(1:ii-1))+classlengths(ii), end) = ii;    
end