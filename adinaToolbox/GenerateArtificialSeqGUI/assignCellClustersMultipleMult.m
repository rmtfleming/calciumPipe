function [flagNeurons, assemblies] = assignCellClustersMultipleMult(numGroups,X,A,percentMult,maxMult)

assemblies = cell(1,numGroups);
nDicts = size(X,3);
numMult = round(percentMult*nDicts);
flagNeurons = zeros(nDicts,maxMult+1);
flagNeurons(:,maxMult+1) = 1;

% here we compute the multiplicities of the groups. The last column of
% flagNeurons indicates the multiplicity of that cell
multiplicity = zeros(numMult,2);
p = randperm(nDicts); % we randomly select the numMult candidates
p = p(1:numMult);
multiplicity(:,1) = p(:);
multiplicity(:,2) = randi([1 maxMult],[numMult 1]);

for i = 1:numMult
    flagNeurons(multiplicity(i,1),maxMult+1) = multiplicity(i,2);
end

for i = 1:nDicts
    idx = find(A(i,:)>0);
    for j = 1:flagNeurons(i,maxMult+1) % number of multiplicities of cell "multiplicity(i,1)"
        currentGroup = randi([1 numGroups],1);
        cellsInThatGroup = assemblies{currentGroup};
        % setdiff returns the values of idx that are not in cellsInThatGroup.
        % So if length of aux is smaller than length of idx it means that there
        % is some repeated values between idx and cellsInThatGroup, ergo that
        % group is not assignable to that cell.
        aux = setdiff(idx,cellsInThatGroup);
        while length(aux)<length(idx)
            currentGroup = randi([1 numGroups],1);
            cellsInThatGroup = assemblies{currentGroup};
            aux = setdiff(idx,cellsInThatGroup);
        end
        cellsInThatGroup = horzcat(cellsInThatGroup, i);
        assemblies{currentGroup} = cellsInThatGroup;
        flagNeurons(i,j) = currentGroup;
    end
end

       
flagNeurons(:,maxMult+1) = [];

% if there is any empty cluster, we redistribute the cells of the largest
% cluster into it. The next subroutine assures us that there will be at least one
% cell per cluster.

empties = cellfun('isempty',assemblies);

if ~empties
    return
else
    idx = find(empties==1); % idx = array that contains clusters which are empty
    for i = 1:length(idx)
        % largerGroup contains the index of the largest cluster, and len
        % contains the length of this cluster
        [len, largerGroup] = max(cellfun('length',assemblies));
        index = randi([1 len],1);
        exchangeCell = assemblies{largerGroup}(index);
        assemblies{largerGroup}(index) = [];
        assemblies{idx(i)} = exchangeCell;
        index = find(flagNeurons(exchangeCell,:)==largerGroup);
        flagNeurons(exchangeCell,index) = idx(i);
    end
end