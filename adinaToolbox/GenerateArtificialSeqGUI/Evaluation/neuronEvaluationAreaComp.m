function [neuronEval,neuronAssembliesEval,SPSNA_neuron_spatial,SPSNA_assemblies_spatial] = ...
    neuronEvaluationAreaComp(Ucells_GT,X,U_cells_data,D,thr_area)


% Ucells_GT:    matrix of dimension #frames x #cells which contains the
%               activation of each cell along the time axis (frames)
% X:            GT dictionary. 3D array of dimensions sizePicture x sizePicture x numCells

% U_cells_data: matrix of dimension numCells x #frames which contains the
%               predicted activation of each cell along the time axis (frames)
% D:            predicted dictionary. 2D matrix of dimensions sizePicture^2 x numCells

%thr_area = 0.75; % threshold of the similarity measure

neuronAssembliesEval = zeros(size(Ucells_GT,1),4); 
neuronEval = zeros(size(Ucells_GT,1),4);

maxGT = max(max(max(X(:,:,:))));

for i = 1:size(Ucells_GT,1),
    % find the cluster(s) which is (are) active
    idxNeurons = find(Ucells_GT(i,:)>0);
    % idxNeurons expresses the idx of the GT cells involved
    if ~isempty(idxNeurons)
        X_tmp = X(:,:,idxNeurons)>0.15*maxGT; %%%%%%%%%%%%%%%%% depends on the sensibility of the extract algorithm
        idxNeuron_data  = find(U_cells_data(i,:)>0);
        % idxNeuron_data expresses the idx of the DATA cells involved
        if ~isempty(idxNeuron_data)
            numCellsInvolved_data = length(idxNeuron_data);
            evalMat = zeros(numCellsInvolved_data,length(idxNeurons));
            countLabels = 0;
            for q = 1:numCellsInvolved_data
                imD = imopen(reshape(D(:,idxNeuron_data(q))>0,[size(X,1) size(X,2)]),ones(4));
                [imLabel,numLabels] = bwlabel(imD);
                countLabels = countLabels+numLabels;
                for k = 1:numLabels
                    imLbk = (imLabel==k);
                    for r = 1:length(idxNeurons)
                        intersection = X_tmp(:,:,r).*imLbk;
                        area_intersection = length(find(intersection>0)); %sum(intersection(:)>0));
                        union = X_tmp(:,:,r)+imLbk;
                        area_union = length(find(union>0));
                        measure = area_intersection/area_union;
                        if measure>=thr_area
                            evalMat(q,r) = evalMat(q,r)+1;
                        end
                    end
                end
            end
                
            
            resultAssemblie = sum(evalMat,2); % the values for each data cell: [a;b;...;n]
            resultCells = sum(evalMat,1); % the values for each GT cell: [I II ... N]
            
            zeroElem = find(resultCells==0);
            % every 0 in the vector resultCells means a False Negative FN
            if ~isempty(zeroElem)
                neuronEval(i,4) = neuronEval(i,4)+length(zeroElem);
            end
            
            % every 0 in the vector resultAssemblie means a False Positive
            % (assuming that there are no empty dictionaries).
            zeroElem = find(resultAssemblie==0);
            if ~isempty(zeroElem)
                neuronEval(i,3) = neuronEval(i,3)+length(zeroElem);
            end
            
            % the difference between the number of labels and the real
            % number of cells (GT) are considered as FP.
            if countLabels-length(idxNeurons)>0
                neuronEval(i,3) = neuronEval(i,3)+(countLabels-length(idxNeurons));
            end
            
            % we always work with the shortest vector. However, the FNs are
            % always computed over resultCells
            if size(evalMat,1) >= size(evalMat,2)
                vectorSum = resultCells;
            else
                vectorSum = resultAssemblie;
            end
            
            for zz = 1:(length(vectorSum))
                if vectorSum(zz) == 1
                    neuronEval(i,1) = neuronEval(i,1)+1;
                elseif vectorSum(zz) > 1
                    neuronEval(i,3) = neuronEval(i,3)+abs(vectorSum(zz))-1;
                    neuronEval(i,1) = neuronEval(i,1)+1;
                end
            end
        else
            % if idxNeurons is not empty and idxNeuron_data is empty, it means that
            % the system has ignored "length(idxNeurons)" cells, so this algorithm have
            % to report "length(idxNeurons)" False Negatives (FN):
            neuronEval(i,4) = neuronEval(i,4)+length(idxNeurons);
        end
        %check for true positives, false negatives, and true negatives
    else
        idxNeuron_data  = find(U_cells_data(i,:)>0);
        if ~isempty(idxNeuron_data),
            for q = 1:length(idxNeuron_data),
                imD = imerode(reshape(D(:,idxNeuron_data(q))>0,[size(X,1) size(X,2)]),ones(4));
                [~,numLabels] = bwlabel(imD);
                neuronEval(i,3) = neuronEval(i,3)+numLabels;
            end
        else
            % if the 2 arrays are empty, we have a True Negative TN
            neuronEval(i,2) = neuronEval(i,2)+1;
        end
    end
    if neuronEval(i,3)~=0 || neuronEval(i,4)~=0
        if neuronEval(i,3) > neuronEval(i,4)
            neuronAssembliesEval(i,3) = 1;
        else
            neuronAssembliesEval(i,4) = 1;
        end
    elseif neuronEval(i,1)~=0
        neuronAssembliesEval(i,1) = 1;
    elseif neuronEval(i,2)~=0
        neuronAssembliesEval(i,2) = 1;
    end
end

SPSNA_assemblies_spatial = zeros(1,5);
[SPSNA_assemblies_spatial(1),SPSNA_assemblies_spatial(2),SPSNA_assemblies_spatial(3), ...
    SPSNA_assemblies_spatial(4),SPSNA_assemblies_spatial(5)] = extractValues(neuronAssembliesEval);

SPSNA_neuron_spatial = zeros(1,5);
[SPSNA_neuron_spatial(1),SPSNA_neuron_spatial(2),SPSNA_neuron_spatial(3),SPSNA_neuron_spatial(4),SPSNA_neuron_spatial(5)] = ... 
    extractValues(neuronEval);
    