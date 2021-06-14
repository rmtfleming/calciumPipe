function [neuronTempEval,neuronTempEvalParam,SPSNA_temporal] = neuronTemporalEvaluation(reassignment,Ucells_data,Ucells_GT,temporalPattern)


Ucells_data_reord = zeros(size(Ucells_GT,1),size(Ucells_GT,2));
sum_cols = zeros(1,size(Ucells_GT,2));
for i = 1:length(reassignment)
    if reassignment(i)~=0
        Ucells_data_reord(:,i) = Ucells_data(:,reassignment(i));
    end
    sum_cols(i) = sum(Ucells_data_reord(:,i));
end

active_tracks = find(sum_cols>0);
neuronTempEvalParam.tempCorresp = (length(active_tracks)/size(Ucells_GT,2))*100;
% percentage of correspondence between the GT temporal tracks and the
% "useful" detected by the algorithm.

neuronTempEvalParam.detectedTracks = size(Ucells_data,2);
neuronTempEvalParam.correctTracks = length(active_tracks);
neuronTempEvalParam.rejectedTracks = size(Ucells_data,2)-length(active_tracks);

% SHIFTS matrix:
% the 1st column contains the difference between the init. of the GT and
% the init. of the data
% the 2nd column contains the difference between the maximum locations
% the 3rd column contains the difference between the end locations
% the 4th column contains the t_rise of each predicted activation

neuronTempEval = false(size(Ucells_GT,1),4,size(Ucells_GT,2));
labelsGT = 0;
labelsData = 0;
correspLabel = 0;
idx = 1;

for i = 1:size(Ucells_GT,2)
    if sum_cols(i)~=0
        [imLabel_data,numLabels_data] = bwlabel(Ucells_data_reord(:,i));
        labelsData = labelsData+numLabels_data;
        [imLabel_GT,numLabels_GT] = bwlabel(Ucells_GT(:,i));
        labelsGT = labelsGT+numLabels_GT;
        for j = 1:numLabels_GT
            imLbk = (imLabel_GT==j);
            corr = logical(imLbk.*imLabel_data);
            if sum(corr)~=0
                correspLabel = correspLabel+1;
                numLab = imLabel_data(find(corr,1,'first'));
                imLbkdata = (imLabel_data==numLab);
                init_data = find(imLbkdata,1,'first');
                end_data = find(imLbkdata,1,'last');
                if isequal(temporalPattern,'slope')
                    max_data = init_data;
                else
                    max_data = find(imregionalmax(Ucells_data_reord((init_data:end_data),i))==1);
                    max_data = max_data(1)+init_data-1;
                end
                init_GT = find(imLbk,1,'first');
                end_GT = find(imLbk,1,'last');
                max_GT = find(imregionalmax(Ucells_GT((init_GT:end_GT),i))==1);
                max_GT = max_GT(1)+init_GT-1;
                shifts(idx,1) = init_data-init_GT;
                shifts(idx,2) = max_data-max_GT;
                shifts(idx,3) = end_data-end_GT;
                shifts(idx,4) = max_data-init_data;
                idx = idx+1;
            end
        end
        for j = 1:size(Ucells_GT,1)
            if Ucells_GT(j,i)>0 && Ucells_data_reord(j,i)>0
                neuronTempEval(j,1,i) = 1;
            elseif Ucells_GT(j,i)==0 && Ucells_data_reord(j,i)==0
                neuronTempEval(j,2,i) = 1;
            elseif Ucells_GT(j,i)==0 && Ucells_data_reord(j,i)>0
                neuronTempEval(j,3,i) = 1;
            elseif Ucells_GT(j,i)>0 && Ucells_data_reord(j,i)==0
                neuronTempEval(j,4,i) = 1;
            end
        end
    else
        neuronTempEval(:,:,i) = -1;
    end
end

if idx==1
    shifts = [NaN NaN NaN NaN];
end

% sensibility shows the percentage of predicted excitations w.r.t. the GT
neuronTempEvalParam.excitationSensibility = (labelsData/labelsGT)*100;

% accuracy shows the percentage of corresponding predicted excitations
% w.r.t. the total number of predicted excitations
neuronTempEvalParam.excitationAccuracy = (correspLabel/labelsData)*100;      


means = mean(shifts,1);
stds = std(shifts,0,1);


neuronTempEvalParam.meanShiftInit = means(1);
neuronTempEvalParam.meanShiftMax = means(2);
neuronTempEvalParam.meanShiftEnd = means(3);
neuronTempEvalParam.meanT_rise = means(4);
neuronTempEvalParam.stdShiftInit = stds(1);
neuronTempEvalParam.stdShiftMax = stds(2);
neuronTempEvalParam.stdShiftEnd = stds(3);
neuronTempEvalParam.stdT_rise = stds(4);

SPSNA_temporal = zeros(size(Ucells_GT,2),5);
% SPSNA_temporal contains in each row the Sensitivity, Precision, Specificity,
% Negative and Accuracy temporal values for each cell

for i = 1:size(Ucells_GT,2)
    if sum_cols(i)==0
        SPSNA_temporal(i,:) = -1;
    else
        [SPSNA_temporal(i,1),SPSNA_temporal(i,2),SPSNA_temporal(i,3),SPSNA_temporal(i,4),SPSNA_temporal(i,5)] = ... 
            extractValues(neuronTempEval(:,:,i));
    end
end

                
