function Ucells_GT = fuseMaximums(Ucells_GT,frameRate,t_rise,type,clipping)

% This function checks if there is any label with more than one local
% maxima. If it happens, the two (or more) maximas are fused into a single
% excitation.

if isequal(type,'slope')
    return
end

for i = 1:size(Ucells_GT,2)
    [imLabel,numLabels] = bwlabel(Ucells_GT(:,i));
    for j = 1:numLabels
        activation = find(imLabel==j,1,'first');
        finish = find(imLabel==j,1,'last');
        excitation_label = Ucells_GT((activation:finish),i);
        maxima = imregionalmax(excitation_label);
        if sum(maxima) > 1
            duration = finish-activation; %in frames
            rand_t_rise = randi([round(t_rise(1)) round(t_rise(2))],1);
            lengthRise = round((rand_t_rise/1000)*frameRate); %in frames
            decay_duration = duration-lengthRise; %in frames
            tau = decay_duration/(3*frameRate);
            activationPattern = generateActivationPattern(type,tau,rand_t_rise,frameRate,clipping);
            len = length(activationPattern);
            Ucells_GT((activation:activation+len-1),i) = activationPattern;
        end
    end
end