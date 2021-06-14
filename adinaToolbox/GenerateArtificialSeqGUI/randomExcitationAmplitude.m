function Ucells_GT = randomExcitationAmplitude(Ucells_GT,range)

precision = 10000;
for i = 1:size(Ucells_GT,2)
    [imLabel,numLabels] = bwlabel(Ucells_GT(:,i));
    for j = 1:numLabels
        activation = find(imLabel==j,1,'first');
        finish = find(imLabel==j,1,'last');
        excitation_label = Ucells_GT((activation:finish),i);
        randAmplitude = randi([min(range)*precision max(range)*precision],1)/precision;
        excitation_label = randAmplitude*excitation_label/max(excitation_label);
        Ucells_GT((activation:finish),i) = excitation_label;
    end
end