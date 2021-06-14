function [Ucells_GT] = assignActivationPatternToCells(type,t_rise,frameRate,clipping,Ucells_temporalPattern,activationPattern)

switch(type)
    
    case 'slope'
        Ucells_GT = Ucells_temporalPattern;
        
    otherwise
        Ucells_GT = zeros(size(Ucells_temporalPattern,1),size(Ucells_temporalPattern,2));
        normFactor = max(activationPattern);
        for i = 1:size(Ucells_temporalPattern,2)
            [imLabel,numLabels] = bwlabel(Ucells_temporalPattern(:,i));
            for j = 1:numLabels
                imLbk = (imLabel==j);
                activation = find(imLbk==1,1,'first');
                Uexp = zeros(1,size(Ucells_temporalPattern,1));
                Uexp(activation) = 1;
                if sum(imLbk)<=length(activationPattern)
                    Uexp = conv(Uexp,activationPattern);%%%%%%%%%
                    len = sum(imLbk);
                    init = find(Uexp>0,1,'first');
                    Ucells_GT((activation:(activation+len-1)),i) = Uexp(init:(init+len-1));
                else
                    lengthRise = round((t_rise/1000)*frameRate); %length (in frames) of the rise period
                    timeFramesDecay = (sum(imLbk)+1)-(lengthRise-1); % we remove a
                    % sample from lengthRise because when constructing the
                    % activation pattern, we eliminate a sample of the expRise.
                    % We also add a sample on sum(imLbk) because the imLbk has
                    % been estimated taking into account that the final sample
                    % of the activation pattern has been eliminated.
                    timeSecondsDecay = timeFramesDecay/frameRate;
                    newTau = timeSecondsDecay/3;
                    
                    newActPattern = generateActivationPattern(type,newTau,t_rise,frameRate,clipping);
                    
                    newActPattern = newActPattern*normFactor/max(newActPattern);
                    newActPattern(newActPattern < (max(newActPattern)*exp(-3))) = [];
                    
                    Uexp = conv(Uexp,newActPattern);
                    len = sum(imLbk);
                    init = find(Uexp>0,1,'first');
                    Ucells_GT((activation:(activation+len-1)),i) = Uexp(init:(init+len-1));
                end
            end
        end
end
