function [Ucells_GT] = assignActivationPatternToCells_v2(type,tau,t_rise,frameRate,clipping,Ucells_temporalPattern,activationPattern)

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
                lengthRise = round((t_rise(2)/1000)*frameRate); %length (in frames) of the rise period
                rand_activation = activation + 0*round(lengthRise*randn(1));
                
                while rand_activation<=0 || rand_activation>(size(Ucells_temporalPattern,1)-length(activationPattern))
                    rand_activation = activation + round(lengthRise*randn(1));
                end
                
                Uexp = zeros(1,size(Ucells_temporalPattern,1));
                Uexp(rand_activation) = 1;
                rand_tau = randi([round(tau(1)*100) round(tau(2)*100)],1)/100;
                rand_t_rise = randi([round(t_rise(1)) round(t_rise(2))],1);
                
                if sum(imLbk)<=length(activationPattern)
                    rand_activationPattern = generateActivationPattern(type,rand_tau,rand_t_rise,frameRate,clipping);
                    rand_activationPattern = rand_activationPattern*normFactor/max(rand_activationPattern);%%%%%%%%
                    Uexp = conv(Uexp,rand_activationPattern);
                    len = length(rand_activationPattern);
                    init = find(Uexp>0,1,'first');
                    Ucells_GT((rand_activation:(rand_activation+len-1)),i) = Uexp(init:(init+len-1));
                else
                    lengthRise = round((rand_t_rise/1000)*frameRate); %length (in frames) of the rise period
                    timeFramesDecay = (sum(imLbk)+1)-(lengthRise-1); % we remove a
                    % sample from lengthRise because when constructing the
                    % activation pattern, we eliminate a sample of the expRise.
                    % We also add a sample on sum(imLbk) because the imLbk has
                    % been estimated taking into account that the final sample
                    % of the activation pattern has been eliminated.
                    timeSecondsDecay = timeFramesDecay/frameRate;
                    newTau = timeSecondsDecay/3;
                    span_tau = tau(2)-tau(1);
                    rand_newTau = randi([round((newTau-span_tau)*100) round(newTau*100)],1)/100;
                    newActPattern = generateActivationPattern(type,rand_newTau,rand_t_rise,frameRate,clipping);
                    
                    newActPattern = newActPattern*normFactor/max(newActPattern);
                    newActPattern(newActPattern < (max(newActPattern)*exp(-3))) = [];
                    
                    Uexp = conv(Uexp,newActPattern);
                    len = length(newActPattern);
                    init = find(Uexp>0,1,'first');
                    Ucells_GT((rand_activation:(rand_activation+len-1)),i) = Uexp(init:(init+len-1));
                end
            end
        end
end
