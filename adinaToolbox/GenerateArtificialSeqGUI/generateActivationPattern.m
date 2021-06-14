function [activationPattern] = generateActivationPattern(type,tau_maxmin,t_rise_maxmin,frameRate,clipping)

% Given the type and the temporal parameters, the function returns the
% activation pattern

tau = max(tau_maxmin);
t_rise = max(t_rise_maxmin);

lengthAssemSec = 10*tau; % lenght (in seconds) of the decay period
lengthAssemFr = round(lengthAssemSec*frameRate); %length in frames of the decay period
lengthRise = round((t_rise/1000)*frameRate); %length (in frames) of the rise period


switch(type)
    
    case 'exp_exp'
        time_rise = 0:1/frameRate:(lengthRise-1)/frameRate;
        time_decay = 0:1/frameRate:(lengthAssemFr-1)/frameRate;
        
        exponentialFunc = exp(-time_decay/tau);
        exponentialRise = exp(-time_rise/(t_rise/3000)); % t_rise = 3*tau_rise (in milliseconds)
        
        activationPattern = imfilter([1-exponentialRise(1:end-1),(1-min(exponentialRise))*exponentialFunc],[1 1 1]/3);
        activationPattern = activationPattern(1:end-1)/max(activationPattern); 
        % here we normalize and refuse the last sample, because is
        % corrupted due to the filtering.
        
        % here the samples which are below the clipping percentage w.r.t.
        % the maxima are processed. The first values below the threshold
        % are eliminated, while the lasts are matched with the last value
        % over the threshold
        activationPattern(activationPattern < (max(activationPattern)*clipping)) = 0;
        init = find(activationPattern>0,1,'first');
        if init~=1
            for i = 1:(init-1)
                activationPattern(1) = [];
            end
        end
        fin = find(activationPattern>0,1,'last');
        for i = (fin+1):length(activationPattern)
            activationPattern(i) = activationPattern(fin);
        end
        
        activationPattern(activationPattern < (max(activationPattern)*exp(-3))) = [];
        
    case 'gauss_exp'
        gaussFunc = fspecial('gaussian', [lengthRise*2 1], lengthRise/3);
        
        time = 0:1/frameRate:(lengthAssemFr-1)/frameRate;
        
        exponentialFunc = exp(-time/tau);
        activationPattern = conv(exponentialFunc, gaussFunc);
        
        activationPattern(activationPattern < (max(activationPattern)*clipping)) = 0;
        init = find(activationPattern>0,1,'first');
        if init~=1
            for i = 1:(init-1)
                activationPattern(1) = [];
            end
        end
        fin = find(activationPattern>0,1,'last');
        for i = (fin+1):length(activationPattern)
            activationPattern(i) = activationPattern(fin);
        end
        
        activationPattern(activationPattern < (max(activationPattern)*exp(-3))) = [];
        
    case 'slope'
        activationPattern = ones(1,(lengthAssemFr+lengthRise));
        
    otherwise
        msgbox('Unknown activation pattern','Error','error','modal')
        error('Unknown activation pattern');
        
end

