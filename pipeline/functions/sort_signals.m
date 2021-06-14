function [F_tab_sorted, spikes_sorted, firing_rate_sorted, cNeurons_sorted, th_percent] = ...
    sort_signals(F_tab, cNeurons, fps,  method, signalClass)

for i = 1 : size(F_tab,1)    

    F = F_tab(i,:);
    
    % spike detection
    [~, n] = findpeaks(F,'MinPeakHeight', mean(F), 'MinPeakDistance',5);

    spike_times(i,n) = 1;
    nSpikes(i) = length(n);
    
    % firing rate estimation
    [firing_rates(i,:),~] = FiringRate(nonzeros(n),fps,size(F_tab,2));    
end

th_percent = 0;
if strcmp(method,'freq')
    
    % order the signals by number of events
    [~, sig_index] = sort(mean(firing_rates'));
    
else if strcmp(method,'number_spikes')
        
        % order the signals by number of events
        [~, sig_index] = sort(nSpikes);
        
    else if strcmp(method,'th')
            
            for i=1:length(cNeurons)
                if ~isempty(cNeurons(i).th_pos)
                    th_cN(i) = cNeurons(i).th_pos;
                else
                    th_cN(i)=0;
                end
            end
            [th_intensity, ~] = sort(-th_cN);
            th_intensity = -th_intensity;
            th_percent = ((th_intensity - min(th_intensity))*100)./(max(th_intensity)-min(th_intensity));
        
        else if strcmp(method,'isi')
                [~, sig_index] = sort(signalClass);
                sig_index = fliplr(sig_index);                
            end
        end

    end    
end

F_tab_sorted = F_tab(sig_index,:);
spikes_sorted = spike_times(sig_index,:);
firing_rate_sorted = firing_rates(sig_index,:);
cNeurons_sorted = cNeurons(sig_index);