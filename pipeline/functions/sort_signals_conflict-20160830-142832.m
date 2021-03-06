function [F_tab_sorted, spikes_sorted, firing_rate_sorted, cNeurons_sorted, th_percent] = ...
    sort_signals(F_tab, spike_times, cNeurons,  method)

for i = 1 : size(F_tab,1)    

    F = F_tab(i,:);
    
    % spike detection
    [m, n] = findpeaks(F,'MinPeakHeight', mean(F), 'MinPeakDistance',5);

    spike_times(i,n) = 1;
    nSpikes(i) = length(n);
    
    % firing rate estimation
    [firing_rates(i,:),t] = FiringRate(nonzeros(n),fps,Im.T);    
end

th_percent = 0;
if strcmp(method,'freq')
    
    % order the signals by number of events
    [nb_peaks, sig_index] = sort(nSpikes);

    for i = 1 : length(sig_index)
        F_tab_sorted(i,:) = F_tab(sig_index(i),:);
        spikes_sorted(i,:) = spike_times(sig_index(i),:);
        firing_rate_sorted(i,:) = firing_rates(sig_index(i),:);
        cNeurons_sorted(i) = cNeurons(sig_index(i));
    end
else if strcmp(method,'th')
        for i=1:length(cNeurons)
            if ~isempty(cNeurons(i).th_pos)
                th_cN(i) = cNeurons(i).th_pos;
            else
                th_cN(i)=0;
            end
        end
        [th_intensity, neuron_index] = sort(-th_cN);
        for i = 1 : length(neuron_index)
            F_tab_sorted(i,:) = F_tab(neuron_index(i),:);
            spikes_sorted(i,:) = spike_times(neuron_index(i),:);
            firing_rate_sorted(i,:) = firing_rates(neuron_index(i),:);
            cNeurons_sorted(i) = cNeurons(neuron_index(i));
        end
        th_intensity = -th_intensity;
        th_percent = ((th_intensity - min(th_intensity))*100)./(max(th_intensity)-min(th_intensity));
    end
end