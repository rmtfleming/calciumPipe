function n_bestt = spikeTrainInference_oopsi(V, F_tab)

% This function infers spike trains from fluorescence traces using fast
% oopsi algorithm
% INPUTS:
%     F: fluorescence traces of size timepoints x numberROI
%     V: parameters for fast_oopsi
%     
% OUTPUTS:
%     n_bestt: inferred spike trains for each fluorescence trace


for i = 1:size(F_tab,2)
    
    F = F_tab(:,i);
    
    V.F = F;
    [n_best(:,i) Phat V C] = fast_oopsi(V.F,V);    
    n_best(1,i) = 0;
        
    [row col] = find(n_best(:,i)<(max(n_best(2:end,i))+min(n_best(2:end,i)))/1.4);
    n_bestt(:,i) = n_best(:,i);
    n_bestt(row,i) = 0;

    [row col] = find(n_bestt(:,i)~=0);
    n_bestt(row,i) = 1.2;
    
end

end