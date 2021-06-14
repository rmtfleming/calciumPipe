function badSigIdx = removeNoisySignals(dF_F_all_z1, sizeGauss, sigma)

% This function removes signals that represent just noise by comparing the 
% number of peaks detected in raw fluorescence signals and the ones detected 
% in the Gaussian smoothed signals
%     INPUTS:
%         dF_F_all_z1 : matrix of fluorescence signals
%         sizeGauss   : the size of the gaussian filter
%         sigma       : sigma for the gaussian filter
%             
%     OUTPUTS:
%         badSigIdx   : indexes of the niosy signals

    dF_F_all_z1_gauss = gaussianFilt(dF_F_all_z1, sizeGauss, sigma);
    % calculate the energy and range fo each signal
    energy = sum(dF_F_all_z1'.*dF_F_all_z1');
    energyFilt = find(energy < mean(energy));
    ranges = range(dF_F_all_z1');
    rangesFilt = find(ranges<mean(ranges));
    % 
    badSigIdx = intersect(energyFilt,rangesFilt);


    % [up,lo] = envelope(dF_F_all',10,'peak');
    % areaBtweenEnv = sum(up-lo);
    % areaBtweenEnvFilt = find(areaBtweenEnv<(std(areaBtweenEnv)/3));

    % Remove signal that represent just noise by comparing the number of
    % detected peaks in the raw signal and the Gaussian smoothed signal
    for j = 1 : length(badSigIdx)

        i = badSigIdx(j);
        % Find peaks in the dF/F signals
        [~, n] = findpeaks(dF_F_all_z1(i,:));
        n_filt = dF_F_all_z1(i,n)< 2*std(dF_F_all_z1(i,n));
        n(n_filt)=[];
        spks = zeros(1,400);
        spks(n)=1;

        % Find peaks in the Gaussian filtered signals
        [~, n] = findpeaks(dF_F_all_z1_gauss(i,:));
        n_filt = dF_F_all_z1_gauss(i,n)<2*std(dF_F_all_z1_gauss(i,n));
        n(n_filt)=[];
        spksGauss = zeros(1,400);
        spksGauss(n)=1;

        % ratio_peaks(j) = length(nonzeros(spks))/length(nonzeros(spksGauss));

        % remove the signals
        if length(nonzeros(spks))/length(nonzeros(spksGauss))<=3
            badSigIdx(j)=0;
        end

    end
    
    badSigIdx(badSigIdx==0)=[];
end
