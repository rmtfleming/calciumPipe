function paramsCaSA = get_mean_std_params(signals, nb_peak_thr)

% This function calculates the mean and std of the signal parameters
% INPUT: 
%     signals: containing the parameters
% OUTPU:
%     paramsCaSA: data + mean and std of each parameter
%     

fileIdx1              = [];
sigID                 = [];
isi.data              = [];
amp.data              = [];
SWidth.data           = [];
SArea.data            = [];
timeToPeak.data       = [];
releasingRate.data    = [];
removingRate.data     = [];

for i = 1 : length(signals)
    fileIdx1             = [fileIdx1; repmat(i,length(signals(i).isi),1)];
    sigID               = [sigID; signals(i).sigID];
    isi.data            = [isi.data; signals(i).isi];
    amp.data            = [amp.data; signals(i).amp];
    SWidth.data         = [SWidth.data; signals(i).SWidth];
    SArea.data          = [SArea.data; signals(i).SArea];
    timeToPeak.data     = [timeToPeak.data; signals(i).timeToPeak];
    releasingRate.data  = [releasingRate.data; signals(i).releasingRate];
    removingRate.data   = [removingRate.data; signals(i).removingRate];
    
end

% separate signal IDs
[n, sigID2] = find_nb_occurences(sigID(:,end));

% remove signals with less than 10 peaks
sig_10 = find(n < nb_peak_thr);
kk = [];
for i = 1 : length(sig_10)
    % find signals idx to remove from sigID
    a = sum(n(1 : sig_10(i)-1)) +1 : sum(n(1:sig_10(i)-1)) + n(sig_10(i));
    kk = [kk, a];
end

n(sig_10)       = [];
sigID2(sig_10)  = [];

fileIdx1(kk)            = [];
sigID(kk)               = [];
isi.data(kk)            = [];
amp.data(kk)            = [];
SWidth.data(kk)         = [];
SArea.data(kk)          = [];
timeToPeak.data(kk)     = [];
releasingRate.data(kk)  = [];
removingRate.data(kk)   = [];

% get ISI for single signals and calculate the mean and std for each
for i = 1 : length(n)
    fileIdx(i) = fileIdx1(sum(n(1:i-1))+1);
    
    isi.mean(i) = mean(isi.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    isi.std(i) = std(isi.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    
    amp.mean(i) = mean(amp.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    amp.std(i) = std(amp.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    
    SWidth.mean(i) = mean(SWidth.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    SWidth.std(i) = std(SWidth.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    
    SArea.mean(i) = mean(SArea.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    SArea.std(i) = std(SArea.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    
    timeToPeak.mean(i) = mean(timeToPeak.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    timeToPeak.std(i) = std(timeToPeak.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    
    releasingRate.mean(i) = mean(releasingRate.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    releasingRate.std(i) = std(releasingRate.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    
    removingRate.mean(i) = mean(removingRate.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));
    removingRate.std(i) = std(removingRate.data(sum(n(1:i-1))+1 : sum(n(1:i-1))+n(i)));

end

% remove outliers
X = isi.mean;
Y = isi.std;
potential_outlier = removeOutliers(X, Y, 0, 4);


fileIdx(potential_outlier) = [];

sigID2(potential_outlier) = [];

isi.mean(potential_outlier) = [];
isi.std(potential_outlier) = [];

amp.mean(potential_outlier) = [];
amp.std(potential_outlier) = [];

SWidth.mean(potential_outlier) = [];
SWidth.std(potential_outlier) = [];

SArea.mean(potential_outlier) = [];
SArea.std(potential_outlier) = [];

timeToPeak.mean(potential_outlier) = [];
timeToPeak.std(potential_outlier) = [];

releasingRate.mean(potential_outlier) = [];
releasingRate.std(potential_outlier) = [];

removingRate.mean(potential_outlier) = [];
removingRate.std(potential_outlier) = [];

% paramsCaSA
paramsCaSA.sigFileId      = [fileIdx', sigID2'];
paramsCaSA.isi            = isi;
paramsCaSA.amp            = amp;
paramsCaSA.SWidth         = SWidth;
paramsCaSA.SArea          = SArea;
paramsCaSA.timeToPeak     = timeToPeak;
paramsCaSA.releasingRate  = releasingRate;
paramsCaSA.removingRate   = removingRate;
