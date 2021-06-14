function Z = fftFiltering(time_traces, m)

X = time_traces'; 
N = length(X);
k = repmat(ifftshift(-floor(N/2):floor((N-1)/2)),size(time_traces,1),1);
Y = fft(X); %//perform fft
Y = Y.*(abs(k')<=m); %// zero out all frequencies larger than 'm'
Z = ifft(Y,'symmetric');

end

% X = time_traces(9,:); N = length(X);
% k = ifftshift(-floor(N/2):floor((N-1)/2)); %//compute the frequency vector
% Y = fft(X); %//perform fft
% Y = Y.*(abs(k)<=50); %// zero out all frequencies larger than 'm'
% Z = ifft(Y,'symmetric');
