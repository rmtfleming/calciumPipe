function filtSig = gaussianFilt(sig, sizeFilt, sigma)

    x = linspace(-sizeFilt/2, sizeFilt/2, sizeFilt);
    gaussFilter = exp(-x .^ 2/(2 * sigma ^ 2));
    gaussFilter = gaussFilter/sum (gaussFilter); % normalize
%     filtSig = filter (gaussFilter,1, sig);
    filtSig = conv2(sig, gaussFilter, 'same');
end