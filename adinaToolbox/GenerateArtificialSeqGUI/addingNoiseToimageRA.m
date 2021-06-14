function [noisyImage] = addingNoiseToimageRA(input,RA)

% noiseimage = randn(size(input));
% Pinputimage=sum(sum(input.^2));                             
% Pnoiseimage=sum(sum((noiseimage).^2));       
% K = (Pinputimage/ Pnoiseimage);                                
% noise = sqrt(K/RA)*(randn(size(input))); 

sigma = (max(input(:)) - mean(input(:)))/RA;


noisyImage = input + sigma*randn(size(input));


%version 2

% input_noisy = input + 0.1*max(input(:))*randn(size(input));
% MaxMinDifftime = max(input_noisy,[],2) - mean(input_noisy,2);
% 
% sigma_pixel = MaxMinDifftime/RA;
% noisyImage = input_noisy;
% for i = 1:size(input,1),
%     noisyImage(i,:) = input_noisy(i,:) + sigma_pixel(i)*randn(size(input_noisy(i,:)));
% end



