function denoised_dict = denoiseGaussPow(dict,gaussdim,power)
% DESCRIPTION:
%     denoise method based on a gaussian power law
%
% USAGE:   
%       denoised_dict = denoiseGaussPow(dict,gaussdim,power)
%
% ARGUMENTS:
%         dict:  input image
%         gaussiandim: size of the gaussian filter
%         power: order of the power law
%       
% OUTPUT: 
%         denoised_dict:  denoised image
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".
if nargin == 1
    gaussdim = [round(size(dict,1)/50) round(size(dict,2)/50)];
    power = 30;
    
elseif nargin == 2
    power = 30;
    
elseif nargin > 3 || nargin < 1
    error('Incorrect number of parameters');
    
end

denoised_dict = zeros(size(dict));
sigma = gaussdim(1)/3;
h = fspecial('gaussian',gaussdim,sigma);

for i = 1:size(dict,3)
    dict(:,:,i) = dict(:,:,i)/max(max(dict(:,:,i)));
    
    %%%
    tmp = reshape(dict(:,:,i),[],1);
    mean_image = mean(tmp);
    %%%
    
    noiseLev = noiseLevelExtraction(dict(:,:,i));
    
    %%%
    noiseLev = 2*noiseLev + mean_image;     %%% sigma multiplicity
    %%%
    
    %im_filtered = medfilt2(dict(:,:,i),gaussdim);
    im_filtered = imfilter(dict(:,:,i),h);
    im_filtered = im_filtered/(max(max(im_filtered)));
    impow = zeros(size(dict,1),size(dict,2));
    for k=1:size(im_filtered,1)
        for j = 1:size(im_filtered,2)
            impow(k,j) = (power^im_filtered(k,j)-1)/(power-1);
        end
    end
    impow(impow<(noiseLev)) = 0;
    impow = imopen(impow,ones(4));
    
    denoised_dict(:,:,i) = impow;
end

