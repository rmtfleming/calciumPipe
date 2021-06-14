%**************************************************************************
%     Copyright (C) Oct. 24, 2008  By YU Weimiao @ BII Singapore
%
%     This file is part of NeuronCyto Version 1.0.
% 
%     NeuronCyto is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     NeuronCyto is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with NeuronCyto.  If not, see <http://www.gnu.org/licenses/>.
%
%**************************************************************************
function im= function_ImPreProcessing(im, SizeOfTopHat, CellORNucleus)
% [im]= ImPreProcessing(im, SizeOfTopHat, SizeOfFilter).
% Pre-process of the image.
%   The input arguments:
%           im: Image for pre-process
%           SizeOfTopHat: Size of the TopHar operation.(remove non-uniform 
%                               background)
%           SizeOfFilter: Neiborhood of filter size. (medfile2, wiener2)
%           Diffusion: Apply LPA-ICI approach to improve the qulity of the
%                          image, also return the diffusion image.
%   The output arguments:
%           Image: the image loaded according to given path and file name
%------------------------------------------------------------------------
%   Class Support
%   -------------
%     TInput image should be uint8. Output image is the same class.
%-------------------------------------------------------------------------
% By YU Weimiao @ BII Singapore
% Date: Sept. 11, 2007.

% input image is uint16
% output image is uint8
im= im2double(im);

if strcmp(CellORNucleus, 'Nucleus')
    
    % remove the noise by wavelet decomposition;
    [thr,sorh,keepapp] = ddencmp('den','wp',im);
    im = wdencmp('gbl',im,'sym4',2,thr,sorh,keepapp);
    
    %remove the non-uniform background, if any;
    im = imtophat(im,strel('disk', SizeOfTopHat));
    
    % normlise the image graylevel value to [0,1];
    im = (im - min(im(:)))./(max(im(:))-min(im(:)));
    
    % estimate the noise level(std) of the image (after remove the noise);
    dev = function_stdEst2D(im, 2);
    
    % addjust the graylevel value of image for better dynamics;
    %im = imadjust(im,[min(0.1, 2.*dev) max(0.1, 2.*dev+0.05)],[0 1]);2008_04
    im = imadjust(im,[0.05 0.1],[0 1]);
    
    % remove any artifacts introduce by the above processing
    im = imtophat(im,strel('disk', SizeOfTopHat));
    
    % normlise the image graylevel value to [0, 255];
    im = uint8(255.*(im - min(im(:)))./(max(im(:))-min(im(:))));

elseif strcmp(CellORNucleus, 'Cell')
    
   % remove the noise by wavelet decomposition;
    [thr,sorh,keepapp] = ddencmp('den','wv',im);
    im1 = wdencmp('gbl',im,'sym4',2,thr,sorh,keepapp);
    
    %remove the non-uniform background, if any;
    im1 = imtophat(im1,strel('disk', SizeOfTopHat));
    im1=(im1-min(im1(:)))./(max(im1(:))-min(im1(:)));

    % estimate the noise level(std) of the image (after remove the noise);
    dev = function_stdEst2D(im1 ,2);
    
    % addjust the graylevel value of image for better dynamics;
    im1 = imadjust(im1, [ min(0.05, 3.*dev)  1], [0 1]); 

    %im1 = histeq(im1);

    % normlise the image graylevel value to [0, 255];
    im =uint8(255.*(im1-min(im1(:)))./(max(im1(:))-min(im1(:))));
else
    error('input/output arguments of ImPreProcessing are not correct\n');
    error('Please use "help ImPreProcessing"\n');
end





