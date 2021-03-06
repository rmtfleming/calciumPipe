%   =======================================================================================
%   Copyright (C) 2013  Erlend Hodneland
%   Email: erlend.hodneland@biomed.uib.no 
%
%   This program is free software: you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation, either version 3 of the License, or
%   (at your option) any later version.
% 
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
% 
%   You should have received a copy of the GNU General Public License
%   along with this program.  If not, see <http://www.gnu.org/licenses/>.
%   =======================================================================================

clear all;
close all;

load ../data/nucleistain_3D.mat;
imsegm = imsegm(150:400,50:300,6);

% smoothing with edge enhancing diffusion
prm.eed.kappa = 10;
prm.eed.maxniter = 300;
imsm = cellsegm.smoothim(imsegm,'eed','prm',prm);

cellsegm.show(imsegm,1);axis off;axis image;title('Raw image');
cellsegm.show(imsm,4);axis off;axis image;title('Edge enhancing diffusion');