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

function  [SeedNumber, BW, Segment, Boundary, Seeds_x, Seeds_y] = function_NucleusSegment(cellim, im, approach)
% [SeedNumber, BW, Segment, Boundary, Seeds_x, Seeds_y] = NucleusSeg(cellim, im, approach)
% The function is developed to segment the image of Nuclues and find the seeds.
%   The input arguments:
%           cellim: Cell image;
%           im Nucleus image;
%           approach: Apporach to segment nucleus.
%                                   1. 'Normal', Thresholding
%                                   2. 'Watershed', Watershed apporach
%                                   3. ''Levelset'', Leveset approach
% Notes: * caps sensitive!
%           * ''Levelset'' is recommended
%   The output arguments:
%           See the help of VariableInitilize for the data structure of
%           Nucleus
%-----------------------------------------------------------
%   Class Support
%   -------------
%     Input image should be uint8.
%-----------------------------------------------------------
% By YU Weimiao @ BII Singapore
% Date: Sept. 11, 2007.

if nargin==2
    approach='Levelset'; % by default use levelset
end

% get the size of the nucleus image
[LX,LY] = size(im);


if strcmp(approach,'Normal')

    % convert the image to binary image according to Otus's method;
    BW = im2bw(im, graythresh(im));

    % remove the thin structure which could be caused by staining.
    BW= imopen(BW,strel('disk',3));

    % Remove the holes in the Nucleus.
    % Since we set the seeds in the center of the seeds, if there is any holes
    % in the nucleus, this will cause errors in the later processing
    BW = imfill(BW,'holes');

    %----------------------------------------------------------------------------------------------------------
    % Eliminate the Seeds without Cell. The false Seeds will be removed, if any
    [Boundary, Segment, SeedNumber] = bwboundaries(BW,'nohole');

    % Convert the cell image to binary image;
    cellBW = im2bw(cellim, 0.02);

    % remove the small regions of the cell binary image, which might be
    % cuased by the noise.
    cellBW = bwareaopen(cellBW,100);

    % fill in the holes of the cell binary image.
    cellBW = imfill(cellBW,'holes');

    % get the segments of the cell images.
    [cellBoundary, cellSegment, cellNumber] = bwboundaries(cellBW,'nohole');

    % convert the cell segments and the nucleus segments to column.
    % (don't use pixel operation, since it will be very slow.)
    cellSegment_Col =im2col(cellSegment,[1 1]);
    Segment_Col =im2col(Segment, [1 1]);

    % check which segments of the nucleus do not contain cell imformation.
    % And remove them by set them to "0"(background).
    for i =1:1:SeedNumber
        y=find (Segment_Col==i);
        if sum(cellSegment_Col(y))==0
            Segment_Col(y)=0;
        end
    end

    % convert the segments column back to image and then get the binary
    % image.
    Segment = col2im(Segment_Col,[1 1],[LX LY]);
    BW = (Segment>0);

    % Remove the holes in the Nucleus.
    BW = imfill(BW,'holes');

    %----------------------------------------------------------------------------------------------------------
    % re-calculate the sgements after remove the false seeds
    [Boundary, Segment, SeedNumber] = bwboundaries(BW,'nohole');

    Seeds_x =zeros(1,SeedNumber);
    Seeds_y =zeros(1,SeedNumber);
    % get the center of the seeds.
    for i =1:1:SeedNumber
        [x y] = find(Segment==i);
        Seeds_x(i) =round(mean(x));
        Seeds_y(i) =round(mean(y));
    end

% Watershed segmentation of the nucleus image recommended for the study of
% mitosis. This approach can separate a dividing nuleis into two. If you
% consider a dividing nuclei is still one object, you should not use this
% approach.
elseif strcmp(approach,'Watershed')

    % convert the image to binary image according to Otus's method;
    BW = im2bw(im, graythresh(im));

    % remove the thin structure which could be caused by staining.
    BW= imopen(BW,strel('disk',3));

    % Remove the holes in the Nucleus.
    % Since we set the seeds in the center of the seeds, if there is any holes
    % in the nucleus, this will cause errors in the later processing
    BW = imfill(BW,'holes');

    %----------------------------------------------------------------------------------------------------------
    % Eliminate the Seeds without Cell. The false Seeds will be removed, if any
    [Boundary, Segment, SeedNumber] = bwboundaries(BW,'nohole');

    % Convert the cell image to binary image;
    cellBW = im2bw(cellim, 0.02);

    % remove the small regions of the cell binary image, which might be
    % cuased by the noise.
    cellBW = bwareaopen(cellBW,100);

    % fill in the holes of the cell binary image.
    cellBW = imfill(cellBW,'holes');

    % get the segments of the cell images.
    [cellBoundary, cellSegment, cellNumber] = bwboundaries(cellBW,'nohole');

    % convert the cell segments and the nucleus segments to column.
    % (don't use pixel operation, since it will be very slow.)
    cellSegment_Col =im2col(cellSegment,[1 1]);
    Segment_Col =im2col(Segment, [1 1]);

    % check which segments of the nucleus do not contain cell imformation.
    % And remove them by set them to "0"(background).
    for i =1:1:SeedNumber
        y=find (Segment_Col==i);
        if sum(cellSegment_Col(y))==0
            Segment_Col(y)=0;
        end
    end

    % convert the segments column back to image and then get the binary
    % image.
    Segment = col2im(Segment_Col,[1 1],[LX LY]);
    BW = (Segment>0);

    % Remove the holes in the Nucleus.
    BW = imfill(BW,'holes');

    % compute the watershed line.
    D = bwdist(~BW,'euclidean');
    SeedsBW  = imextendedmax(D,1);

    D2 = bwdist(SeedsBW,'euclidean');
    L = watershed(D2);
    SeedsEdges = (L==0);
    SeedsEdges = bwmorph(SeedsEdges,'thicken' );
    BW = BW - (SeedsEdges&BW);

    [Boundary, Segment, SeedNumber] = bwboundaries(BW,'nohole');

    Seeds_x =zeros(1,SeedNumber);
    Seeds_y =zeros(1,SeedNumber);

    for i =1:1:SeedNumber
        [x y] = find(Segment==i);
        Seeds_x(i) =round(mean(x));
        Seeds_y(i) =round(mean(y));
    end


elseif strcmp(approach,'Levelset')

    % convert the nucleus image to double;
    im = double(im);
    im = (im-min(im(:)))./(max(im(:))-min(im(:)));
    [bg std] = function_BackGroundFinder(im);
        

    % initialise the levelset function phi;
    phi=im-1;


    ITERATIONS = 50;
    delta_t = 10;
    lambda1 = 1;
    lambda2 = 10;
    nu = 0;
    h = 1; h_sq = h^2;
    epsilon = 10;
    END =0;
    c_in=[mean2(im)];
    c_out=[mean2(im)];
    iterations =0;
    TE = [];
    
    while iterations<ITERATIONS & END ==0;

        iterations = iterations+1;

        % compute dirac(phi) * delta_t factor
        dirac_delta_t = delta_t * function_dirac( phi, epsilon );

        % calculate inside and outside curve terms
        [c1,c2, inside, outside ] = function_calcenergy( im, phi, epsilon );
        energy_term = (-nu - lambda1 .* inside + lambda2 .* outside);
        % total_energy = lambda1 .* sum(sum(inside)) + lambda2 .*sum(sum(outside));
        total_energy = lambda1 .* sum(sum((im-c1).*(im-c1))) + lambda2 .*sum(sum((im-c2).*(im-c2)));
        TE = [TE  total_energy];

        c_in=[c_in c1];
        c_out=[c_out c2];

        phi = phi + dirac_delta_t .* ( energy_term );
%         figure(102);mesh(phi);
%         figure(103);imagesc(phi>0)
%         CCC =0;
        if iterations >+2
            if (abs(c_out(iterations)-bg )<3.*std)
                END =1;
            end
        end
    end
    
%     figure(4);plot(c_in);
%     figure(5);plot(c_out);
%     figure(6);plot(TE);
    % overlay the front on the original image
    BW = (phi>0);
    %BW = bwareaopen(BW,30); % remove area which is too small or caused by noise;
    BW= imopen(BW,strel('disk', 5)); % remove the thin structure which could be caused by staining.
    %------------------------------------------------------
    % Eliminate the Seeds without Cell
    % The false Seeds will be removed, if any
    [Boundary, Segment, SeedNumber] = bwboundaries(BW,'nohole');

    cellBW = im2bw(cellim, 0.02);

    cellBW = bwareaopen(cellBW,100);
    cellBW = imfill(cellBW,'holes');
    [cellBoundary, cellSegment, cellNumber] = bwboundaries(cellBW,'nohole');
    cellSegment_Col =im2col(cellSegment,[1 1]);
    Segment_Col =im2col(Segment, [1 1]);

    for i =1:1:SeedNumber
        y=find (Segment_Col==i);
        if sum(cellSegment_Col(y))==0
            Segment_Col(y)=0;
        end
    end

    Segment = col2im(Segment_Col,[1 1],[LX LY]);
    BW = (Segment>0);
    % Remove the holes in the Nucleus. Since we set the seeds in the center
    % of the seeds, if there is any holes in the nucleus, this will cause
    % errors in the later processing
    BW = imfill(BW,'holes');
    %------------------------------------------------------
    [Boundary, Segment, SeedNumber] = bwboundaries(BW,'nohole');

    Seeds_x =zeros(1,SeedNumber);
    Seeds_y =zeros(1,SeedNumber);

    for i =1:1:SeedNumber
        [x y] = find(Segment==i);
        Seeds_x(i) =round(mean(x));
        Seeds_y(i) =round(mean(y));
    end 
    %______________________________________________

else
    error('Wrong input in function SeedsFiner')
end
% display the Nuclei and their Centers
figure(3);hold on;

for j =1:1:SeedNumber
        plot(Seeds_y(j), Seeds_x(j),'r.','MarkerSize',15);
        h1=plot( Boundary{j}(:,2),Boundary{j}(:,1), 'c','LineWidth',2);
end