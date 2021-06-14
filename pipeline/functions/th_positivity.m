function [th_pos] = th_positivity(cNeurons, imDAPI, imTH)

% segment DAPI image
[nucleus_mask] = nucleus_segm(imDAPI);

% estimate TH positivity for each neuron
[nucleus, Seeds_x, Seeds_y, th_int] = calculate_th_positivity(imTH, nucleus_mask);

% image of nucleus centoids
im_NucSeed = zeros(size(imDAPI,1),size(imDAPI,2));
for i =1:length(nucleus)
    im_NucSeed(Seeds_x(i), Seeds_y(i))=1;
end

im_NeuSeed = zeros(size(imDAPI,1),size(imDAPI,2));
for i =1:length(cNeurons)
    im_NeuSeed(round(cNeurons(i).Centroid(1)), round(cNeurons(i).Centroid(2)))=1;
end

for i =1:length(cNeurons)
    neu_seed = cNeurons(i).Centroid;
    
    if round(neu_seed(1))-50 >= 0 && round(neu_seed(1))+50 <= size(imDAPI,1) ...
            && round(neu_seed(2))-50 >= 0 && round(neu_seed(2))+50 <= size(imDAPI,2)
        
        % find seeds in a 10 x 10 neighborhood 
        [row,col] = find((im_NucSeed(round(neu_seed(1))-50:round(neu_seed(1))+50,...
            round(neu_seed(2))-50:round(neu_seed(2))+50))==1);
        neu_seed_small = [51,51]; % center of the neighborhood
        % calculate Euclidean distance and keep the smallest (closest) one
        dist = 1000000;
        for j = 1 : length(row)
            D  = sqrt(sum((neu_seed_small - [row(j),col(j)]) .^ 2));
            if D < dist
                dist = D;
                index = j;
            end
        end
        if ~isempty(row)
            % find the position of the closest seed in the entire image
            row_idx = round(neu_seed(1))-51 + row(j);
            col_idx = round(neu_seed(2))-51 + col(j);
            % search for the seed's index 
            indx = strmatch([row_idx,col_idx],[Seeds_x;Seeds_y]');
            th_pos(i) = th_int(indx);
        else
            th_pos(i) = 0;
        end
    else
        disp('out of bound')
        th_pos(i) = 0;
    end
end