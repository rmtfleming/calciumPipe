
function X = generateSamples2(nDicts,sizeIm,shape,dimension_max,dimension_min)

%%%% Var log:
% 
% generateSamples: first created
% 
% generateSamples2: includes the option of generate realistic form cells
% 
%%%%

X = zeros(sizeIm,sizeIm,nDicts);

switch(shape)
    case 'octogonal'
        for q = 1:nDicts,
            im = zeros(sizeIm);
            i = randi([round(dimension_max(1)*1.1) sizeIm-round(dimension_max(1)*1.1)],1);
            j = randi([round(dimension_max(2)*1.1) sizeIm-round(dimension_max(2)*1.1)],1);
            im(i,j) = 1;
            dimension_rand = [randi([dimension_min(1) dimension_max(1)],1) randi([dimension_min(2) dimension_max(2)],1)];
            sigma = (dimension_rand(1)+dimension_rand(2))/2;
            sigma = round(sigma/(3*2));
            sigma = 3*sigma;
            im = imdilate(im,strel('octagon',sigma));
            h = fspecial('gaussian',[3 3],1);
            im = imfilter(im,h,'same');
            im = im/norm(im,'fro');
            X(:,:,q) = im;
        end
    case 'gaussian'
        for q = 1:nDicts,
            im = zeros(sizeIm);
            i = randi([round(dimension_max(1)*1.1) sizeIm-round(dimension_max(1)*1.1)],1);
            j = randi([round(dimension_max(2)*1.1) sizeIm-round(dimension_max(2)*1.1)],1);
            im(i,j) = 1;
            dimension_rand = [randi([dimension_min(1) dimension_max(1)],1) randi([dimension_min(2) dimension_max(2)],1)];
            sigma_gau(1) = dimension_rand(1)/3;
            sigma_gau(2) = dimension_rand(2)/3;
            h1 = fspecial('gaussian',[dimension_rand(1) 1],sigma_gau(1));
            h2 = fspecial('gaussian',[1 dimension_rand(2)],sigma_gau(2));
            h = conv2(h1,h2);
            im = conv2(im,h,'same');
            im = im/norm(im,'fro');
            X(:,:,q) = im;
        end
    case 'real'
        load cellOutNormalized.mat
        maxSize = max(cellfun('length',cellOutNew));
        for q = 1:nDicts
            im = zeros(512);
            i = randi([ceil(maxSize/2)+1 512-ceil(maxSize/2)+1],1);
            j = randi([ceil(maxSize/2)+1 512-ceil(maxSize/2)+1],1);
            im(i,j) = 1;
            randomCell = cellOutNew{randi(length(cellOutNew),1)};
            im = conv2(im,randomCell,'same');
            if sizeIm~=512
                im = imresize(im,[sizeIm sizeIm]);
            end
            X(:,:,q) = im;
        end
        
        
    otherwise
        errordlg('Unknown shape','Error','modal')
        error('Unknown shape')
end
        