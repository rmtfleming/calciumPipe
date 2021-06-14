
% EMD in y direction (rows)
j = 1;
for i = 1 : size(data3D,1)
    [imf{i},ort,nbits] = emd(data3D(i,:,1));
    imff = imf{i};
    data_reconst3(j,:) = imff(3,:);
    j = j+1;
end

%EMD
F = F_tab(:,55);
[imf,ort,nbits] = emd(F');
% calculate HHT
[A,f,tt] = hhspectrum(imf);
% Display HHT
[im,tt,ff] = toimage(A,f);
disp_hhs(im);