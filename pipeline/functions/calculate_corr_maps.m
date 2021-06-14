function cMaps = calculate_corr_maps(data3D, F)

    dim = size(data3D);
    nRois = size(F,1);
    
    traces = zeros(dim(4),nRois);
    for ii = 1:nRois
        t = F(ii,:)';
        t = t - mean(t(:));
        traces(:,ii) =  t/norm(t);
    end
        
    cMaps = zeros(dim(1)*dim(2), dim(3), nRois, 'single');
    for ii = 1:dim(3)
        s = reshape(squeeze(data3D(:,:,ii,:)),[dim(1)*dim(2) dim(4)]);
        cMaps(:,ii,:) = s * traces;
    end
    cMaps = reshape(cMaps, [dim(1:3) nRois]);
        
end