function Im = load_calcium_data(datadir, imname, dataName, fname)

% This function loads the calcium time-series and saves it in mat format

    MovInf  = imfinfo(dataName);
    Im.T  = numel(MovInf);                                                  % get number of frames
    Im.h  = MovInf(1).Height;                                               % get height
    Im.w  = MovInf(1).Width;                                                % get width
    Im.Np = Im.w*Im.h;
    Im.fname = fname;
    
    Im.DataMat = zeros(Im.w*Im.h,Im.T);                                     % initialize mat to store movie
    for j=1:Im.T
        X = imread(dataName,j);
        Im.DataMat(:,j) = X(:);
    end
    Im.MeanFrame = mean(Im.DataMat,2);
    save([datadir imname],'Im','-v7.3')
