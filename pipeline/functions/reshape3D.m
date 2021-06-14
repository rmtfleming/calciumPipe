function data3D = reshape3D(Im)

% This function reshapes a m x n dataset into 3D dataset h x w x T

data3D = zeros(Im.h, Im.w, Im.T);

for i=1:Im.T
    img = reshape(Im.DataMat(:,i), Im.h, Im.w);
    data3D(:,:,i) = img(:,:);
end