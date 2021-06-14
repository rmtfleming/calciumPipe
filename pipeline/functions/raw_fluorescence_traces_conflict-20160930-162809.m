function time_trace = raw_fluorescence_traces(data3D, cNeurons)

dim = size(data3D);
if length(dim)==4
    data3D = reshape(data3D,[],dim(4));
else if length(dim)==3
        data3D = reshape(data3D,dim(1)*dim(2),dim(3));
end

for i = 1 : length(cNeurons)
    time_trace(i,:) = nanmean(data3D(cNeurons(i).imMask(:)==1,:),1);
end