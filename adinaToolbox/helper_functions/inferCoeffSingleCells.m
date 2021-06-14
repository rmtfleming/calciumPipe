function [Coeffm,Dicts,ReconsDict,ReconsSeq] = inferCoeffSingleCells(data,cNeurons,param)
% data should be preprocessed using preprocessingData

if nargin <3,
    param = initial_parameter_SPAMS;
elseif nargin == 3 && ~isstruct(param),
    param_init = initial_parameter_SPAMS;
    param = set_parameters_SPAMS(param_init,param);
elseif ~isstruct(param)
    error('Incorrect number of input parameters')
end

[nrows,ncols,nFrames] = size(data);
data       =  reshape(data,[nrows*ncols nFrames]);
Dicts = [];
nBasis = length(cNeurons);

for i = 1:nBasis,
    Dicts(:,i) = preprocessingData(cNeurons(i).imOrig(:),'non','non');
end

if sum(class(Dicts) ~= class(data))
    Dicts = cast(Dicts,class(data));
end

Coeffm     =  full(mexLasso(data,Dicts,param))';

if nargout >= 3,
    for i = 1:nBasis,
        ReconsDict(:,:,i) = reshape(Dicts(:,i),[nrows ncols]);
    end
    if nargout >= 4,
        ReconsSeq = reshape(Dicts*Coeffm',[nrows ncols nFrames]);
    end
end