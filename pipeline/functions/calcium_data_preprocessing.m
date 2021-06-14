function [data_preproc, param_preproc] = calcium_data_preprocessing(data3D)

% this function uses wavelet transform (from the ADINA toolbox) to deniose  calcium time-series

    % Parameters for PreProcessing Step
   param_preproc = initial_parameter_preprocessing;

    % run preprocessing
    data_tmp = extractingNormalizedSeq(data3D(:,:,:),param_preproc);
    data_preproc = data_tmp.data_denoise;