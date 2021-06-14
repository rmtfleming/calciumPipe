function [cNeurons, param, Coeff_s, param_preproc] = detectNeurons(data)

% This function detects individual neurons from calcium imaging data using
% a matrix factorization method. Sparse dictionary learning (ADINA toolbox)
% or Non-negative matrix factorisation
% 
% INPUTS:
%     data: calcium imaging data
%     method: matrix factorization method to be used. 'spams' or 'nmf'
%     
% OUTPUTS:
%     cNeurons: detected neurons


% data preprocessing
param_preproc = initial_parameter_preprocessing;

% run preprocessing
data_tmp = extractingNormalizedSeq(data(:,:,:),param_preproc);
data_preproc = data_tmp.data_denoise;
X = data_preproc;

% SPAMS decomposition
param_spams = initial_parameter_SPAMS;
[Coeff,~,ReconsDict] = infer_spams(X,1,param_spams);
    
% exract single cells
param = initial_parameter_Segmentation; 
param = param.all;

cNeurons = [];
if ~isempty(data)
    if ~isempty(ReconsDict)
        disp('Extracting single cells ... '),
        if ~isempty(cNeurons) && param.resetNeurons
             % ask if they want to update the current extracted cells                 
             [cNeurons] = extractSingleCells(Coeff,ReconsDict,cNeurons,param);
        else                
             [cNeurons] = extractSingleCells(Coeff,ReconsDict,[],param);
        end

        if ~isempty(cNeurons)                
            disp('Extracted Single Cells is done.')
            disp('Computing Temporal Evolution for the extracted single Cells')
            [Coeff_single,Dicts_single,ReconsDict_single] = inferCoeffSingleCells(X,cNeurons,param);

            %% filter neurons that are not used for reconstructing
            [cNeurons,Coeff_single,Dicts_single,ReconsDict_single] = proneNeurons(cNeurons,Coeff_single,Dicts_single,ReconsDict_single);

%                 cNeurons = identifyDistinctCells(cNeurons, cNeurons,param);

            disp('Temporal Evolution is computed')    
        else
            disp('No found single cells.')
        end

    else
         msgbox('Please run Matrix Factorization algorithm.','Error','error','modal');
    end
else
    msgbox('Please load a sequence before setting the number of basis functions.','Error','error','modal');
end

Coeff_s = Coeff_single;

end