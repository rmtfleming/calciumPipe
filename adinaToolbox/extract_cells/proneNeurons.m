function [cNeurons,Coeff_single,Dicts_single,ReconsDict_single] = proneNeurons(cNeurons,Coeff_single,Dicts_single,ReconsDict_single)
% DESCRIPTION:
%     remove noisy neurons
%
% USAGE:   
%       [cNeurons,Coeff_single,Dicts_single,ReconsDict_single] = proneNeurons(cNeurons,Coeff_single,Dicts_single,ReconsDict_single)
%
% ARGUMENTS:
%         cNeurons : current set of single and unique neurons
%         Coeff_single : activation of the previous neurons
%         Dicts_single : matrix representation of the neurons
%         ReconsDict_single : image representation of the neurons
%       
% OUTPUT: 
%         cNeurons : updated set of single and unique neurons
%         Coeff_single : updated activation of the previous neurons
%         Dicts_single : updated matrix representation of the neurons
%         ReconsDict_single : updated image representation of the neurons
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

nActivations = sum(Coeff_single>0,1);
idx = find(nActivations <= 3); % param.actThr

cNeurons(idx) =[];
Coeff_single(:,idx) = [];
Dicts_single(:,idx) = [];
ReconsDict_single(:,:,idx) = [];



