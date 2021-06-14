function [neuron] = fuseNeurons(cNeuron,eNeuron)
% DESCRIPTION:
%     fuse two identical neurons
%
% USAGE:   
%       [neuron] = fuseNeurons(cNeuron,eNeuron)
%
% ARGUMENTS:
%        cNeuron: current neuronal
%        eNeuron: detected approx. identical neuron
%       
% OUTPUT: 
%        neuron: updated neuron
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;

%     [nrows,ncols] = size(neuron.im);
    tmpImg = mean(cat(3,cNeuron.imOrig,eNeuron.imOrig),3);
    imMask_tmp = tmpImg>(0.01*max(tmpImg(:)));
    [I,J] = find(imMask_tmp>0);
    X = [[J,I]];
    obj = gmdistribution.fit(X,1);
    
    neuron.imMask = double(imMask_tmp);
    neuron.im = double(tmpImg.*imMask_tmp);
    neuron.imOrig =  single(neuron.im);
    neuron.Area = sum(neuron.imMask(:)>0);
    neuron.maxF = max(cNeuron.maxF,eNeuron.maxF);
    neuron.X = X;
    neuron.obj.mu(1)=obj.mu(2);
    neuron.obj.mu(2)=obj.mu(1);
    neuron.obj.Sigma=obj.Sigma;
%     neuron.X_v =  tmpImg(imMask_tmp>0);
    neuron.Centroid = neuron.obj.mu;
    neuron.time = [];
               
end
