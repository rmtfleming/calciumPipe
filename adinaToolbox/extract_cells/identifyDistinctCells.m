function cNeurons = identifyDistinctCells(cNeurons, eNeurons,param)
% DESCRIPTION:
%     check whether two neurons are the same o not
%
% USAGE:   
%       cNeurons = identifyDistinctCells(cNeurons, eNeurons,param)
%
% ARGUMENTS:
%        cNeurons: current neurons
%        eNeurons: detected neurons
%        params  : parameters of similarity rules
%       
% OUTPUT: 
%        cNeurons : update list of detected neurons
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
    thr = param.SimThr;

    nNewNeurons = length(eNeurons);
    nCurrentNeurons = length(cNeurons);
    Cost = zeros(nNewNeurons,nCurrentNeurons);
    Cost2 = zeros(nNewNeurons,nCurrentNeurons);
    Cost3 = zeros(nNewNeurons,nCurrentNeurons);
    for i = 1:nNewNeurons,
        for j = 1:nCurrentNeurons,          
            switch lower(param.Tsimilarity)
                case 'jointdistance'
                    Cost(i,j) = jointDistance(eNeurons(i).obj.mu',eNeurons(i).obj.Sigma,cNeurons(j).obj.mu',cNeurons(j).obj.Sigma);
                case 'mahalanobisdistance'
                    Cost(i,j) = MahalanobisDistance(eNeurons(i).obj.mu',eNeurons(i).obj.Sigma,cNeurons(j).obj.mu',cNeurons(j).obj.Sigma);
                case 'crosscorrelation'
                    Cost(i,j) =  cNeurons(j).imOrig(:)'*eNeurons(i).imOrig(:)/(norm(cNeurons(j).imOrig(:))*norm(eNeurons(i).imOrig(:)));
                    Cost2(i,j) = sum(cNeurons(j).imMask(:).*eNeurons(i).imMask(:))/sum(cNeurons(j).imMask(:)); 
                    Cost3(i,j) = sum(cNeurons(j).imMask(:).*eNeurons(i).imMask(:))/sum(eNeurons(i).imMask(:));
                otherwise
                    error('Error')
            end
            % Cost(i,j) = feval(str2func(['@' param.Tsimilarity]),eNeurons(i).obj.mu',eNeurons(i).obj.Sigma,cNeurons(j).obj.mu',cNeurons(j).obj.Sigma);
        end
    end
    
    for i = 1:nNewNeurons,
        [val,pos] = max(Cost(i,:));
        [val2,pos2] = max(Cost2(i,:));
        [val3,pos3] =  max(Cost3(i,:));
        if val > thr
            cNeurons(pos) = fuseNeurons(cNeurons(pos),eNeurons(i));
        else
            if val3 > 0.6,
                cNeurons(pos3) = fuseNeurons(cNeurons(pos3),eNeurons(i));
            elseif val2 >0.6
                cNeurons(pos2) = fuseNeurons(cNeurons(pos2),eNeurons(i));             
            else
                cNeurons = [ cNeurons eNeurons(i)];
            end
        end
    end
end


