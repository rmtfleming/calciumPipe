function displayPerformanceResults(name,SPSNA)
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
        disp(['Results of ' name])
        disp(['Sensitivity         : ' num2str(SPSNA(1)) '%'])
        disp(['Precision           : ' num2str(SPSNA(2)) '%'])
        disp(['Specifity           : ' num2str(SPSNA(3)) '%'])
        disp(['Negative Predictive : ' num2str(SPSNA(4)) '%'])
        disp(['Accuracy            : ' num2str(SPSNA(5)) '%'])
        disp('')