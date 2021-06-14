function [data_bis] = removingLinesMotionArtifacts(data)
% DESCRIPTION:
%    remove line artifacts introduced by line-to-line motion compensation
%    algorithms
%
% USAGE:   
%       [data_bis] = removingLinesMotionArtifacts(data)
%
%
% ARGUMENTS:
%         data : raw data
%
% OUTPUT: 
%         data_bis : corrected data
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".
[nrows,ncols,nFrames] = size(data);
[X,Y] = meshgrid(1:ncols,1:nrows);
warning off
data_bis = zeros(nrows,ncols,nFrames);
for i = 1:nFrames,
   tmpImg =data(:,:,i);
   maskLines = imerode(tmpImg==0,ones(1,21));
   Z =tmpImg; 
   intPoint = (imdilate(maskLines,ones(3,1))-maskLines);
   F = TriScatteredInterp(X(intPoint==1),Y(intPoint==1),tmpImg(intPoint==1));
   try
        Z(maskLines~=0) = F(X((maskLines~=0)),Y((maskLines~=0)));
        Z(isnan(Z))=0; Z(isinf(Z))=0;
   catch   	
        tt = imdilate(Z,ones(5,1));
        Z(maskLines~=0) = tt(maskLines~=0);

   end
   data_bis(:,:,i) = Z;  
end
warning in
