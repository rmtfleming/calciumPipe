function [] = displayExploreCalciumImageTwo_bis(data,data1,image,param,idFig)
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
if nargin <4;
    initF = 1;
    endF  = size(data,3);
    step = 1;
elseif nargin == 4
    lengthSeq = size(data,3);
    initF = min(max(param.initF,1),lengthSeq);
    endF  = min(max(param.endF,initF),lengthSeq);
    step  = min(max(param.step,1),endF-initF+1);
end


[nrows,ncols,nFrames] = size(data);

if nargin == 3 || nargin == 4
    idFig = figure; 
elseif nargin ~=5,
    error('Number of input parameters is incorrect');
end

figure(idFig), imagesc(image); axis off equal


params = [];
params.data       = data;
params.data1       = data1;
params.rangeData  = [min(data(:)) max(data(:))];
params.initF       = initF;
params.endF        = endF;
params.step        = step;

set(idFig,'UserData',params);
set(idFig,'WindowButtonDownFcn',@ButtonDownCallback);

function ButtonDownCallback(src,evedentdata)

   persistent hTime
   
   f1 = src;
   params = get(f1,'UserData');
   if (isempty(hTime))
        hTime = figure;        
   end 
   [f_cp, a1_cp] = pointer2d(f1);
   col = round(a1_cp(1)); %dictionary element
   row = round(a1_cp(2)); %frame
   
   
   initF = params.initF;
   endF = params.endF;
   step  = params.step;
   
   figure(hTime), subplot(3,1,1),
   plot(initF:step:endF, shiftdim(params.data(row,col,initF:step:endF),1),'b')
   %subplot(2,1,2),
   hold on, 
   plot(initF:step:endF, shiftdim(params.data1(row,col,initF:step:endF),1),'r'), hold off
   subplot(3,1,2)
   plot(initF:step:endF, shiftdim(params.data(row,col,initF:step:endF),1)-shiftdim(params.data1(row,col,initF:step:endF),1),'r')
   subplot(3,1,3)
   plot(initF:step:endF, shiftdim(params.data(row,col,initF:step:endF),1)./shiftdim(params.data1(row,col,initF:step:endF),1),'r')
   
    
  