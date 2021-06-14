function [] = displayExploreCells(data,Coeff,cNeurons,param,idFig)
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
[nFramesU,nBasisU] = size(Coeff);
nBasisN = length(cNeurons);

if nFrames ~= nFramesU,
    error('Temporal dimension between data and Coeff are not the same')
end
if nBasisN ~= nBasisU,
    error('Number of Basis between Coeff and cNeurons are not the same')
else
    nBasis = nBasisN;
end

if nargin == 3 || nargin == 4
    idFig = figure;
elseif nargin ~=5,
    error('Number of input parameters is incorrect');
end

str_title = 'Identified Cells';
[idFig,dispImg] = displaySingleCellsIdentification(data,cNeurons,idFig);
figure(idFig), title(str_title)

params = [];
params.str_title  = str_title;
params.Coeff      = Coeff;
params.cNeurons   = cNeurons;
params.data       = data;
params.rangeCoeff = [min(Coeff(:)) max(Coeff(:))];
params.rangeData  = [min(data(:)) max(data(:))];
params.dispImg     = dispImg;
params.initF       = initF;
params.endF        = endF;
params.step        = step;
maskCells = [];
for i = 1:nBasis
    maskCells(:,:,i) = cNeurons(i).imMask;
end

params.maskCells = maskCells;

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
   cNeurons = params.cNeurons;
   
   initF = params.initF;
   endF = params.endF;
   step  = params.step;
   
   idx = find(params.maskCells(row,col,:)>0);
   
   nRows = length(idx);
   nCols = 2;
   
   figure(hTime),
   t = 1;
   for i = 1:nRows,
        subplot(nRows,nCols,t), imagesc(params.dispImg),
        
        hold on,plot(cNeurons(idx(i)).obj.mu(2),cNeurons(idx(i)).obj.mu(1),'.k','MarkerSize',7)
        hold on, contour(cNeurons(idx(i)).imOrig,1,'Color','w');
        text(cNeurons(idx(i)).obj.mu(1)+2,cNeurons(idx(i)).obj.mu(2)+2,num2str(idx(i)),'FontSize',12,'Color','black','FontWeight','bold');
        axis off, axis square, axis equal
        
        subplot(nRows,nCols, t+1),
       
        plot(initF:step:endF,cNeurons(idx(i)).imOrig(:)'*reshape(params.data(:,:,initF:step:endF),[],length(initF:step:endF))/sum(cNeurons(idx(i)).imOrig(:)),'r')
        %plot(initF:step:endF, shiftdim(params.data(round(cNeurons(idx(i)).obj.mu(1)),round(cNeurons(idx(i)).obj.mu(2)),initF:step:endF)),'r')
        hold on,
        plot(initF:step:endF,params.Coeff(initF:step:endF,idx(i)),'b'),         
        xlabel('time stamps'), axis([0 size(params.Coeff,1) 0 max(params.rangeData(2),params.rangeCoeff(2))]); hold off;
        legend({'Extracted from Data','Inferred from Basis'})
        t= t+2;
   end
