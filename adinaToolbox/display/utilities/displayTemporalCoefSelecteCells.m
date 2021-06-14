function [] = displayTemporalCoefSelecteCells(data,Coeff,cNeurons,idxNeurons,idFig)

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

if max(idxNeurons) > nBasis | min(idxNeurons)<1
    error(['Selected indexes are incorrect. They should be inside [1,' num2str(max(idxNeurons)) ']' ]) 
end

if nargin == 4
    idFig = figure;
elseif nargin ~=5,
    error('Number of input parameters is incorrect');
end

str_title = 'Identified Cells';
[idFig,dispImg] = displaySelectedCellsIdentification(data,cNeurons,idxNeurons,idFig);
figure(idFig), title(str_title)


figure,
nRows = length(idxNeurons);
nCols = 2;

maxRange = max([data(:); Coeff(:)]);
t = 1;
for i = 1:nRows
    subplot(nRows,nCols,t), imagesc(dispImg),
    hold on,plot(cNeurons(idxNeurons(i)).obj.mu(2),cNeurons(idxNeurons(i)).obj.mu(1),'.k','MarkerSize',7)
    hold on, contour(cNeurons(idxNeurons(i)).imMask,'Color','w');
    text(cNeurons(idxNeurons(i)).obj.mu(2)+2,cNeurons(idxNeurons(i)).obj.mu(1)+2,num2str(idxNeurons(i)),'FontSize',12,'Color','black','FontWeight','bold');
    axis off, axis square

    subplot(nRows,nCols, t+1),
    plot(shiftdim(data(round(cNeurons(idxNeurons(i)).obj.mu(1)),round(cNeurons(idxNeurons(i)).obj.mu(2)),:)),'r')
    hold on,
    plot(Coeff(:,idxNeurons(i)),'b'),         
    xlabel('time stamps'), axis([0 nFrames 0 maxRange]);
    legend({'Extracted from Data','Inferred from Basis'})
    t = t+2;
end
