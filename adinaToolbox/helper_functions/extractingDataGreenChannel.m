function [data] = extractingDataGreenChannel(dirBase,files,step,nrows,ncols,frames)
% DESCRIPTION:
%    extract the functional image sequence
%
% USAGE:   
%    [data ] = extractingDataGreenChannel(dirBase,files,step,nrows,ncols,frames)
%
%
% ARGUMENTS:
%         dirBase : path of the image files
%         files   : list of files to load
%         step    : number of channels per frame
%         nrows   : number of rows
%         ncols   : number of columns
%         frames  : number of frames
%
% OUTPUT: 
%         data: raw data (nrows x ncols x frames)
%
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".
if nargin ==3
    fname = files(1).name;
    info = imfinfo([dirBase fname]); 
    nrows = info(1).Height;
    ncols = info(1).Width;
    frames = 0;
    for i = 1:length(files),
        fname = files(i).name;
        info = imfinfo([dirBase fname]); 
        frames = frames+length(info)/step;
    end
elseif nargin~=6
    error('Incorrect number of parameters');
end
    
%% Reading the files of the directory (Green/Red Channel)
data = single(zeros(nrows,ncols,frames));
t = 1;
for i = 1:length(files),
    fname = files(i).name;
    info = imfinfo([dirBase fname]); 
    for j = 1:step:length(info),
        data(:,:,t) = single(imread([dirBase files(i).name],j));        
        t = t+1;
    end
end
