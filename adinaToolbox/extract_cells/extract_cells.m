function [this] = extract_cells(varargin)

%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".;
  error(nargchk(0, 2, length(varargin)));

  
  switch (length(varargin))
   case 0
    tmp = initial_parameter_Segmentation;
    this = extract_cells(tmp.type,tmp);
    %this = class(this, 'extract_cells');
    
   case 1
    if isa(varargin{1}, 'extract_cells')
      this = other;
    else
      switch (varargin{1})
       case 'wavelet_denoising'
        this.type = 'wavelet_denoising';
       case 'image_smoothing'
        this.type = 'image_smoothing';
        case 'image_smoothing_powerlaw'
            this.type = 'image_smoothing_powerlaw';
       otherwise
        error('Invalid extract cell function type.');
      end     
      this.param  = initial_parameter_Segmentation;
      %this = class(this, 'extract_cells');
    end

   case 2
    switch (varargin{1})
       case 'wavelet_denoising'
        this.type = 'wavelet_denoising';
       case 'image_smoothing'
        this.type = 'image_smoothing';
        case 'image_smoothing_powerlaw'
            this.type = 'image_smoothing_powerlaw';
       otherwise
        error('Invalid extract cell function type.');
    end  

    this.param = varargin{2};
    %this = class(this, 'extract_cells');
  end