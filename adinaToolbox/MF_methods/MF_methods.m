function [this] = MF_methods(method)

  if nargin~=1
      error('MF_methods: Number of input parameters is incorrect.')
  end

  
  switch (lower(method))
   case 'spams'
    this.List = {'temporal_sparsity','spatial_sparsity','spatiotemporal_sparsity'};
    for i = 1:length(this.List)
        this.params{i} = feval(this.List{i}); %get_params('spams',this.List{i});
        this.params{i}.type = 'infer_spams';
        %this.params{i} = class(this.params{i},'MF_methods');
    end   
   
    
   case 'nmf'
    this.List = {'nmfnnls','nmfrule','sparsenmfnnls','seminmfnnls','seminmfrule','convexnmfrule',...
                 'orthnmfrule','wnmfrule'};
    for i = 1:length(this.List)
        this.params{i}      = get_params('nmf',this.List{i});
        this.params{i}.type = 'infer_nmf';
        %this.params{i} = class(this.params{i},'MF_methods');
    end      
       
   otherwise
      error('MF_methods: Incorrect Matrix Factorization method.')
  end