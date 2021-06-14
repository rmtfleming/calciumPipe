function [param_out] = initial_parameter_ArtificialSeq()

param = [];
param.nDicts = 10; % number of cells
param.nrows = 512; % image size
param.frameRate = 30; % frame rate in fps
param.seqDuration = 15; % sequence duration in seconds
param.t_rise_min = 200; %400; % neuron activation rise time in milliseconds
param.t_rise_max = 200;
param.tau_min = 0.5; %0.8; % decay parameter of the cells --> exp(-t/tau)
param.tau_max = 0.5;
param.thr_cells = 0.3; % threshold of overlapping between cells
param.thr_seq = 0.2; % threshold of overlapping between time sequences
param.numGroups = 5; % number of cell clusters. MUST BE LOWER OR EQUAL THAN "NDICTS"
param.percentMult = 0.3; % the percentage of cells that will be multiple assigned
param.maxMult = 3; % maximum number of multiplicities per cell
param.shape = 'real'; % shape of the cells. Options: 'octogonal', 'gaussian', 'real'
param.temporal_pattern = 'exp_exp'; % temporal activation pattern. Options: 'slope', 'gauss_exp', 'exp_exp'
param.rho_x_min = 0.04; % rho_x determines the min and max x-size of cells w.r.t. the image size. (0<ro<0.2. optimal 0.05)
param.rho_x_max = 0.05;
param.rho_y_min = 0.04; % rho_y determines the min and max y-size of cells w.r.t. the image size. (0<ro<0.2. optimal 0.05)
param.rho_y_max = 0.05;
param.RA = 1; % Relative Amplitude of the sequence generated
param.clipping = 0; % clipping percentage w.r.t. the maximum value. The values of the activation pattern below max*clipping will be clipped to 0.
param.min_prob_candidate = 1; %0.6; % minimum probability for the cells onset in a group
param.time_fuse = 1; % the cell excitations which are closer than time_fuse will be fused into a single one.
param.prct_artifacts = 0.01; % percentage of artifacts w.r.t. the number of cells (nDicts)
param.ampl_artifacts_min = 0.2; % min and max amplitude of artifacts. Range = [0 1]
param.ampl_artifacts_max = 0.5;

param_out.all = param;
param_out.ListComplete = {'nDicts','nrows','frameRate','seqDuration','t_rise_min','t_rise_max','tau_min','tau_max',...
    'thr_cells','thr_seq','numGroups','percentMult','maxMult','shape','temporal_pattern',...
    'rho_x_min','rho_x_max','rho_y_min','rho_y_max','RA','clipping','min_prob_candidate','time_fuse','prct_artifacts',...
    'ampl_artifacts_min','ampl_artifacts_max'};
param_out.List = {'nrows','frameRate','t_rise_min','t_rise_max','tau_min','tau_max','thr_cells','thr_seq','percentMult',...
    'maxMult','rho_x_min','rho_x_max','rho_y_min','rho_y_max','RA','clipping','min_prob_candidate','time_fuse','prct_artifacts',...
    'ampl_artifacts_min','ampl_artifacts_max'};
