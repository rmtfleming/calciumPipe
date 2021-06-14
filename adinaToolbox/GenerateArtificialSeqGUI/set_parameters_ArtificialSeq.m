function [this] = set_parameters_ArtificialSeq(this,param)

if nargin == 2,
    if mod(length(param), 2) ~=0
        error('Parse_input_parameter: Input parameters must be given in pairs (name and value)!');
    end;
else
    error('Incorrect number of input parameters')
end


flag1 = false;
flag2 = false;
    
i = 1;

while (i <= length(param))
    
    %     if ischar(param{i+1})
    %         param{i+1} = str2double(param{i+1});
    %     end;
    
    switch lower(param{i})
        
        case 'ndicts'
            this.nDicts             =  max(param{i+1},1); % min: 1 dict
            [this] = set_parameters_ArtificialSeq(this,{'numGroups',this.numGroups});
            
        case 'nrows'
            this.nrows              =  max(param{i+1},1);
            
        case 'framerate'
            this.frameRate          =  max(param{i+1},1);
            
        case 'seqduration'
            this.seqDuration        =  max(param{i+1},0);
            
        case 't_rise_min'
            this.t_rise_min         =  max(param{i+1},0);
            
        case 't_rise_max'
            this.t_rise_max         =  max(param{i+1},0);
            
        case 'tau_min'
            this.tau_min            =  max(param{i+1},0);
            
        case 'tau_max'
            this.tau_max            =  max(param{i+1},0);
            
        case 'thr_cells'
            this.thr_cells          =  min(max(param{i+1},0),1);
            
        case 'thr_seq'
            this.thr_seq            =  min(max(param{i+1},0),1);
            
        case 'numgroups'
            this.numGroups          =  min(max(param{i+1},1),this.nDicts);
            flag1 = true;
            numGroups_tmp = param{i+1};
            
        case 'percentmult'
            this.percentMult        =  min(max(param{i+1},0),1);
            
        case 'maxmult'
            this.maxMult            =  min(max(param{i+1},1),this.numGroups);
            flag2 = true;
            maxMult_tmp = param{i+1};
            
        case 'shape'
            this.shape              =  lower(param{i+1});
            if ~isequal(this.shape,'gaussian') && ~isequal(this.shape,'octogonal') ...
                    && ~isequal(this.shape,'real')
                this.shape = 'gaussian';
            end
            
        case 'temporal_pattern'
            this.temporal_pattern   =  lower(param{i+1});
            if ~isequal(this.temporal_pattern,'gauss_exp') && ~isequal(this.temporal_pattern, ...
                    'slope') && ~isequal(this.temporal_pattern,'exp_exp')
                this.temporal_pattern = 'exp_exp';
            end
            
        case 'rho_x_min'
            this.rho_x_min          =  min(max(param{i+1},0),0.2);
%             if this.roxmin>param{i+1}
%                 this.roxmax = min(max(param{i+1},0.0001),0.2);
%                 this.roxmin = this.roxmax;
%             end
%             
        case 'rho_x_max'
            this.rho_x_max          =  min(max(param{i+1},0),0.2);
%             if param{i+1}>this.roxmax
%                 this.roxmin = min(max(param{i+1},0.0001),0.2);
%                 this.roxmax = this.roxmin;
%             end
%             
        case 'rho_y_min'
            this.rho_y_min          =  min(max(param{i+1},0),0.2);
%             if this.roymin>param{i+1}
%                 this.roymax = min(max(param{i+1},0.0001),0.2);
%                 this.roymin = this.roymax;
%             end
%             
        case 'rho_y_max'
            this.rho_y_max          =  min(max(param{i+1},0),0.2);
%             if param{i+1}>this.roymax
%                 this.roymin = min(max(param{i+1},0.0001),0.2);
%                 this.roymax = this.roymin;
%             end
            
        case 'ra'
            this.RA                 =  max(param{i+1},0.0001);
            
        case 'clipping'
            this.clipping           =  min(max(param{i+1},0),1);
            
        case 'min_prob_candidate'
            this.min_prob_candidate =  min(max(param{i+1},0),1);
            
        case 'time_fuse'
            this.time_fuse          =  max(param{i+1},0);
            
        case 'prct_artifacts'
            this.prct_artifacts     =  max(param{i+1},0);
            
        case 'ampl_artifacts_min'
            this.ampl_artifacts_min =  min(max(param{i+1},0),1);
            
        case 'ampl_artifacts_max'
            this.ampl_artifacts_max =  min(max(param{i+1},0),1);
            
        otherwise
            error('The variable %s does not exist\n',param{i});
            
    end
    
    i = i+2;
end

if flag1
    this.numGroups = min(max(numGroups_tmp,1),this.nDicts);
    this.maxMult = min(max(this.maxMult,1),this.numGroups);
end

if flag2
    this.maxMult = min(max(maxMult_tmp,1),this.numGroups);
end
