function [X_new,X_new_clean,assemblies,X,Ucells_GT] = generateArtificialSequence_v6(varargin)

% function "generateArtificialSequences" returns:
% 
%   X_new: nrows x nrows x nSamples matrix. It contains the noisy video
% sequence of the neurons.
%   X_new_clean: nrows x nrows x nSamples matrix. It contains the clean video
% sequence of the neurons.
%   assemblies: cell array which contains the cells that belong to each
% group.
%   X: dictionary of cells
%   Ucells_GT: temporal evolution of the cells
% 
% Its input arguments can be:
% generateArtificialSequence(nDicts): the user specifies the number of cells
% generateArtificialSequence(nDicts,seqDuration): the user also specifies the sequence duration (in seconds)
% generateArtificialSequence(nDicts,seqDuration,{'key',value,(...)}): the user also specifies some parameters


%%%% Var log:
% 
% generateArtificialSequence: creation
% 
% generateArtificialSequence_v2: the decay function is computed for each cell
%     to avoid cases in which a cell is excited in contiguous groups and therefore
%     it presents more than one local maxima. So in v2 if a cell is excited
%     in two or more contiguous groups, the decay factor is reestimated in 
%     order to cover all the excitation time.
%     In addition, the values of the activation pattern which are below a
%     percentage (w.r.t. the maxima) chosen by the user are supressed.
%     The output variable "Uexp" is substituted by "Ucells_GT", which contains
%     the temporal behavior of each cell, instead of the group one.
% 
% generateArtificialSequence_v3: the decayCt variable has been redefined as tau,
%     which expresses the decay time in seconds.
%     The clipping of low activation pattern values is now computed automatically
%     for these values which are below 3*tau. These values are not eliminated,
%     they are equal to the last non-eliminated value.
%     
% generateArtificialSequence_v4: the activation pattern is made by concatenate two
%     exponential functions, instead of a gaussian and an exponential.
% 
% generateArtificialSequence_v5: some code has been packed into independent
%     functions.
%     
% generateArtificialSequence_v6: the initial time of the excitation pattern is
%     randomly shifted for each individual cell fitting a normal distribution.
%     The t_rise and tau are randomly selected between a max and a min
%     especified by the user.
%     The different cells in a group have an onset probability depending on
%     the number of cells in that group.
%     The length of the sequence is subdivided into k "subsequences" if
%     #numgroups*(t_rise+tau) < length_seq.
%     If two excit. patterns are closer than time_fuse, this two
%     excitations are fused into a single one.
%%%%

% DEFAULT PARAMETERS SETTINGS:
% 
% nDicts             = 10;         % number of cells
% nrows              = 512;        % image size
% frameRate          = 30;         % frame rate in fps
% seqDuration        = 15;         % sequence duration in seconds
% t_rise             = [200 400];  % neuron activation rise time in milliseconds
% tau                = [0.5 0.8];  % decay parameter of the cells --> exp(-t/tau)
% thr_cells          = 0.5;        % threshold of overlapping between cells
% thr_seq            = 0.2;        % threshold of overlapping between time sequences
% numGroups          = 5;          % number of cell clusters. MUST BE LOWER OR EQUAL THAN "NDICTS"
% percentMult        = 0.4;        % the percentage of cells that will be multiple assigned
% maxMult            = 3;          % maximum number of multiplicities per cell
% shape              = 'real';     % shape of the cells. Options: 'octogonal', 'gaussian', 'real'
% temporal_pattern   = 'exp_exp';  % temporal activation pattern. Options: 'slope', 'gauss_exp', 'exp_exp'
% rho_x              = [0.04 0.05];% rho_x determines the min and max x-size of cells w.r.t. the image size. (0<ro<0.2. optimal 0.05)
% rho_y              = [0.04 0.05];% rho_y determines the min and max y-size of cells w.r.t. the image size. (0<ro<0.2. optimal 0.05)
% RA                 = 10;         % Relative Amplitude of the sequence generated
% clipping           = 0;          % clipping percentage w.r.t. the maximum value. The values of the activation pattern below 
%                                  % max*clipping will be clipped to 0.
% min_prob_candidate = 0.3;        % minimum probability for the cells onset in a group
% time_fuse          = 1;          % the cell excitations which are closer than time_fuse will be fused into a single one.
% prct_artifacts     = 1;          % percentage of artifacts w.r.t. the number of cells (nDicts)



%% USER PARAMETERS SECTION
tic

%%%%%%%%%%%%%%%%%%%%%%%%%% USER PARAMETERS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[param_out] = initial_parameter_ArtificialSeq();
param = param_out.all;
fields = fieldnames(param);

for i = 1:length(fields),
eval([fields{i} ' = param.' fields{i} ';']);
end

if nargin == 1
    nDicts = varargin{1};
    if ~isa(nDicts,'double')
        msgbox('Incorrect argument types!','Error','error','modal')
        return
    end

elseif nargin == 2
    nDicts = varargin{1};
    seqDuration = varargin{2};
    if ~isa(nDicts,'double') || ~isa(seqDuration,'double')
        msgbox('Incorrect argument types!','Error','error','modal')
        return
    end

    
elseif nargin == 3
    
    param.nDicts = varargin{1};
    param.seqDuration = varargin{2};
    if ~isa(param.nDicts,'double') || ~isa(param.seqDuration,'double')
        msgbox('Incorrect argument types!','Error','error','modal')
        return
    end
    
    if ~isstruct(varargin{3})
        param = set_parameters_ArtificialSeq(param,varargin{3});
    else
        param = varargin{3};
    end    
    
    for i = 1:length(fields),
        eval([fields{i} ' = param.' fields{i} ';']);
    end
       
    
elseif nargin > 3
    msgbox('Too many input arguments!','Error','error','modal')
    return
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% WARNING SECTION
rho_x = [rho_x_min rho_x_max];
rho_y = [rho_y_min rho_y_max];
dimension_max = [ceil(nrows*max(rho_x)) ceil(nrows*max(rho_y))]; % maximum dimension [x y] of the cells in pixels
dimension_min = [ceil(nrows*min(rho_x)) ceil(nrows*min(rho_y))]; % minimum dimension [x y] of the cells in pixels
aux = [t_rise_min t_rise_max];
t_rise(1) = min(aux);
t_rise(2) = max(aux);
aux = [tau_min tau_max];
tau(1) = min(aux);
tau(2) = max(aux);
clear aux
nSamples = frameRate*seqDuration; % temporal samples

%%%%%% Spatial warning section %%%%%%%


sizeCells = dimension_max(1)*dimension_max(2);

areaCovered = sizeCells*nDicts*(1-thr_cells) + (thr_cells^2)*sizeCells;

if areaCovered > (nrows^2)
    errordlg('The spatial parameters do not match. Check whether the number of cells fits into the spatial samples!','Error','modal')
    error('Spatial parameters do not match');
elseif 2*areaCovered > (nrows^2)
    %msgbox('Your spatial parameters just fit. It migth take several time to converge...','Warning','warn','modal')
end

clear areaCovered;
clear sizeCells;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%% Temporal warning section %%%%%%%

lenDec = 3*tau(2); % lenght (in seconds) of the decay period
lenDecFr = round(lenDec*frameRate); %length in frames of the decay period
lenRise = round((t_rise(2)/1000)*frameRate); %length (in frames) of the rise period

lenSeq = lenDecFr+lenRise;
areaCovered = lenSeq*(1-thr_seq)*numGroups + thr_seq*lenSeq;

if areaCovered > nSamples
    errordlg('The temporal parameters do not match. Check whether the number of clusters and its active period duration fits into the temporal samples!'...
        ,'Error','modal')
    error('Temporal parameters do not match');
elseif 2*areaCovered > nSamples
    %msgbox('Your temporal parameters just fit. It migth take several time to converge...','Warning','warn','modal')
end

clear lenDec;
clear lenDecFr;
clear lenRise;
clear lenSeq;
clear areaCovered;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SAMPLES GENERATION
X = generateSamples2(nDicts,nrows,shape,dimension_max,dimension_min); %generate cells arbitrarily allocated
A = zeros(nDicts); 
final = 1;

while final,
    for q = 1:nDicts,
        for z = 1:nDicts,
            % Sum of all the cells obtained by direct product between cell q and z, divided
            % by the sum of all the cells of q, in order to obtain the overlapping
            % percentage between cells q and z
            A(q,z) = sum(sum(X(:,:,q).*X(:,:,z)>0))/sum(sum(X(:,:,q)>0));
        end
        A(q,q) = 0; % q cell is 100% overlapped by itself. We do not want to eliminate it!
    end
    
    % delete the cells whose overlapping level exceeds 'thr_cells'
    [ii,~] = find(A>thr_cells);
    if ~isempty(ii)        
        listDel = unique(ii);
        X(:,:,listDel) = [];       
        % generate as many new cells as eliminated
        Xnew = generateSamples2(length(listDel),nrows,shape,dimension_max,dimension_min);
        X = cat(3,X,Xnew);       
    else
        final = 0;
        for q = 1:nDicts,
            for z = 1:nDicts,
                A(q,z) = sum(sum(X(:,:,q).*X(:,:,z)>0))/sum(sum(X(:,:,q)>0));
            end
        end
        A(A>=thr_cells) = 1;
        A(A<thr_cells) = 0;

        % all the pair of cells which are overlapped have a "1" in the
        % matrix A
    end
end

disp('Samples generated')

%% CLUSTERS DEFINITION
% the following function returns flagNeurons, which is a matrix with the
% (multiple) assigned group(s) for each cell; and assemblies cell array,
% which contains the cells that belong to each group. The function assures
% at least one cell per assembly.
[~,assemblies] = assignCellClustersMultipleMult(numGroups,X,A,percentMult,maxMult);
disp('Cells assigned to clusters')


%% TEMPORAL DIMENSION

beta = 0.05; % security margin factor
reps_double = seqDuration/((t_rise(2)/1000+3*tau(2))*numGroups);
reps = round(reps_double);
if reps*(1+beta)>reps_double
    reps = reps-1;
end
Ucells_GT_merged = [];
fprintf('Sequence cell split in %d subsequences\n',reps);

activationPattern = generateActivationPattern(temporal_pattern,tau,t_rise,frameRate,clipping);
disp('Activation pattern generated')

lengthSequence = length(activationPattern);
nSamples_repint = round((seqDuration/reps)*frameRate);

for k = 1:reps
    U = zeros(nSamples_repint,numGroups);
    for i = 1:numGroups,
        posDelta = randi([1 nSamples_repint-lengthSequence],1);
        U((posDelta:(posDelta+lengthSequence-1)),i) = 1;
        % random distribution of the sequences along the time axis
    end
    
    final = 1;
    A = zeros(numGroups);
    while final,
        % Sum of all the sequences obtained by direct product between seq. q and z, divided
        % by the sum of seq. q, in order to obtain the overlapping
        % percentage between sequences q and z
        for q = 1:numGroups,
            for z = 1:numGroups,
                A(q,z) = sum(sum(U(:,q).*U(:,z)>0))/sum(sum(U(:,q)>0));
            end
            A(q,q) = 0;
        end
        [ii,~] = find(A>thr_seq);
        if ~isempty(ii)
            listDel = unique(ii);
            for q = 1:length(listDel)-1,
                posDelta = randi([1 nSamples_repint-lengthSequence],1);
                U(:,listDel(q)) = 0;
                U((posDelta:(posDelta+lengthSequence-1)),listDel(q)) = 1;
            end
        else
            final = 0;
        end
    end
    
    %% Activation pattern assignment
    
    Ucells_temporalPattern = zeros(nSamples_repint,size(X,3));
    for g = 1:length(assemblies)
        selected_cells = selectCandidates(assemblies{g},param.min_prob_candidate);
        active = find(U(:,g)>0);
        Ucells_temporalPattern(active,selected_cells) = 1;
    end
    
    
    Ucells_GT = assignActivationPatternToCells_v2(temporal_pattern,tau,t_rise,frameRate,clipping,Ucells_temporalPattern,activationPattern);
    Ucells_GT_merged = [Ucells_GT_merged; Ucells_GT];
end
Ucells_GT = Ucells_GT_merged;
Ucells_GT = fuseMaximums(Ucells_GT,frameRate,t_rise,temporal_pattern,clipping); %%%%%%%
Ucells_GT = randomExcitationAmplitude(Ucells_GT,[1 1.5]); %%%%%%%

[Ucells_GT,cnt] = fuseNearbyExcitations_v3(Ucells_GT,time_fuse,t_rise,frameRate,clipping,temporal_pattern);
fprintf('%d fusions have been carried out\n',cnt);

disp('Activation pattern assigned to cells')

%% DELETE EMPTY GROUPS

empty_cols = [];
for i = 1:size(Ucells_GT,2)
    if ~any(Ucells_GT(:,i))
        empty_cols = [empty_cols i];
    end
end

if ~isempty(empty_cols)
    Ucells_GT(:,empty_cols) = [];
    X(:,:,empty_cols) = [];
    % eliminate those cells that has not been activated during the sequence
    % from their cell clusters, and reassign cell numbers.
    for i = 1:numGroups
        for j = 1:length(empty_cols)
            aux = assemblies{i};
            aux(aux==empty_cols(j)) = [];
            assemblies{i} = aux;
        end
    end
    for k = 1:numGroups
        aux = assemblies{k};
        subtract = find(aux>empty_cols(j));
        for l = 1:length(subtract)
            aux(subtract(l)) = aux(subtract(l))-1;
        end
        assemblies{k} = aux;
    end
end

%% NOISE ADDITION AND VIDEO SEQUENCE CREATION
vecDict = zeros(size(X,1)*size(X,2),size(Ucells_GT,2));
for i = 1:size(Ucells_GT,2)
    vecDict(:,i) = reshape(X(:,:,i),[],1);
end

X_clean = vecDict*Ucells_GT';
X_artifacts = addArtifacts(X_clean,param); %%%%%%%%%%
disp('Adding some noise to make it realistic...')
X_noisy = addingNoiseToimageRA(X_artifacts,RA);

X_new = reshape(X_noisy,[nrows nrows size(Ucells_GT,1)]);
X_new_clean = reshape(X_clean,[nrows nrows size(Ucells_GT,1)]);
toc
disp('Ready!')

X_new = single(X_new);
X_new_clean = single(X_new_clean);
X = single(X);
Ucells_GT = single(Ucells_GT);

