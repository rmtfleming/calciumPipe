function [X_new,X_new_clean,assemblies,X,Ucells_GT] = generateArtificialSequence_v5(varargin)

% function "generateArtificialSequences" returns:
% 
%   X_new: nrows x nrows x nSamples matrix. It contains the noisy video
% sequence of the neurons.
%   X_new_clean: nrows x nrows x nSamples matrix. It contains the noisy video
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
%%%%


%% USER PARAMETERS SECTION
tic

%%%%%%%%%%%%%%%%%%%%%%%%%% USER PARAMETERS: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[param_out] = initial_parameter_ArtificialSeq();
fields = param_out.List;
param = param_out.all;

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
    
    if ~isstruct(varargin{3})
        param = set_parameters_ArtificialSeq(param,varargin{3});
    else
        param = varargin{3};
    end    
    
    for i = 1:length(fields),
        eval([fields{i} ' = param.' fields{i} ';']);
    end
    
    nDicts = varargin{1};
    seqDuration = varargin{2};    
    if ~isa(nDicts,'double') || ~isa(seqDuration,'double')
        msgbox('Incorrect argument types!','Error','error','modal')
        return
    end
    
elseif nargin > 3
    msgbox('Too many input arguments!','Error','error','modal')
    return
end

% DEFAULT SETTINGS:
% 
% nDicts = 10; % number of cells
% nrows = 512; % image size
% frameRate = 30; % frame rate in fps
% seqDuration = 15; % sequence duration in seconds
% t_rise = 300; % neuron activation rise time in milliseconds
% tau = 0.7; % decay parameter of the cells --> exp(-t/tau)
% thr_cells = 0.5; % threshold of overlapping between cells
% thr_seq = 0.2; % threshold of overlapping between time sequences
% numGroups = 5; % number of cell clusters. MUST BE LOWER OR EQUAL THAN "NDICTS"
% percentMult = 0.4; % the percentage of cells that will be multiple assigned
% maxMult = 3; % maximum number of multiplicities per cell
% shape = 'gaussian'; % shape of the cells. Options: 'octogonal', 'gaussian', 'real'
% temporal_pattern = 'exp_exp'; % temporal activation pattern. Options: 'slope', 'gauss_exp', 'exp_exp'
% roxmax = 0.05; % roxmax determines the maximum x-size of cells w.r.t. the image size. (0<ro<0.2. optimal 0.05)
% roxmin = 0.04; % roxmin determines the minimum x-size of cells w.r.t. the image size
% roymax = 0.05; % roymax determines the maximum y-size of cells w.r.t. the image size. (0<ro<0.2. optimal 0.05)
% roymin = 0.04; % roymin determines the minimum y-size of cells w.r.t. the image size
% RA = 10; % Relative Amplitude of the sequence generated
% clipping = 0; % clipping percentage w.r.t. the maximum value. The values of the activation pattern below max*clipping will be clipped to 0.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% WARNING SECTION
dimension_max = [ceil(nrows*roxmax) ceil(nrows*roymax)]; % maximum dimension [x y] of the cells in pixels
dimension_min = [ceil(nrows*roxmin) ceil(nrows*roymin)]; % minimum dimension [x y] of the cells in pixels
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

lenDec = 3*tau; % lenght (in seconds) of the decay period
lenDecFr = round(lenDec*frameRate); %length in frames of the decay period
lenRise = round((t_rise/1000)*frameRate); %length (in frames) of the rise period

lenSeq = lenDecFr+lenRise;
areaCovered = lenSeq*(1-thr_seq)*numGroups + thr_seq*lenSeq;

if areaCovered > nSamples
    errordlg('The temporal parameters do not match. Check whether the number of clusters and its active period duration fits into the temporal samples!','Error','modal')
    error('Temporal parameters do not match');
elseif 2*areaCovered > nSamples
    msgbox('Your temporal parameters just fit. It migth take several time to converge...','Warning','warn','modal')
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
        A(A>0) = 1;

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

activationPattern = generateActivationPattern(temporal_pattern,tau,t_rise,frameRate,clipping);
disp('Activation pattern generated')

lengthSequence = length(activationPattern);


U = zeros(nSamples,numGroups);
for i = 1:numGroups,
    posDelta = randi([1 nSamples-lengthSequence],1);
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
            posDelta = randi([lengthSequence nSamples-lengthSequence],1);
            U(:,listDel(q)) = 0;
            U((posDelta:(posDelta+lengthSequence-1)),listDel(q)) = 1;
        end
    else
        final = 0;
    end
end

disp('Sequence temporal position determined')

%% ACTIVATION PATTERN ASSIGNMENT

Ucells_temporalPattern = zeros(size(U,1),size(X,3));
for i = 1:size(U,1)
    groups = find(U(i,:)>0);
    if isempty(groups)
        Ucells_temporalPattern(i,:)=0;
    else
        for j = 1:length(groups)
            cells = assemblies{groups(j)};
            for k = 1:length(cells)
                Ucells_temporalPattern(i,cells(k))=U(i,groups(j));
            end
        end
    end
end

Ucells_GT = assignActivationPatternToCells(temporal_pattern,t_rise,frameRate,clipping,Ucells_temporalPattern,activationPattern);

disp('Activation pattern assigned to cells')

%% NOISE ADDITION AND VIDEO SEQUENCE CREATION

% % images
% imNeuron = zeros(size(X,1),size(X,2),numGroups);
% vecNeuron = zeros(size(X,1)*size(X,2),numGroups);
% for i = 1:numGroups,
%     idx = assemblies{i};
%     % imNeuron = 3D matrix which contains #numGroups clusters. Each
%     % image contains the whole i'th cluster
%     imNeuron(:,:,i) = max(X(:,:,idx),[],3);
%     % vecNeuron = 2D matrix which contains the clusters reordered by
%     % columns
%     vecNeuron(:,i) = reshape(imNeuron(:,:,i),[],1);
% end

vecDict = zeros(size(X,1)*size(X,2),nDicts);
for i = 1:nDicts
    vecDict(:,i) = reshape(X(:,:,i),[],1);
end

X_clean = vecDict*Ucells_GT';
disp('Adding some noise to make it realistic...')
X_noisy = addingNoiseToimageRA(X_clean,RA);

X_new = reshape(X_noisy,[nrows nrows nSamples]);
X_new_clean = reshape(X_clean,[nrows nrows nSamples]);
toc
disp('Ready!')


