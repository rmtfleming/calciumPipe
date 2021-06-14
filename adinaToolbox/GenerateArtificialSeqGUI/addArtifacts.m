function X_artifacts = addArtifacts(X_clean,param)

fields = fieldnames(param);

for i = 1:length(fields),
eval([fields{i} ' = param.' fields{i} ';']);
end

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
nArtifacts = round(prct_artifacts*nDicts)

%% SAMPLES GENERATION
X = generateSamples2(nArtifacts,nrows,shape,dimension_max,dimension_min); %generate cells arbitrarily allocated
if isempty(X)
    X_artifacts = X_clean;
    return
end

A = zeros(nArtifacts); 
final = 1;

while final,
    for q = 1:nArtifacts,
        for z = 1:nArtifacts,
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
        for q = 1:nArtifacts,
            for z = 1:nArtifacts,
                A(q,z) = sum(sum(X(:,:,q).*X(:,:,z)>0))/sum(sum(X(:,:,q)>0));
            end
        end
%         A(A>=thr_cells) = 1;
%         A(A<thr_cells) = 0;
% 
%         % all the pair of cells which are overlapped have a "1" in the
%         % matrix A
    end
end

h = fspecial('gaussian',[15 15],3);
for i = 1:size(X,3)
    X(:,:,i) = imfilter(X(:,:,i),h);
    %X(:,:,i) = X(:,:,i)./norm(X(:,:,i),'fro');
end

%% TEMPORAL DIMENSION
activationPattern = generateActivationPattern(temporal_pattern,tau/2,t_rise/2,frameRate,clipping);

lengthSequence = length(activationPattern);

U = zeros(nSamples,nArtifacts);
for i = 1:nArtifacts,
    posDelta = randi([1 nSamples-lengthSequence],1);
    U((posDelta:(posDelta+lengthSequence-1)),i) = 1;
    % random distribution of the sequences along the time axis
end

Ucells_GT = assignActivationPatternToCells_v2(temporal_pattern,tau/3,t_rise/3,frameRate,clipping,U,activationPattern);
Ucells_GT = randomExcitationAmplitude(Ucells_GT,[ampl_artifacts_min ampl_artifacts_max]); %%%%

vecDict = zeros(size(X,1)*size(X,2),size(Ucells_GT,2));
for i = 1:size(Ucells_GT,2)
    vecDict(:,i) = reshape(X(:,:,i),[],1);
end

X_art = vecDict*Ucells_GT';
X_artifacts = zeros(size(X_art));

for i = 1:size(X_art,2)
    X_artifacts(:,i) = X_art(:,i)+X_clean(:,i);
end


