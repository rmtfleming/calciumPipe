% clear,

tic,

load('sigClassTracesAll.mat')

sigClassTracesAll = cell2mat(sigClassTracesAll');
b_val = cell(1,size(sigClassTracesAll,1));

for jj = 1 : size(sigClassTracesAll,1)

    F = sigClassTracesAll(jj,:);
    F = nonzeros(F);
    [dF_F, ~] = deltaF_F(F', 10, 0.2, 2, 3); 
    dF_F = dF_F';

    for i = 1 : size(dF_F,1)
        events{i} = EventDetection2(dF_F(i,:));
    end

    if length(events{:}) > 1
        s.dF_cell       = dF_F;
        s.Spikes_cell   = events;
        s.fps           = 2;
        [A,G,fit,gof,rise_time,fall_time,CV] = getTransients(s);
        fit = fit(~cellfun('isempty',fit));
        b_val{jj} = 1./cellfun(@(x) x.b, fit);

        X = [1 : length(b_val{jj})];
        Y = b_val{jj};
        stats = regstats(Y,X,'linear');

    % if Cook's Distance > n/4 is a typical treshold that is used to suggest
    % the presence of an outlier
        potential_outlier = stats.cookd > 0.3/length(X);
%         figure,
%         h = scatter(X,Y, 'r.');
%         hold on        
%         h = scatter(X(potential_outlier),Y(potential_outlier), 'b.');
%         grid on
% 
%         
    b_val{jj}(potential_outlier) = [];
    end

end

% bClassAll{ii} = b_val;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc,