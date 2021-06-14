% clear,

tic,

load('sigClassTracesAll.mat')

for ii = 3 : length(sigClassTracesAll)
    
    b_val={};
    for jj = 1 : size(sigClassTracesAll{1,ii},1)
        
        F = sigClassTracesAll{1,ii}(jj,:);
        F = nonzeros(F);
        [dF_F Fbase] = deltaF_F(F', 10, 0.2, 2, 3); 
        dF_F = dF_F';

        for i = 1 : size(dF_F,1)
            events{i} = EventDetection(dF_F(i,:));
        end
        
        if length(events{:}) > 1
            s.dF_cell       = dF_F;
            s.Spikes_cell   = events;
            s.fps           = 10;
            [A,G,fit,gof,rise_time,fall_time,CV] = getTransients(s);
            fit = fit(~cellfun('isempty',fit));
            b_val{jj} = 1./cellfun(@(x) x.b, fit);

            X = [1 : length(b_val{jj})];
            Y = b_val{jj};
            potential_outlier = removeOutliers(X, Y, 0, 1);
            
        end
        
    end
    
    bClassAll{ii} = b_val;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc,