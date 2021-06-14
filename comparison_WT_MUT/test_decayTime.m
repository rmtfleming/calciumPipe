
load('/Volumes/images/SBG/experimentalProgram/Experiments2016/ELM-Exp-46_analysis/analysis/controlplate12/decayTime/bClassAll.mat')

bAll.data = [bClassAll{:}];
% bAll.data = b_val;
bAll.mean = cellfun(@mean, bAll.data);
bAll.std  = cellfun(@std, bAll.data);

% bAll.data   = bAll.data*2.5;
% bAll.mean   = bAll.mean*2.5;
% bAll.std    = bAll.std*2.5;

% find WT and MUT indexes
sigClassAll = cell2mat(sigClass');
colWT = find(sigClassAll(:,end)==1);
colMUT = find(sigClassAll(:,end)==2);

tauWT.mean  = bAll.mean(colWT);
tauMUT.mean = bAll.mean(colMUT);
tauWT.std   = bAll.std(colWT);
tauMUT.std  = bAll.std(colMUT);

% find and remove outliers
potential_outlier = removeOutliers(1:length(tauWT.mean), tauWT.mean, 0, 4);
tauWT.mean(potential_outlier) = [];
tauWT.std(potential_outlier)   =[];

tauWT.mean(isnan(tauWT.mean))=[];
tauWT.std(isnan(tauWT.std))=[];

potential_outlier = removeOutliers(1:length(tauMUT.mean), tauMUT.mean, 0, 4);
tauMUT.mean(potential_outlier)  = [];
tauMUT.std(potential_outlier)   =[];

tauMUT.mean(isnan(tauMUT.mean)) =[];
tauMUT.std(isnan(tauMUT.std))   =[];

% plot distribution
% [fWT,xiWT] = ksdensity(tauWT,'support','positive');
% [fMUT,xiMUT] = ksdensity(tauMUT,'support','positive');
% figure,
% h1 = plot(xiWT,fWT,'Color', [0.85 0.325 0.098]);
% hold on
% h2 = plot(xiMUT,fMUT, 'Color', [0 0.44 0.74]);
% hold off
% alpha(.8)
% set(h1,'LineWidth',2)
% set(h2,'LineWidth',2)
% set(gca,'FontSize', 24,'FontWeight','bold')
% legend([h1 h2],' WT',' LRRK2');
% xlabel('\tau (s)')
% ylabel('Probability')
% 
% dist = KLDivergence_random_distribution(tauWT', tauMUT', 1000, {'\tau'}, 0, pathKLDivPercent);

% dist = KLDivergence_random_distribution(tauWT.mean', tauMUT.mean', 1000,  {'\tau'}, 0, pathKLDivPercent);

