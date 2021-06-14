pathK7 = '/Users/siham.hachi/Desktop/all_CaSA_analysis/TavFilter10/k7/';
pathIMM = '/Users/siham.hachi/Desktop/all_CaSA_analysis/TavFilter10/IMM/';
pathRes = '/Users/siham.hachi/Desktop/all_CaSA_analysis/TavFilter10/results/';
listK7 = dir(pathK7);
listIMM = dir(pathIMM);

%% K7

j=1;
for i=4 : length(listK7)-1
    peakNumK7{j} = xlsread([pathK7 listK7(i).name],'J:J');    
    [num TavTxtK7{j}] = xlsread([pathK7 listK7(i).name],'N:N');    
    [num TavSDTxtK7{j}]= xlsread([pathK7 listK7(i).name],'O:O');
    [num SWTxtK7{j}] = xlsread([pathK7 listK7(i).name],'Q:Q');   
    [num SWSDTxtK7{j}] = xlsread([pathK7 listK7(i).name],'R:R');
    [num ampTxtK7{j}] = xlsread([pathK7 listK7(i).name],'S:S');
    [num energyTxtK7{j}] = xlsread([pathK7 listK7(i).name],'U:U');
    [num powerTxtK7{j}] = xlsread([pathK7 listK7(i).name],'V:V');
    [num spikeAreaMeanTxtK7{j}] = xlsread([pathK7 listK7(i).name],'Y:Y');
    [num spikeAreaSDTxtK7{j}] = xlsread([pathK7 listK7(i).name],'Z:Z');
    [num timeToPeakTxtK7{j}] = xlsread([pathK7 listK7(i).name],'AC:AC');
    [num AVCalciumReleasingRateTxtK7{j}] = xlsread([pathK7 listK7(i).name],'AD:AD');
    [num AVCalciumRemovingRateTxtK7{j}] = xlsread([pathK7 listK7(i).name],'AE:AE');
   
    j=j+1;
end

% convert strings to numbers
TavTxtK7 = cellfun(@str2double,TavTxtK7, 'UniformOutput', false);
TavSDTxtK7 = cellfun(@str2double,TavSDTxtK7, 'UniformOutput', false);
SWTxtK7 = cellfun(@str2double,SWTxtK7, 'UniformOutput', false);
SWSDTxtK7 = cellfun(@str2double,SWSDTxtK7, 'UniformOutput', false);
ampTxtK7 = cellfun(@str2double,ampTxtK7, 'UniformOutput', false);
energyTxtK7 = cellfun(@str2double,energyTxtK7, 'UniformOutput', false);
powerTxtK7 = cellfun(@str2double,powerTxtK7, 'UniformOutput', false);
spikeAreaMeanTxtK7 = cellfun(@str2double,spikeAreaMeanTxtK7, 'UniformOutput', false);
spikeAreaSDTxtK7 = cellfun(@str2double,spikeAreaSDTxtK7, 'UniformOutput', false);
timeToPeakTxtK7 = cellfun(@str2double,timeToPeakTxtK7, 'UniformOutput', false);
AVCalciumReleasingRateTxtK7 = cellfun(@str2double,AVCalciumReleasingRateTxtK7, 'UniformOutput', false);
AVCalciumRemovingRateTxtK7 = cellfun(@str2double,AVCalciumRemovingRateTxtK7, 'UniformOutput', false);

% concatenate cells from each dataset in one vector
peakNumK7 = vertcat(peakNumK7{:});
TavTxtK7 = vertcat(TavTxtK7{:}); TavTxtK7(isnan(TavTxtK7(:,1)),:)=[];
TavSDTxtK7 = vertcat(TavSDTxtK7{:}); TavSDTxtK7(isnan(TavSDTxtK7(:,1)),:)=[];
SWTxtK7 = vertcat(SWTxtK7{:}); SWTxtK7(isnan(SWTxtK7(:,1)),:)=[];
SWSDTxtK7 = vertcat(SWSDTxtK7{:}); SWSDTxtK7(isnan(SWSDTxtK7(:,1)),:)=[];
ampTxtK7 = vertcat(ampTxtK7{:}); ampTxtK7(isnan(ampTxtK7(:,1)),:)=[];
energyTxtK7 = vertcat(energyTxtK7{:}); energyTxtK7(isnan(energyTxtK7(:,1)),:)=[];
powerTxtK7 = vertcat(powerTxtK7{:}); powerTxtK7(isnan(powerTxtK7(:,1)),:)=[];
spikeAreaMeanTxtK7 = vertcat(spikeAreaMeanTxtK7{:}); spikeAreaMeanTxtK7(isnan(spikeAreaMeanTxtK7(:,1)),:)=[];
spikeAreaSDTxtK7 = vertcat(spikeAreaSDTxtK7{:}); spikeAreaSDTxtK7(isnan(spikeAreaSDTxtK7(:,1)),:)=[];
timeToPeakTxtK7 = vertcat(timeToPeakTxtK7{:}); timeToPeakTxtK7(isnan(timeToPeakTxtK7(:,1)),:)=[];
AVCalciumReleasingRateTxtK7 = vertcat(AVCalciumReleasingRateTxtK7{:}); AVCalciumReleasingRateTxtK7(isnan(AVCalciumReleasingRateTxtK7(:,1)),:)=[];
AVCalciumRemovingRateTxtK7 = vertcat(AVCalciumRemovingRateTxtK7{:}); AVCalciumRemovingRateTxtK7(isnan(AVCalciumRemovingRateTxtK7(:,1)),:)=[];

%% IMM

j=1;
for i=3 : length(listIMM)
    peakNumIMM{j} = xlsread([pathIMM listIMM(i).name],'J:J');    
    [num TavTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'N:N');    
    [num TavSDTxtIMM{j}]= xlsread([pathIMM listIMM(i).name],'O:O');
    [num SWTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'Q:Q');    
    [num SWSDTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'R:R');
    [num ampTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'S:S');
    [num energyTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'U:U');
    [num powerTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'V:V');
    [num spikeAreaMeanTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'Y:Y');
    [num spikeAreaSDTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'Z:Z');
    [num timeToPeakTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'AC:AC');
    [num AVCalciumReleasingRateTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'AD:AD');
    [num AVCalciumRemovingRateTxtIMM{j}] = xlsread([pathIMM listIMM(i).name],'AE:AE');
   
    j=j+1;
end

TavTxtIMM = cellfun(@str2double,TavTxtIMM, 'UniformOutput', false);
TavSDTxtIMM = cellfun(@str2double,TavSDTxtIMM, 'UniformOutput', false);
SWTxtIMM = cellfun(@str2double,SWTxtIMM, 'UniformOutput', false);
SWSDTxtIMM = cellfun(@str2double,SWSDTxtIMM, 'UniformOutput', false);
ampTxtIMM = cellfun(@str2double,ampTxtIMM, 'UniformOutput', false);
energyTxtIMM = cellfun(@str2double,energyTxtIMM, 'UniformOutput', false);
powerTxtIMM = cellfun(@str2double,powerTxtIMM, 'UniformOutput', false);
spikeAreaMeanTxtIMM = cellfun(@str2double,spikeAreaMeanTxtIMM, 'UniformOutput', false);
spikeAreaSDTxtIMM = cellfun(@str2double,spikeAreaSDTxtIMM, 'UniformOutput', false);
timeToPeakTxtIMM = cellfun(@str2double,timeToPeakTxtIMM, 'UniformOutput', false);
AVCalciumReleasingRateTxtIMM = cellfun(@str2double,AVCalciumReleasingRateTxtIMM, 'UniformOutput', false);
AVCalciumRemovingRateTxtIMM = cellfun(@str2double,AVCalciumRemovingRateTxtIMM, 'UniformOutput', false);

peakNumIMM = vertcat(peakNumIMM{:});
TavTxtIMM = vertcat(TavTxtIMM{:}); TavTxtIMM(isnan(TavTxtIMM(:,1)),:)=[];
TavSDTxtIMM = vertcat(TavSDTxtIMM{:}); TavSDTxtIMM(isnan(TavSDTxtIMM(:,1)),:)=[];
SWTxtIMM = vertcat(SWTxtIMM{:}); SWTxtIMM(isnan(SWTxtIMM(:,1)),:)=[];
SWSDTxtIMM = vertcat(SWSDTxtIMM{:}); SWSDTxtIMM(isnan(SWSDTxtIMM(:,1)),:)=[];
ampTxtIMM = vertcat(ampTxtIMM{:}); ampTxtIMM(isnan(ampTxtIMM(:,1)),:)=[];
energyTxtIMM = vertcat(energyTxtIMM{:}); energyTxtIMM(isnan(energyTxtIMM(:,1)),:)=[];
powerTxtIMM = vertcat(powerTxtIMM{:}); powerTxtIMM(isnan(powerTxtIMM(:,1)),:)=[];
spikeAreaMeanTxtIMM = vertcat(spikeAreaMeanTxtIMM{:}); spikeAreaMeanTxtIMM(isnan(spikeAreaMeanTxtIMM(:,1)),:)=[];
spikeAreaSDTxtIMM = vertcat(spikeAreaSDTxtIMM{:}); spikeAreaSDTxtIMM(isnan(spikeAreaSDTxtIMM(:,1)),:)=[];
timeToPeakTxtIMM = vertcat(timeToPeakTxtIMM{:}); timeToPeakTxtIMM(isnan(timeToPeakTxtIMM(:,1)),:)=[];
AVCalciumReleasingRateTxtIMM = vertcat(AVCalciumReleasingRateTxtIMM{:}); AVCalciumReleasingRateTxtIMM(isnan(AVCalciumReleasingRateTxtIMM(:,1)),:)=[];
AVCalciumRemovingRateTxtIMM = vertcat(AVCalciumRemovingRateTxtIMM{:}); AVCalciumRemovingRateTxtIMM(isnan(AVCalciumRemovingRateTxtIMM(:,1)),:)=[];
   
%% Plot results

xlabels = {'K7'; 'IMM'};

% Peak number
figure
bar(1:2,[mean(peakNumK7), mean(peakNumIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(peakNumK7), mean(peakNumIMM)], ...
    [std(peakNumK7), std(peakNumIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('Peak Number')
grid on
figname=[pathRes 'peakNumber'];
print('-dtiff',figname);

% Tav
figure
bar(1:2,[mean(TavTxtK7), mean(TavTxtIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(TavTxtK7), mean(TavTxtIMM)], ...
    [std(TavTxtK7), std(TavTxtIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('Tav')
grid on
figname=[pathRes 'Tav'];
print('-dtiff',figname);

% Spike Width
figure
bar(1:2,[mean(SWTxtK7), mean(SWTxtIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(SWTxtK7), mean(SWTxtIMM)], ...
    [std(SWTxtK7), std(SWTxtIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('Spike Width')
grid on
figname=[pathRes 'SpikeWidth'];
print('-dtiff',figname);

% Mean Spike Area
figure
bar(1:2,[mean(spikeAreaMeanTxtK7), mean(spikeAreaMeanTxtIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(spikeAreaMeanTxtK7), mean(spikeAreaMeanTxtIMM)], ...
    [std(spikeAreaMeanTxtK7), std(spikeAreaMeanTxtIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('Mean Spike Area')
grid on
figname=[pathRes 'MeanSpikeArea'];
print('-dtiff',figname);

% Time To Peak
figure
bar(1:2,[mean(timeToPeakTxtK7), mean(timeToPeakTxtIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(timeToPeakTxtK7), mean(timeToPeakTxtIMM)], ...
    [std(timeToPeakTxtK7), std(timeToPeakTxtIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('Time To Peak (s)')
grid on
figname=[pathRes 'TimeToPeak'];
print('-dtiff',figname);

% Amplitude
figure
bar(1:2,[mean(ampTxtK7), mean(ampTxtIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(ampTxtK7), mean(ampTxtIMM)], ...
    [std(ampTxtK7), std(ampTxtIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('Amplitude')
grid on
figname=[pathRes 'Amplitude'];
print('-dtiff',figname);

% Energy
figure
bar(1:2,[mean(energyTxtK7), mean(energyTxtIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(energyTxtK7), mean(energyTxtIMM)], ...
    [std(energyTxtK7), std(energyTxtIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('Energy')
grid on
figname=[pathRes 'Energy'];
print('-dtiff',figname);

% Power
figure
bar(1:2,[mean(powerTxtK7), mean(powerTxtIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(powerTxtK7), mean(powerTxtIMM)], ...
    [std(powerTxtK7), std(powerTxtIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('Power')
grid on
figname=[pathRes 'Power'];
print('-dtiff',figname);

% AV Calcium Releasing Rate
figure
bar(1:2,[mean(AVCalciumReleasingRateTxtK7), mean(AVCalciumReleasingRateTxtIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(AVCalciumReleasingRateTxtK7), mean(AVCalciumReleasingRateTxtIMM)], ...
    [std(AVCalciumReleasingRateTxtK7), std(AVCalciumReleasingRateTxtIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('AVCalciumReleasingRate')
grid on
figname=[pathRes 'AVCalciumReleasingRate'];
print('-dtiff',figname);

% AV Calcium Removing Rate
figure
bar(1:2,[mean(AVCalciumRemovingRateTxtK7), mean(AVCalciumRemovingRateTxtIMM)], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mean(AVCalciumRemovingRateTxtK7), mean(AVCalciumRemovingRateTxtIMM)], ...
    [std(AVCalciumRemovingRateTxtK7), std(AVCalciumRemovingRateTxtIMM)],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])
ylabel('AVCalciumRemovingRate')
grid on
figname=[pathRes 'AVCalciumRemovingRate'];
print('-dtiff',figname);

%% Plot points

% peak number
figure,
for i=1 : length(peakNumK7)
    plot(i,peakNumK7(i),'*r'), hold on
end
hold on,
for i=1 : length(peakNumIMM)
    plot(i,peakNumIMM(i),'*b'), hold on
end
set(gca, 'FontSize', 14,'FontWeight','bold', 'YTickLabel', [])    
ylabel('Peak Number')
grid on

% Tav
TavTxtK7 = sort(TavTxtK7);
TavTxtIMM = sort(TavTxtIMM);
figure,
for i=1 : length(TavTxtK7)
    plot(i,TavTxtK7(i),'*r'), hold on
end
hold on,
for i=1 : length(TavTxtIMM)
    plot(i,TavTxtIMM(i),'*b'), hold on
end

% spike Width
SWTxtK7 = sort(SWTxtK7);
SWTxtIMM = sort(SWTxtIMM);
figure,
for i=1 : length(SWTxtK7)
    plot(i,SWTxtK7(i),'*r'), hold on
end
hold on,
for i=1 : length(SWTxtIMM)
    plot(i,SWTxtIMM(i),'*b'), hold on
end

% spike area
figure,
for i=1 : length(spikeAreaMeanTxtK7)
    plot(i,spikeAreaMeanTxtK7(i),'*r'), hold on
end
hold on,
for i=1 : length(spikeAreaMeanTxtIMM)
    plot(i,spikeAreaMeanTxtIMM(i),'*b'), hold on
end

% time To Peak
figure,
for i=1 : length(timeToPeakTxtK7)
    plot(i,timeToPeakTxtK7(i),'*r'), hold on
end
hold on,
for i=1 : length(timeToPeakTxtIMM)
    plot(i,timeToPeakTxtIMM(i),'*b'), hold on
end

% amplitude
figure,
for i=1 : length(ampTxtK7)
    plot(i,ampTxtK7(i),'*r'), hold on
end
hold on,
for i=1 : length(ampTxtIMM)
    plot(i,ampTxtIMM(i),'*b'), hold on
end

% energy
figure,
for i=1 : length(energyTxtK7)
    plot(i,energyTxtK7(i),'*r'), hold on
end
hold on,
for i=1 : length(energyTxtIMM)
    plot(i,energyTxtIMM(i),'*b'), hold on
end

% power 
figure,
for i=1 : length(powerTxtK7)
    plot(i,powerTxtK7(i),'*r'), hold on
end
hold on,
for i=1 : length(powerTxtIMM)
    plot(i,powerTxtIMM(i),'*b'), hold on
end

% AV Calcium Releasing Rate
figure,
for i=1 : length(AVCalciumReleasingRateTxtK7)
    plot(i,AVCalciumReleasingRateTxtK7(i),'*r'), hold on
end
hold on,
for i=1 : length(AVCalciumReleasingRateTxtIMM)
    plot(i,AVCalciumReleasingRateTxtIMM(i),'*b'), hold on
end

% 
figure,
for i=1 : length(AVCalciumRemovingRateTxtK7)
    plot(i,AVCalciumRemovingRateTxtK7(i),'*r'), hold on
end
hold on,
for i=1 : length(AVCalciumRemovingRateTxtIMM)
    plot(i,AVCalciumRemovingRateTxtIMM(i),'*b'), hold on
end

%%

% peak number
mk7 = mean([mean(peakNumK7{1}),mean([mean(peakNumK7{2}),mean(peakNumK7{3}),mean(peakNumK7{4})])]);
sk7 = std([mean(peakNumK7{1}),mean([mean(peakNumK7{2}),mean(peakNumK7{3}),mean(peakNumK7{4})])]);

mIMM = mean([mean(peakNumIMM{1}),mean([mean(peakNumIMM{2}),mean(peakNumIMM{3}),mean(peakNumIMM{4})])]);
sIMM = std([mean(peakNumIMM{1}),mean([mean(peakNumIMM{2}),mean(peakNumIMM{3}),mean(peakNumIMM{4})])]);


figure
bar(1:2,[mk7, mIMM], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mk7, mIMM], ...
    [sk7, sIMM],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold')
ylabel('peak number')
grid on


% amp
mk7 = mean([mean(peakNumK7{1}),mean([mean(peakNumK7{2}),mean(peakNumK7{3}),mean(peakNumK7{4})])]);
sk7 = std([mean(peakNumK7{1}),mean([mean(peakNumK7{2}),mean(peakNumK7{3}),mean(peakNumK7{4})])]);

mIMM = mean([mean(peakNumIMM{1}),mean([mean(peakNumIMM{2}),mean(peakNumIMM{3}),mean(peakNumIMM{4})])]);
sIMM = std([mean(peakNumIMM{1}),mean([mean(peakNumIMM{2}),mean(peakNumIMM{3}),mean(peakNumIMM{4})])]);


figure
bar(1:2,[mk7, mIMM], 0.5, 'FaceColor', [.2 .5 .7])
hold on, errorbar(1:2,[mk7, mIMM], ...
    [sk7, sIMM],'.', 'LineWidth', 3, 'Color', 'k')
set(gca,'xticklabel',xlabels, 'FontSize', 14,'FontWeight','bold')
ylabel('peak number')
grid on