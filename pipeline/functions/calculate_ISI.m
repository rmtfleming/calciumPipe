function [signalClass, isiMean] = calculate_ISI(F_tab, fps)

load('Z:\SBG\experimentalProgram\Experiments2016\ELM-Exp-46_analysis\analysis\controlplate12\LRRK2_paper\mat\sigClass.mat')

T = size(F_tab,2);
tvec(1,1) = 0;
for i = 2:T+1
    tvec(1,i-1)=(1/fps)*i; 

end
% find max of each class
classIntervals = cellfun(@(x) max(x(:,1)), sigClass);

% spike detection
for i = 1 : size(F_tab,1)
    [~, events{i}] = findpeaks(F_tab(i,:),'MinPeakHeight',...
        mean(F_tab(i,:)), 'MinPeakDistance',5);
    
%     events{i} = EventDetection2(F_tab(i,:));

    isiMean(i) = mean(diff(tvec(events{i})));
    col = find(classIntervals > isiMean(i));
    if isempty(col)
        signalClass(i) = length(classIntervals);
    else
        signalClass(i) = min(col);
    end
end