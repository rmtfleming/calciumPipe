
pathh = 'Z:\SBG\experimentalProgram\Experiments2016\ELM-Exp-46_analysis\analysis\controlplate12\';
files = dirrec(pathh,'.xls');

j=1;
for i=1:length(files)
    if strfind(files{i},'_SigalProfile.xls')
        fileIdx(j) = i;
        j=j+1;
    end
end

%%
for i=1:length(fileIdx)    
    load(files{fileIdx(i)});
    T = array2table(time_traces');
    time = [1:size(time_traces,2)]';
    T1 = table(time);
    T = [T1 T];
    [pathstr,name,ext] = fileparts(files{fileIdx(i)}); 
    writetable(T,[pathstr filesep name '.xls']); 
    writetable(T,[pathstr filesep name '.xlsx']);
end


%% change fileNames
for i = 1:length(fileIdx)
    splitt =  strsplit(files{fileIdx(i)},filesep);
    newFileName = [splitt{end-3} '_' splitt{end}];
    [pathstr,name,ext] = fileparts(files{fileIdx(i)});
    movefile(files{fileIdx(i)}, [pathstr filesep newFileName]);
end

%% Find time_traces.mat for K7 and K7M

p = '/Volumes/images/SBG/experimentalProgram/Experiments2016/ELM-Exp-46_analysis/analysis/controlplate12/';
files = dirrec(p,'.mat');

j=1;
for i=1:length(files)
    if ~isempty(strfind(files{i},'time_traces'))...
            && (~isempty(strfind(files{i},'K7')) || ~isempty(strfind(files{i},'K7M')~=[]))
        fileIdx(j) = i;
        j=j+1;
    end
end
