function [jF, jFN, jFIm] = find_trace_files(files)

% find the time traces files for WT and MUT
jF      = [];
jFN     = [];
jFIm    = [];
for i=1:length(files)
    if isempty(strfind(files{i},'old'))...
            && (~isempty(strfind(files{i},'K7')) || ~isempty(strfind(files{i},'K7M')))
        if ~isempty(strfind(files{i},'time_traces'))
            jF = [jF i];
        elseif  ~isempty(strfind(files{i},'bits_cNeurons.'))
            jFN = [jFN i];
            elseif  ~isempty(strfind(files{i},'8bits.'))
                jFIm = [jFIm i];
            end
        end

    end
end