function [dF_F Fbase] = deltaF_F(time_traces, halfwindow, quantile, step, sigma)

% This function calculates the relative changes in fluorescence intensity
% (F-F0)/F0 for each time trace

% INPUTS:
%     traces     : a matrix n*m of fluorescence traces where each row is a
%                   time series
%     halfwindow : half of the size of the sliding windows
%     quantile   : which quantile to use for f0 computation
%     step      : to speed up f0 estimation, f0 is only computed every this many frames


% Created by: Siham Hachi, 01/10/2016


    h_3(1,1,:) = fspecial('gaussian',[1 5*sigma],sigma);
    Wgausstime = imfilter(time_traces, h_3,'corr','symmetric','same');
    reshape(double(Wgausstime),[],size(time_traces,2));
    [Fbase] = fastExtractingQuantilesSeq(reshape(double(Wgausstime),[],size(time_traces,2)),halfwindow,quantile,step);
    Fbase_smooth = imfilter(Fbase,h_3,'corr','symmetric','same');
    dF_F = (time_traces'-Fbase_smooth)./Fbase_smooth;

end