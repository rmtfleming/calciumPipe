function sigs = load_signal_param(pathh)

% This function loads data from an excel sheet
%
% INPUT:
%     pathh   : folder containing the data in excel sheets
% OUTPUT:
%     sigs:

fileList = dir(pathh);

j=1;
for i=4 : length(fileList)
    sigs(j).fileName = fileList(i).name;       

    [num sigs(j).isi] = xlsread([pathh fileList(i).name],'I:I');
    sigs(j).isi = str2double(sigs(j).isi);
    sigs(j).isi(1) = [];
    zeroIdx = find(sigs(j).isi==0);
    sigs(j).isi(zeroIdx) = [];

    sigs(j).sigID = xlsread([pathh fileList(i).name],'A:A');
    sigs(j).sigID(zeroIdx) = [];

    [num sigs(j).amp] = xlsread([pathh fileList(i).name],'G:G');
    sigs(j).amp = str2double(sigs(j).amp);
    sigs(j).amp(1) = [];
    sigs(j).amp(zeroIdx) = [];

    [num sigs(j).SWidth] = xlsread([pathh fileList(i).name],'F:F');
    sigs(j).SWidth = str2double(sigs(j).SWidth);
    sigs(j).SWidth(1) = [];
    sigs(j).SWidth(zeroIdx) = [];

    [num sigs(j).SArea] = xlsread([pathh fileList(i).name],'J:J');
    sigs(j).SArea = str2double(sigs(j).SArea);
    sigs(j).SArea(1) = [];
    sigs(j).SArea(zeroIdx) = [];

    [num sigs(j).timeToPeak] = xlsread([pathh fileList(i).name],'N:N');
    sigs(j).timeToPeak = str2double(sigs(j).timeToPeak);
    sigs(j).timeToPeak(1) = [];
    sigs(j).timeToPeak(zeroIdx) = [];

    [num sigs(j).releasingRate] = xlsread([pathh fileList(i).name],'O:O');
    sigs(j).releasingRate = str2double(sigs(j).releasingRate);
    sigs(j).releasingRate(1) = [];
    sigs(j).releasingRate(zeroIdx) = [];

    [num sigs(j).removingRate] = xlsread([pathh fileList(i).name],'P:P');
    sigs(j).removingRate = str2double(sigs(j).removingRate);
    sigs(j).removingRate(1) = [];
    sigs(j).removingRate(zeroIdx) = [];

    j=j+1;
end

