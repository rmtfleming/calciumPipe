function aci
% Activity correlation imaging is a tool for the visualization of neuronal
% morphology based on time-lapsed imaging data (Junek et al., Biophys J 2009).
% The time traces extracted from “reference regions” (e.g. neuronal somata or
% neuropile structures, such as glomeruli) are used to generate correlation
% maps by calculating the correlation coefficient between the reference trace
% and the time course of each pixel (or voxel) of the observed volume. 
% This algorithm has been implemented in a Matlab function for convenient use. 
% Note that PC systems with memory limitations (in particular older 32-bit 
% systems with little memory) might run out of memory. In this case try to 
% reduce the size of your data set (e.g. by reducing the number of time points 
% or spatially binning the data). Read the documentation for further
% information about the use of this function.
%
% Comments and questions to: sjunek1@gwdg.de
%
% Copyright S. Junek, 2009
% _________________________________________________________________________


scrsz = get(0,'ScreenSize');
guiPos = [100 scrsz(4)-150 450 100]; 
mainGui = figure('Position',guiPos,'MenuBar','none','Name','ACI - Main Menu',...
        'NumberTitle','off','Resize','off','CloseRequestFcn',@closeAll);
% Button file io
uicontrol('Parent',mainGui,'Style','pushbutton','Position',[35 30 180 50],'String','1. Prepare Data',...
    'FontSize',12,'FontWeight','bold','Callback',@selectData);
% Button Select Rois
uicontrol('Parent',mainGui,'Style','pushbutton','Position',[235 30 180 50],'String','2. Select ROIs',...
    'FontSize',12,'FontWeight','bold','Callback',@selectRoi);
% Button calculate corr Maps
% uicontrol('Parent',mainGui,'Style','pushbutton','Position',[390 30 180 50],'String','3. Correlation Maps',...
%     'FontSize',12,'FontWeight','bold','Callback',@corrMap);
% Status text
txt_status = uicontrol('Parent',mainGui, 'Style','text','Position',[10 3 400 20], 'String','',...
    'BackgroundColor',get(mainGui, 'Color'),'HorizontalAlignment','left');

% Variables
d = [];
oldPath = cd;
varMap = [];

%% Select and prepare data
    function selectData(varargin)
        workspaceVar = evalin('base','whos');
        if isempty(workspaceVar)
           warndlg('The workspace is empty. Load data set as a single variable ([x,y,t] or [x,y,z,t])into the Matlab workspace!');
           return
        end
        for i=1:numel(workspaceVar)
            s{i} = workspaceVar(i).name;  %#ok
        end
        pos_fig_var = [100 100 200 scrsz(4)-350];
        fig_var = figure('Position',pos_fig_var,'MenuBar','none','Name','ACI - Select Data',...
        'NumberTitle','off','Resize','off','CloseRequestFcn',@closeGetData);
        lb_var = uicontrol('Style','listbox','Position',[10 50 pos_fig_var(3)-20 pos_fig_var(4)-60],...
            'String',s,'Value',1);
        uicontrol('Style','pushbutton','Position',[10 10 pos_fig_var(3)/2-15 20],'String','OK',...
            'FontWeight','bold','Callback',@getData);
        uicontrol('Style','pushbutton','Position',[5+pos_fig_var(3)/2 10 pos_fig_var(3)/2-15 20],'String','Cancel',...
            'FontWeight','bold','Callback',@closeGetData);
        function getData(varargin)
            s = get(lb_var, 'String');
            d = evalin('base',s{get(lb_var,'Value')});
            if ~isnumeric(d)
                d = [];
                warndlg('Selected variable is not a numeric array. Select other variable!','Wrong data type');
            elseif numel(size(d))==2 || numel(size(d))>4 
                d = [];
                warndlg('Selected variable has wrong dimension. Select matrix with dimension 3 or 4!','Wrong dimension');
            else
                if numel(size(d)) == 3
                    d = permute(d,[1 2 4 3]);
                end
                closeGetData
%                  [fn pn] = uigetfile;
%                  d = readCustomTif([pn fn]);
                pathName = uigetdir(cd,'Select or create directory for this data set!');
                if ~pathName
                    return
                end
                cd(pathName);
                set(txt_status, 'String','Wait while data are prepared for further analysis!');
                drawnow;
                prepareData
                set(txt_status, 'String','Done.');
            end
        end
        function prepareData 
            dim = size(d);
%             varMap = mean(d,4);
            subDivide = 1;  %subdivide files for bleach correction
            bleachCorrLength = dim(4) / subDivide;
            % Bleach correction
            A=[ones(bleachCorrLength,1),(1:bleachCorrLength)'];
            [A,s]=svd(A,0);
            d = single(reshape(d,[],dim(3),dim(4)));
            for ii=1:dim(3)
                for k=1:subDivide
                    d(:,ii,(k-1)*bleachCorrLength+1:k*bleachCorrLength) = d(:,ii,(k-1)*bleachCorrLength+1:k*bleachCorrLength) - permute((squeeze(d(:,ii,(k-1)*bleachCorrLength+1:k*bleachCorrLength))*A)*A',[1 3 2]);
                end
            end
            d = reshape(d,dim);
            varMap = std(d,[],4);
            for ii=1:dim(3)
%                 load(['slice',num2str(i,'%02d')]);
                s = squeeze(d(:,:,ii,:));
                s = reshape(s, [dim(1)*dim(2),dim(4)]);
                s = s./repmat((sum(s.^2,2)).^0.5,[1 dim(4)]);
                d(:,:,ii,:) = reshape(s,dim([1 2 4]));
                
%                 save(['slice',num2str(i,'%02d')],'s','-v6');
                clear s;
            end
            save aci_data d varMap
            selectRoi('direct');
        end
        function closeGetData(varargin)
            delete(fig_var);
        end
    end
%% ROI Selection
    function selectRoi(varargin)
        if ~strcmp(varargin{1}, 'direct')
            pathName = uigetdir(cd, 'Select directory containing ACI data!');
            if ~pathName
                return;
            end
            cd(pathName);
            try
                varMap = struct2array(load('aci_data','varMap'));
                d    = struct2array(load('aci_data','d'));
                clear w;
            catch %#ok
                warndlg('Directory does not contain ACI data! Select different directory or create ACI data first!');
                return
            end
            try
                roiList = struct2array(load('aci_roiList','roiList'));
                labels = struct2array(load('aci_roiList','labels'));
                clear w;
            catch  %#ok
                roiList = [];
                labels = [];
            end
        else
            roiList = [];
            labels = [];
        end
        set(txt_status, 'String', 'Select ROIs in ACI-Roi Selector. Close ROI Selector to continue.');
        drawnow;
        if isempty(roiList)
            [roiList labels] = aci_roiSelector(varMap,size(d,4),'Data',d);
        else
            [roiList labels] = aci_roiSelector(varMap,size(d,4),'Data',d,'roiList',roiList,'labels',labels);
        end
        if isempty(roiList)
            return
        end
        save aci_roiList roiList labels
        set(txt_status, 'String', 'Wait while correlation maps are calculated!');
        drawnow;
        % Data size
        dim = size(d);
        % Normalize time traces
        nRois = size(roiList,2);
        traces = zeros(dim(4),nRois);
        for ii=1:nRois
            t = roiList(ii).timeTrace';
            t = t - mean(t(:));
            traces(:,ii) =  t/norm(t);
        end

        % Calculate correlation maps
        cMaps = zeros(dim(1)*dim(2), dim(3), nRois, 'single');
        for ii=1:dim(3)
%             load(['slice',num2str(i,'%02d')]);
            s = reshape(squeeze(d(:,:,ii,:)),[dim(1)*dim(2) dim(4)]);
            cMaps(:,ii,:) = s * traces;
        end
        cMaps = reshape(cMaps, [dim(1:3) nRois]);
        save aci_corrMaps cMaps
        set(txt_status, 'String', 'Done.');
        drawnow;
        matVis(cMaps)
        clear d
        clear cMaps
    end

%% Close GUI


    function closeAll(varargin)
        try delete(fig_var); end %#ok
        delete(mainGui);
        cd(oldPath);
    end


end