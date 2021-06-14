function [roiList labels] = roiSelector(varMap, tDim, varargin)  %varMap, tDim, data, fileNameBase, roiList, labels)

%Arguments:
%   varMap       - Variance map
%   tDim          - Number of time points used for analysis

% Exactly one of the following two argument pairs has to be specified
%   'Data', data            - Data to be used for analysis
%   'Files', fileNameBase   - Initial (common) string of filenames of files containing
%                                   data in case data are not included directly

%Optional ('ParName', parValue)
%   'roiList', roiList      - Previously generated list of rois
%   'labels', labels        - Label matrix that was generated together with roiList

maxVal = max(varMap(:));
minVal = min(varMap(:));
dim = [size(varMap), tDim];
selectFromCorr = 2;

%% Read optional Arguments

%Check for optional arguments
numberOptionalArgs = (nargin - 2) / 2;
if numberOptionalArgs 
    optionalArgs = 1:2:numberOptionalArgs*2;
end

withRois = 0;
%Read optional arguments
for i = 1:numberOptionalArgs
    identifier = varargin{optionalArgs(i)};
    val = varargin{optionalArgs(i)+1};
    switch  identifier
        case 'Data'
            withData = 1;
            data = val;
        case 'Files'
            fileNameBase = val;
            withData = 0;
        case 'roiList'
            roiList = val;
            withRois = 1;
        case 'labels'
            labels = val;
    end
end
%% GUI
scrsz = get(0,'ScreenSize');
gui = figure('Units', 'Pixel', 'Position', [100 300 160 scrsz(4)-500], 'Name', 'ACI - RS',...
    'MenuBar', 'none', 'Resize', 'off', 'NumberTitle', 'off',...
    'CloseRequestFcn', @closeGui);

%Status Text
uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [0 15 180 20], 'String','______________________________________________________', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'Bold');
txt_status = uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [3 0 160 20], 'String', '', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'Bold');

%Roi Controls
%Listbox for Rois
uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [45 scrsz(4)-530 180 20], 'String','List of Rois', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'Bold');
roiListbox = uicontrol(gui, 'Style', 'listbox', 'Units', 'Pixel', ...
    'BackgroundColor', get(gui, 'Color'),'Value',0,...
    'Position', [10 190 140 scrsz(4)-720], 'String', '','Max',3, 'Callback', @listboxCallback);
%Delete Button
btDeleteRoi = uicontrol(gui, 'Style', 'pushbutton', 'Units', 'Pixel',...
    'BackgroundColor', get(gui, 'Color'),'Position', [10 155 140 20],...
    'String', 'Delete selected ROI(s)', 'Callback', @deleteRoi);
% btSaveRoi = uicontrol(gui, 'Style', 'pushbutton', 'Units', 'Pixel',...
%     'BackgroundColor', get(gui, 'Color'),'Position', [10 145 65 20],...
%     'String', 'Roi >>', 'Callback', @saveRoi);
% btExportTimeTraces = uicontrol(gui, 'Style', 'pushbutton', 'Units', 'Pixel',...
%     'BackgroundColor', get(gui, 'Color'),'Position', [85 145 65 20],...
%     'String', 'Traces >>', 'Callback', @exportTimeTraces);
%Popmenu for selection method
% pop_selMethod = uicontrol(gui, 'Style', 'popupmenu', 'Units', 'Pixel',...
%     'BackgroundColor', get(gui, 'Color'),'Position', [85 165 65 20],...
%     'String', {'Corr ';'Corr+Erode';'Volume'}, 'Callback', @setSelectionMethod, 'Value', selectFromCorr);
%z Controls
%Text
uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [55 89 60 17], 'String', 'z Position', ...
    'HorizontalAlignment', 'center', 'FontWeight', 'Bold');
%z slider 
zSlider = uicontrol(gui, 'Style', 'Slider', 'Callback', {@zPosCallback,'sld'}, 'Units', 'Pixel', ...
    'Position', [10 75 90 17], 'Min', 1, 'Max', dim(3), ...
    'SliderStep', [min(1,1/(dim(3)-1)) min(1,10/(dim(3)-1))],'Value',1);
if  usejava('awt') % java enabled -> use it to update while dragging
    hListeners = handle.listener(zSlider,'ActionEvent',@zPosCallback);
    set(zSlider, 'Callback', '');
    setappdata(gui,'sliderListeners', hListeners); % -> hListeners can be list of handles
end
%Text slider position
zPosTxt = uicontrol(gui, 'Style', 'Edit', 'Callback', {@zPosCallback,'txt'}, 'Units', 'Pixel', ...
    'Position', [105 75 22 17], 'String', '1');
%Text z max
uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [127 73 20 17], 'String', [' / ', num2str(dim(3))], ...
    'HorizontalAlignment', 'right');

%Colormap and Adjustment
uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [35 53 100 17], 'String', 'Display Controls', ...
    'HorizontalAlignment', 'center', 'FontWeight', 'Bold');
%Colormap
popLut = uicontrol(gui, 'Style', 'popupmenu', 'Callback', @colormapCallback, 'Units', 'Pixel', ...
    'Position', [10 35 80 20], 'String', {'Gray'; 'Jet'; 'HSV'; 'Hot'; 'Cool'},...
    'Value',2, 'TooltipString', 'Choose colormap');  %; 'Rainbow1';'Rainbow2';'Rainbow3';'Rainbow4';'Gray (Range)'}
%Contrast slider min 
contrastMinSlider = uicontrol(gui, 'Style', 'Slider', 'Callback', @updateImageWindow, 'Units', 'Pixel', ...
    'Position', [95 45 55 10], 'Min', minVal, 'Max', maxVal, ...
    'SliderStep', [1/(maxVal-minVal) 10/(maxVal-minVal)],'Value',minVal);
%Contrast slider max
contrastMaxSlider = uicontrol(gui, 'Style', 'Slider', 'Callback', @updateImageWindow, 'Units', 'Pixel', ...
    'Position', [95 34 55 10], 'Min', minVal, 'Max', maxVal, ...
    'SliderStep', [1/(maxVal-minVal) 10/(maxVal-minVal)],'Value',maxVal);

%Interval values for correlation
uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [8 125 200 17], 'String', 'Volume for Corr Analysis', ...
    'HorizontalAlignment', 'left', 'FontWeight', 'Bold');
uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [10 108 10 17], 'String', 'x:');
etxt_xSize = uicontrol(gui, 'Style', 'Edit',  'Units', 'Pixel', ...
    'Position', [25 110 25 17], 'String', '3');
uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [60 108 10 17], 'String', 'y:');
etxt_ySize = uicontrol(gui, 'Style', 'Edit',  'Units', 'Pixel', ...
    'Position', [75 110 25 17], 'String', '3');
uicontrol(gui, 'Style', 'Text', 'Units', 'Pixel', 'BackgroundColor', get(gui, 'Color'),...
    'Position', [110 108 10 17], 'String', 'z:');
etxt_zSize = uicontrol(gui, 'Style', 'Edit',  'Units', 'Pixel', ...
    'Position', [125 110 25 17], 'String', '1');

%% varMap Window
varMapWin = figure('Units', 'Pixel', 'Position', [300 300 700 scrsz(4)-500],'Name', 'ACI - ROI Selector: Variance Map', ...
        'MenuBar', 'none', 'NumberTitle', 'off', 'WindowButtonMotionFcn', @mouseMotion,...
        'WindowButtonDownFcn', @buttonDownCallback, 'CloseRequestFcn', ''); 

%% Plot Window
plotWin = figure('Units', 'Pixel', 'Position', [100 10 900 240],'Name', 'ACI - ROI Selector: Plot Window', ...
        'MenuBar', 'none', 'NumberTitle', 'off', 'WindowButtonMotionFcn', @mouseMotion,...
        'WindowButtonDownFcn', @buttonDownCallback, 'CloseRequestFcn', ''); 

%% Adjust for 2D data
if ndims(varMap)== 2
    set([zSlider etxt_zSize zPosTxt], 'Enable', 'off');
    set(etxt_zSize, 'String', '0');
    dim([3 4]) = [1 dim(3)];
    if withData && ndims(data) == 3
        data = permute(data,[1 2 4 3]);
    end
end

%% Display first image and set initial parameters
set(0, 'CurrentFigure', varMapWin);
varMapIm = imagesc(varMap(:,:,1));
varMapAx = gca;
axis equal off ij;
colormap(gray(255));
set(gca, 'Position', [0 0 1 1]);
%Initial values
zPos = 1;
if withRois
    numberRoi = size(roiList,2);
    if withData
        data = reshape(data,[],dim(4));
    end
    for i=1:numberRoi
        s{i} = [num2str(i,'%03d'),': ',roiList(i).name];
        roiText(i) = text(max([roiList(i).index.x])+1,nanmean([roiList(i).index.y]),num2str(roiList(i).number),'Color','w','HorizontalAlignment','left','Parent',varMapAx);
        if withData
            roiList(i).timeTrace = nanmean(data(labels(:)==i,:),1);
        end
    end
    set(roiListbox,'String',s,'Value',1);
    listboxCallback;
else
    numberRoi = 0;
    roiText = [];
    roiList = [];
    labels = zeros(size(varMap),'uint16');
end
if withData 
    data = reshape(data,dim);
end
cmap = [];
colormapCallback;
updateImageWindow;
uiwait(gui);
%% Callbacks
%z Position
    function zPosCallback(hObject,event,caller)
        if nargin == 2
            caller = 'sld';
        end
        switch caller
            case 'sld'
                zPos = round(get(zSlider, 'Value'));
                set(zPosTxt, 'String', num2str(zPos));
            case 'txt'
                z = str2double(get(zPosTxt,'String'));
                if z < 1 || z > dim(3)
                    set(zPosTxt, 'String', num2str(zPos));
                    return
                end
                zPos = z;
                set(zSlider,'Value',zPos);
        end
       %set(varMapIm, 'CData', varMap(:,:,zPos));       
       updateImageWindow;
    end
%Contrast
    function contrastSliderCallback(varargin)
       cMin = get(contrastMinSlider, 'Value');
       cMax = get(contrastMaxSlider, 'Value');
       if cMax <= cMin
           cMin = cMax * 0.99;
       end
       set(varMapAx, 'CLim', [cMin cMax]);
    end
%Colormap
    function colormapCallback(varargin)
       switch get(popLut, 'Value')
            case 1
                cmap = gray(255);
                cmap(255,:) = [1 0 0];
                cmap(256,:) = [0.5 1 0.5];
            case 2
                cmap = jet(255);
                cmap(255,:) = [0 0 0];
                cmap(256,:) = [0 .6 0];
            case 3
                cmap = hsv(255);
                cmap(255,:) = [0 0 0];
                cmap(256,:) = [1 1 1];
            case 4
                cmap = hot(255);
                cmap(255,:) = [0 0 1];
                cmap(256,:) = [0 .7 0];
            case 5
                cmap = cool(255);
                cmap(255,:) = [0 0 1];
                cmap(256,:) = [0 1 0];
            case {6;7;8;9;10}
                cmap = color_map(255, get(popLut, 'Value')-5);
                cmap(255,:) = [0 0 0];
                cmap(256,:) = [0.9 0.9 0];
       end
       colormap(varMapAx, cmap);
    end
%Roi Selection (Listbox)
    function listboxCallback(varargin)
        roiSel = get(roiListbox, 'Value');
        set(0, 'CurrentFigure', plotWin);
        for i=1:numel(roiSel)
            plot(roiList(roiSel(i)).timeTrace);
            leg{i} = roiList(roiSel(i)).name;
            hold all;
        end
        hold off;
        legend(leg);
        updateImageWindow; 
    end

%Delete Roi
    function deleteRoi(varargin)
        roiSel = get(roiListbox, 'Value');
        roiList(roiSel) = [];
        for i=1:numel(roiSel)
            labels(labels == roiSel(i)) = 0;
        end
        s = get(roiListbox, 'String');
        s(roiSel) = [];
        for i=1:size(s,1)
            s{i}(1:3) = num2str(i,'%03d');
        end
        set(roiListbox,'Value',1);
        set(roiListbox, 'String',s);
        delete(roiText(roiSel));
        roiText(roiSel) = [];
        numberRoi = numberRoi - size(roiSel,2);
        if numberRoi > 0
            set(roiListbox, 'Value',1);
        end
        if numberRoi > 0
            listboxCallback;
        else
            set(roiListbox, 'String', '');
            updateImageWindow;
        end
    end
%Save Rois
    function saveRoi(varargin)
       save('roiList', 'roiList');
       save('labels','labels');
       assignin('base','roiList',roiList);
       assignin('base','labels',labels);
    end
%Export Time Traces
    function exportTimeTraces(varargin)
        for i=1:numberRoi
            tt(i,:) = roiList(i).timeTrace;
        end
        assignin('base','timeTraces',tt);
    end

%% Mouse controls
    function mouseMotion(varargin)
       if gcf == varMapWin
             %Get current position of cursor
             p = round(get(gca, 'CurrentPoint'));
             p = p(1,1:2);
             p(1) = max([min([p(1),dim(2)]),1]);
             p(2) = max([min([p(2),dim(1)]),1]);
             %Display position / value of current point while inside Image
             %window
             set(varMapWin, 'Name', ['varMap - Pos: (', num2str(p(2:-1:1)),'  ',num2str(zPos),')  Val: ',...
                 num2str(varMap(p(2),p(1),zPos))]);
       end
    end

    function buttonDownCallback(varargin)
       if strcmp(get(gcf,'SelectionType'),'normal')  
            %Zoom by drawing zoom region
            p1 = round(get(gca,'CurrentPoint'));
            rbbox;   
            p2 = round(get(gca,'CurrentPoint'));
            if p1 == p2
                return;
            end
            coord(1,:)    = round(p1(1,1:2));
            coord(2,:)    = round(p2(1,1:2));
            if coord(1,1) == coord(2,1) || coord(1,2) == coord(2,2)
                return;
            end
            if coord(1,1) > coord(2,1) 
                coord(:,1) = flipdim(coord(:,1),1);
            end
            if coord(1,2) > coord(2,2) 
                coord(:,2) = flipdim(coord(:,2),1);
            end
            axis(varMapAx, [coord(1,1), coord(2,1),coord(1,2),coord(2,2)]); 
        elseif strcmp(get(gcf,'SelectionType'),'alt')
            %Unzoom
            axis(varMapAx, [1 dim(2) 1 dim(1)]); 
       elseif strcmp(get(gcf,'SelectionType'),'extend')
            %Detect Cell
            p1 = round(get(gca,'CurrentPoint'));
            selectCell(p1);
       end
    end

%% Select cell
    function selectCell(varargin)
        %Get selected point
        p1 = round(get(gca, 'CurrentPoint'));
        p1 = p1(1,1:2);
        %Check if click is inside cell
%         if labels(p1(2),p1(1),zPos) ~= 0
%            set(roiListbox, 'Value',  labels(p1(2),p1(1),zPos));
%            listboxCallback;
%            return;
%         end
        %Set Status Bar
        set(txt_status, 'String', 'Defining ROI ...'); drawnow;
        %Get size of volume
        xSize = str2double(get(etxt_xSize, 'String'));
        ySize = str2double(get(etxt_ySize, 'String'));
        zSize = str2double(get(etxt_zSize, 'String'));
        % Draw rectangle indicating volume
        r = rectangle('Position',[p1(1)-ySize p1(2)-xSize 2*ySize+1 2*xSize+1],'Parent',varMapAx,...
            'EdgeColor','w');
        %Extract volume for corr analysis from complete data set
        %tic;
        if withData
            corrVol = data(max(1,p1(2)-xSize):min(dim(1),p1(2)+xSize),...
                max(1,p1(1)-ySize):min(dim(2),p1(1)+ySize),...
                max(1, zPos-zSize):min(dim(3), zPos+zSize),...
                1:dim(4));
        else
            corrVol = readCustomTifParts(fileNameBase,...
                [max(1,p1(2)-xSize) min(dim(1),p1(2)+xSize)],...
                [max(1,p1(1)-ySize) min(dim(2),p1(1)+ySize)],...
                [max(1, zPos-zSize) min(dim(3), zPos+zSize)],...
                [1 dim(4)],[1 1]);
        end
        %Position of corner of volume in data set
        volPos = [max(1,p1(2)-xSize) max(1,p1(1)-ySize) max(1,zPos-zSize)];
        %Size of extracted volume
        volDim = size(corrVol);
        cellInVol = ones(size(corrVol(:,:,:,1)));
        timeTrace = squeeze(nanmean(nanmean(nanmean(corrVol,1),2),3));
        set(txt_status,'String','Done');
        %Embed in volume of data size
        vol = zeros(dim(1:3),'uint8');
        vol(volPos(1):volPos(1)+volDim(1)-1,...
            volPos(2):volPos(2)+volDim(2)-1,...
            volPos(3):volPos(3)+volDim(3)-1) = cellInVol;
        %Set roi counter
        numberRoi = numberRoi + 1;
        if numberRoi == 1
            roiList(numberRoi).number = 1;
        else
            roiList(numberRoi).number = max([roiList.number])+1;
        end
        %Get surround of Roi using contour function
        hold on
        roiList(numberRoi).mask = vol;
        labels(vol==1) = numberRoi;
        for i=1:dim(3)
            if any(any(vol(:,:,i)))
                [roi,hroi] = contour(vol(:,:,i));
                delete(hroi);
                roiList(numberRoi).index(i).x = roi(1,:);
                roiList(numberRoi).index(i).y = roi(2,:);
                roiList(numberRoi).index(i).x(roiList(numberRoi).index(i).x>dim(2)) = dim(2);
                roiList(numberRoi).index(i).y(roiList(numberRoi).index(i).y>dim(1)) = dim(1);
            else
                roiList(numberRoi).index(i).x = [];
                roiList(numberRoi).index(i).y = [];
            end
        end
        hold off
        %Draw final polygons and numbers in image 
        %windows and delete temporary lines
        roiText(numberRoi) = text(max([roiList(numberRoi).index.x])+1,nanmean([roiList(numberRoi).index.y]),num2str(roiList(numberRoi).number),'Color','w','HorizontalAlignment','left');
        roiList(numberRoi).name = ['Roi',num2str(roiList(numberRoi).number)];
        roiList(numberRoi).timeTrace = timeTrace;
        %Add to listbox 
        s = get(roiListbox, 'String');
        s{size(s,1)+1} = [num2str(size(get(roiListbox,'String'),1)+1,'%03d'),': ',roiList(numberRoi).name];
        if numberRoi == 1
            set(roiListbox,'Value', numberRoi,'String', s);
        else
            set(roiListbox,'String', s, 'Value', numberRoi);
        end
        delete(r);
        updateImageWindow;
        listboxCallback;
    end

%% Update Image Window
    function updateImageWindow(varargin)
       cMin = get(contrastMinSlider, 'Value');
       cMax = get(contrastMaxSlider, 'Value');
       if cMax <= cMin
           cMin = cMax * 0.99;
       end
       im = varMap(:,:,zPos);
       im = (im-cMin)/(cMax-cMin) * 255;
       im(im<0) = 0;
       im(im>254) = 254;
       if numberRoi > 0
           for i=1:numberRoi
               ind=sub2ind(dim(1:2),...
                           ceil(roiList(i).index(zPos).y),...
                           ceil(roiList(i).index(zPos).x));
               if  numel(roiList(i).index(zPos).y) == 0
                   set(roiText(i),'Visible', 'off');
               else
                   set(roiText(i),'Visible', 'on');
               end
               im(ind) = 255;
           end
           selRoi = get(roiListbox,'Value');
           for i=1:numel(selRoi)
               ind=sub2ind(dim(1:2),...
                           ceil(roiList(selRoi(i)).index(zPos).y),...
                           ceil(roiList(selRoi(i)).index(zPos).x));
               im(ind) = 256;
           end
       end
       set(varMapIm, 'CData', im);
       set(varMapAx, 'CLim', [0 257]);
    end

%% Close All Windows
    function closeGui(varargin)
        delete(gui);
        delete(varMapWin);
        delete(plotWin)
    end
end