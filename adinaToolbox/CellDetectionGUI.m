 function varargout = CellDetectionGUI(varargin)
% CELLDETECTIONGUI MATLAB code for CellDetectionGUI.fig
%      CELLDETECTIONGUI, by itself, creates a new CELLDETECTIONGUI or raises the existing
%      singleton*.
%
%      H = CELLDETECTIONGUI returns the handle to a new CELLDETECTIONGUI or the handle to
%      the existing singleton*.
%
%      CELLDETECTIONGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CELLDETECTIONGUI.M with the given input arguments.
%
%      CELLDETECTIONGUI('Property','Value',...) creates a new CELLDETECTIONGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before CellDetectionGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to CellDetectionGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help CellDetectionGUI

% Last Modified by GUIDE v2.5 06-Oct-2014 15:28:13

% Begin initialization code - DO NOT EDI
%   Author:  Ferran Diego, HCI, IWR, University of Heidelberg
%   Contact: ferran.diego@iwr.uni-heidelberg.de
%   $Date: 2014-10-01 $
%   $Revision: 1 $
%   
%   [1] F. Diego, S. Reichinnek, M. Both, F. A. Hamprecht. "Automated Identification 
%    of Neuronal Activity from Calcium Imaging by Sparse Dictionary Learning".
%    ISBI 2013. Proceedings, (2013), 1058-1061".

startup
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @CellDetectionGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @CellDetectionGUI_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT



% --- Executes just before CellDetectionGUI is made visible.
function CellDetectionGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to CellDetectionGUI (see VARARGIN)

global MFmethods;
global cMethod;
global paramgui;

initSegmentationMethods()

cMethod.MethodSelected = 1;
cMethod.subMethodSelected = 1;
cMethod.paramSelected = 1;

MFmethods.List = get_methods(pwd);
for i = 1:length(MFmethods.List)
    MFmethods.Method{i} = MF_methods(MFmethods.List{i});
end

initPreprocessingStep()


% Choose default command line output for CellDetectionGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
set(handles.figure1,'CloseRequestFcn',@closeGUI);

axes(handles.DisplaySequences);
imagesc(ones(512,512,3))
axis off
axis equal

%Initialize methods, submethods, parameters



%MFmethods.Method{2}.List = {'NMF','Sparse NMF'};

% UIWAIT makes CellDetectionGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = CellDetectionGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure


%% All the variables used must be cleaned



% ----------------------------------------------------------------------------------------------------------------------
% Beginning Project Actions (Save/Load)
% ----------------------------------------------------------------------------------------------------------------------
function closeGUI(src,evnt)
selection = questdlg('Do you want to close the GUI?','Close Request Function','Yes','No','Yes');
switch selection
    case 'Yes'
        ListVariables = get_Variablelist();
        for i = 1:length(ListVariables),
                %eval(['global ' ListVariables{i} ';'])
                eval(['clear global ' ListVariables{i} ';'])
        end
        delete(gcf)
    case 'No'
        return
end




% --- Executes during object deletion, before destroying properties.
function ProjectActionsTag_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to LoadButton (see GCBO)
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function filemenu_Callback(hObject, eventdata, handles)
% hObject    handle to filemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function saveproject_Callback(hObject, eventdata, handles)
% hObject    handle to saveproject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
ListVariables = get_Variablelist();

[FileName,PathName] = uiputfile({'*.mat'},'File Selector');
if isequal(FileName,0)
    %disp('User pressed cancel')
else
    disp(['User selected File: ', [PathName FileName]])
    str = ['save('''  PathName FileName ''''];
    for i = 1:length(ListVariables),
        eval(['global ' ListVariables{i} ';'])
        str = [str ',''' ListVariables{i} ''''];
    end
    str = [str ',''-v7.3'');'];
    disp(['Saving file'])
    disp(str)
    eval(str);
    disp('Saved file')
end

% --------------------------------------------------------------------
function loadprojectmenu_Callback(hObject, eventdata, handles)
% hObject    handle to loadprojectmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile({'*.mat'},'File Selector');
if isequal(FileName,0)
    %disp('User pressed cancel')
else
    disp(['User selected File: ', [PathName FileName]])
    S = load([PathName FileName]);
    names = fieldnames(S);
    flag = 0;
    %Extract data
    if ~isempty(names)
        ListVariables =  get_Variablelist();
        for i = 1:length(ListVariables),
            if isfield(S,ListVariables{i})
                 eval([' global ' ListVariables{i} ';'])
                 eval([ 'tmp = S.' ListVariables{i} ';']);
                 if ~isempty(tmp)
                     eval([ListVariables{i} ' = tmp;'])
                 else
                    disp([ ListVariables{i} ' is not updated. There is not inside the matlab file.'])                     
                 end
            else
                disp([ ListVariables{i} ' is not updated. There is not inside the matlab file.'])
            end
        end
%         nFields = length(names);
%         for i = 1:nFields,
%             if fieldname_GUI(names{i})
%                 eval(['global ' names{i} ';'])
%                 eval([names{i} ' = S.' names{i} ';']);
%             end
%         end  
        updateCellExtractionActions(handles);
        updateMatrixFactorizationActions(handles);
        updatePreprocessing(handles);
        updateEnhancement(handles);
        updateSequence(handles)
        updateDisplay(handles)
        updateArtificialSequence(handles)
        axes(handles.DisplaySequences);
        imagesc(max(data_orig,[],3))
        colorbar; axis off,axis equal
        disp(['Load mat file'])
    else
        %disp('Input sequence is empty. Please select a new sequence (nrows x ncols x nFrames)');
        msgbox('Input file is empty. Please select a new file','Error','error','modal');
    end
end



% -----------------------------------------------------------------------------------------------------------------------
% Ending Project Actions (Save/Load)
% -----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Beginning Sequence Actions
% ----------------------------------------------------------------------------------------------------------------------

% --- Executes on button press in loadDir.
function loadDir_Callback(hObject, eventdata, handles)
% hObject    handle to loadDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global paramgui

if ~isfield(paramgui,'dirBase')
    tmpDir = uigetdir();
    if ~isequal(tmpDir,0)	
        paramgui.dirBase = [tmpDir '/'];
        set(handles.SequencePathTag,'String', paramgui.dirBase)
    end
else
    tmpDir = uigetdir(paramgui.dirBase);
    if ~isequal(tmpDir,0)	
        paramgui.dirBase = [tmpDir '/'];
        set(handles.SequencePathTag,'String', paramgui.dirBase)
    end
end



function SequencePathTag_Callback(hObject, eventdata, handles)

global paramgui

PathName = get(hObject,'String');
if isequal(PathName,0)
    %disp('User pressed cancel')
else
    
    if exist(PathName,'dir')
        paramgui.dirBase = PathName;
    else
        msgbox('Incorrect folter. Please select a new folder','Error','error','modal');
    end
end

function selectWildCard_Callback(hObject, eventdata, handles)
% hObject    handle to selectWildCard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of selectWildCard as text
%        str2double(get(hObject,'String')) returns contents of selectWildCard as a double

global paramgui

PathName = get(hObject,'String');
if isequal(PathName,0)
    %disp('User pressed cancel')
else
    files = dir([paramgui.dirBase PathName]);
    if ~isempty(files),
        paramgui.fileExt = PathName;
        disp([ num2str(length(files)) ' files found'])
    else
        msgbox('Incorrect Wild Card. 0 files found','Error','error','modal');
    end
end


% --- Executes during object creation, after setting all properties.
function selectWildCard_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectWildCard (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global paramgui
paramgui.step = 1;


function selectNumChannels_Callback(hObject, eventdata, handles)
% hObject    handle to selectNumChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of selectNumChannels as text
%        str2double(get(hObject,'String')) returns contents of selectNumChannels as a double
global paramgui

step = get(hObject,'String');
if isequal(step,0)
    %disp('User pressed cancel')
else
    step = str2double(step);
    paramgui.step = max(step,1);
    disp(['number of channels: ' num2str(paramgui.step)])
    set(handles.selectNumChannels,'String',num2str(paramgui.step))
end


% --- Executes during object creation, after setting all properties.
function selectNumChannels_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selectNumChannels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in RunSequence.
function RunSequence_Callback(hObject, eventdata, handles)

global paramgui
global data
global displayParam
global data_orig


if isfield(paramgui,'dirBase') && isfield(paramgui,'fileExt')
    files = dir([paramgui.dirBase paramgui.fileExt]);
    disp('Loading image sequence')
    data = extractingDataGreenChannel(paramgui.dirBase,files,paramgui.step);
    data_orig = data;
    [nrows,ncols,frames] = size(data);
    disp(['Image sequence is loaded with ' num2str(nrows) ' rows, ' num2str(ncols) ' cols and ' num2str(frames) ' frames.'])
    paramgui.data = data;
    paramgui.nrows = nrows;
    paramgui.ncols = ncols;
    paramgui.nFrames = frames;
    axes(handles.DisplaySequences);
    imagesc(mean(data,3))
    colorbar; axis off,axis equal
    displayParam.endF = frames;
    updateDisplay(handles)
    paramgui.displayParam = displayParam;
else
    if isfield(paramgui,'dirBase') 
        msgbox('Please select a folder first.','Error','error','modal');
    else
        msgbox('Please select a wild card.','Error','error','modal');
    end
end

function SequencePathTag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% ----------------------------------------------------------------------------------------------------------------------
% Ending Sequence Actions
% ----------------------------------------------------------------------------------------------------------------------


% ----------------------------------------------------------------------------------------------------------------------
% Beginning Preprocessing Actions
% ----------------------------------------------------------------------------------------------------------------------


% --- Executes on selection change in SubtractionListTag.
function SubtractionListTag_Callback(hObject, eventdata, handles)

global preproc_params
num = get(hObject, 'Value');
% disp(['Selected ' num2str(num)])
% disp(preproc_params.ListSubs{num})
preproc_params.selSubs = num;
preproc_params.modeSubs = preproc_params.ListSubs{num};


function RunButtonPreProcessing_Callback(hObject, eventdata, handles)
global X 
global data
global preproc_params

if ~isempty(data)
     disp('Preprocessing data ... '),
     tic
     X = preprocessingData(data,preproc_params.modeSubs,preproc_params.modeNorm);
     time_pre = toc;
     disp(['Preprocessing is done. ' num2str(time_pre) ' sec.'])
else
    msgbox('Please load a sequence before setting the number of basis functions.','Error','error','modal');
end

function updatePreprocessing(handles)
global preproc_params
set(handles.NormalizationListTag,'Value',preproc_params.selNorm);
set(handles.SubtractionListTag,'Value',preproc_params.selSubs);


function NormalizationListTag_Callback(hObject, eventdata, handles)
global preproc_params
num = get(hObject, 'Value');
preproc_params.selNorm = num;
preproc_params.modeNorm = preproc_params.ListNorm{num};




% --- Executes during object creation, after setting all properties.
function NormalizationListTag_CreateFcn(hObject, eventdata, handles)

global preproc_params
preproc_params.ListNorm = {'non','l1-norm','l2-norm','squared l2-norm'};
preproc_params.selNorm = 1;
preproc_params.modeNorm = preproc_params.ListNorm{preproc_params.selNorm};
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject,'String',preproc_params.ListNorm)

function SubtractionListTag_CreateFcn(hObject, eventdata, handles)
global preproc_params
preproc_params.ListSubs = {'non','mean','median','min'};
preproc_params.selSubs = 1;
preproc_params.modeSubs = preproc_params.ListSubs{preproc_params.selSubs};

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String',preproc_params.ListSubs)

% ----------------------------------------------------------------------------------------------------------------------
% Ending Preprocessing Actions
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% Beginning Matrix Factorization Actions
% ----------------------------------------------------------------------------------------------------------------------

% --- Executes on button press in SaveMFButton.
function SaveMFButton_Callback(hObject, eventdata, handles)

ListVariables = get_Variablelist('mf');
[FileName,PathName] = uiputfile({'*.mat'},'File Selector');
if isequal(FileName,0)
    %disp('User pressed cancel')
else
    disp(['User selected File: ', [PathName FileName]])
    str = ['save('''  PathName FileName ''''];
    for i = 1:length(ListVariables),
        eval(['global ' ListVariables{i} ';'])
        str = [str ',''' ListVariables{i} ''''];
    end
    str = [str ');'];
    disp(['Saving Extracted Cells file'])
    disp(str)
    eval(str);
    disp('Saved file')
end


% --- Executes on button press in LoadMFButton.
function LoadMFButton_Callback(hObject, eventdata, handles)
global Dicts


ListVariables = get_Variablelist('mf');
[FileName,PathName] = uigetfile({'*.mat'},'File Selector');
if isequal(FileName,0)
    %disp('User pressed cancel')
else
    disp(['User selected File: ', [PathName FileName]])
    S = load([PathName FileName]);
    names = fieldnames(S);
    %Extract data
    if ~isempty(names)        
        if ~isempty(Dicts)
            selection = questdlg('Do you want to replace the current matrix factorization?','Close Request Function','Yes','No','Yes');
            switch selection
                case 'Yes'            
                    for i = 1:length(ListVariables),
                        if isfield(S,ListVariables{i})
                             eval([' global ' ListVariables{i} ';'])
                             eval([ ListVariables{i} ' = S.' ListVariables{i} ';']);
                        else
                            disp([ ListVariables{i} ' is not updated. There is not inside the matlab file.'])
                        end
                    end
                case 'No'
                    disp(['File has not been loaded: ' PathName FileName])
            end
        else            
            for i = 1:length(ListVariables),
                if isfield(S,ListVariables{i})
                     eval([' global ' ListVariables{i} ';'])
                     eval([ ListVariables{i} ' = S.' ListVariables{i} ';']);
                else
                    disp([ ListVariables{i} ' is not updated. There is not inside the matlab file.'])
                end
            end
        end       
        updateCellExtractionActions(handles);
        updateMatrixFactorizationActions(handles);
        updatePreprocessing(handles);
        updateDisplay(handles)
        disp(['Load mat file.' PathName FileName])        
    else
        %disp('Input sequence is empty. Please select a new sequence (nrows x ncols x nFrames)');
        msgbox('Input file is empty. Please select a new file','Error','error','modal');
    end
end





% --- Executes on button press in RunMatrixFactorizationTag.
function RunMatrixFactorizationTag_Callback(hObject, eventdata, handles)
global cMethod;
global MFmethods;
global X
global data
global Coeff
global Dict
global nGBasis
global ReconsDict

params = MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected};

if ~isempty(data)
    if nGBasis >0
        if isempty(X)            
           X = data; 
        end
        disp('Inferring data ... '),
        [Coeff,Dict,ReconsDict] = feval(params.type,X,nGBasis,params.all);
        disp(['Extracted Grouped Cells is done.'])
    else
         msgbox('Please set the number of basis functions bigger thanb zero.','Error','error','modal');
    end
else
    msgbox('Please load a sequence before setting the number of basis functions.','Error','error','modal');
end


function DisplayMatrixFactorizationTag_Callback(hObject, eventdata, handles)


% --- Executes on selection change in MethodListMF.
function MethodListMF_Callback(hObject, eventdata, handles)
global cMethod;
global MFmethods;
num = get(hObject, 'Value');
cMethod.MethodSelected = num;
updateMatrixFactorizationActions(handles)

function SubMethodListMF_Callback(hObject, eventdata, handles)
global cMethod;
global MFmethods;
num = get(hObject, 'Value');
cMethod.subMethodSelected = num;
updateMatrixFactorizationActions(handles);


function updateMatrixFactorizationActions(handles)
global cMethod;
global MFmethods;
global nGBasis
set(handles.MethodListMF,'String',MFmethods.List,'Value',cMethod.MethodSelected);
set(handles.SubMethodListMF,'String',MFmethods.Method{cMethod.MethodSelected}.List,'Value',cMethod.subMethodSelected );
set(handles.ParamaterListTag,'String',MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected}.subList,'Value',cMethod.paramSelected);
set(handles.nBasisVal,'String',num2str(nGBasis));


param = MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected};
val = getfield(param.all,param.subList{cMethod.paramSelected});
if val <0, tmp = '-';else  tmp = []; end
set(handles.ParameterListValue,'String',[tmp num2str(abs(val))]);

function nBasisVal_Callback(hObject, eventdata, handles)
global nGBasis
global data
num = str2num(get(hObject,'String'));
if ~isempty(data)
    maxBasis = size(data,3);
    if length(num) == 1 & num <=maxBasis & num >=0
       nGBasis = num;
       disp(['The number of basis functions is set to :' num2str(nGBasis) ])
    else
        msgbox(['The value should be a number in the range [1,' num2str(maxBasis) ']'],'Error','error','modal');
        set(hObject,'String','0')
    end  
else
    msgbox('Please load a sequence before setting the number of basis functions.','Error','error','modal');  
    set(hObject,'String','0')  
end


function ParamaterListTag_Callback(hObject, eventdata, handles)
global cMethod;
global MFmethods;

num = get(hObject, 'Value');
cMethod.paramSelected = num;
% param = MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected};
% val = getfield(param.all,param.subList{cMethod.paramSelected});
% if val <0, tmp = '-';else  tmp = []; end
% set(handles.ParameterListValue,'String',[tmp num2str(abs(val))]);
updateMatrixFactorizationActions(handles);



function ParameterListValue_Callback(hObject, eventdata, handles)
global cMethod;
global MFmethods;

value  = str2double(get(hObject,'String'));
paramStruct = MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected};
paramField = paramStruct.subList{cMethod.paramSelected};

newParam = feval(['set_parameters_' paramStruct.Type],paramStruct.all,{paramField,value});
MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected}.all = newParam;
param = MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected};
val = getfield(param.all,param.subList{cMethod.paramSelected});
if val <0, tmp = '-';else  tmp = []; end
set(hObject,'String',[tmp num2str(abs(val))]);



% --- Executes during object creation, after setting all properties.
function ParameterListValue_CreateFcn(hObject, eventdata, handles)
global cMethod;
global MFmethods;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if ~exist('MFmethods.Method','var'),
    initMethods()
end
param = MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected};
val = getfield(param.all,param.subList{cMethod.paramSelected});
if val <0, tmp = '-';else  tmp = []; end
set(hObject,'String',[tmp num2str(abs(val))]);

function ParamaterListTag_CreateFcn(hObject, eventdata, handles)
global cMethod;
global MFmethods;
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if ~exist('MFmethods.Method','var'),
    initMethods()
end
set(hObject,'String',MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected}.subList,'Value',cMethod.paramSelected);


function nBasisVal_CreateFcn(hObject, eventdata, handles)
global nGBasis
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
nGBasis = 0;
set(hObject,'String','0');


function SubMethodListMF_CreateFcn(hObject, eventdata, handles)
global cMethod;
global MFmethods;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if ~exist('MFmethods.Method','var'),
    initMethods()
end
set(hObject,'String',MFmethods.Method{cMethod.MethodSelected}.List,'Value',cMethod.subMethodSelected );

function MethodListMF_CreateFcn(hObject, eventdata, handles)
global cMethod;
global MFmethods;
%MFmethods.MethodList = {'SPAMS','NMF (Not implemented)'};

%mFactor.MethodMode = mFactor.MethodList{mFactor.MethodSelected};

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

if ~exist('MFmethods.Method','var'),
    initMethods()
end
set(hObject,'String',MFmethods.List,'Value',cMethod.MethodSelected);


% ----------------------------------------------------------------------------------------------------------------------
% Ending Matrix Factorization Actions
% ----------------------------------------------------------------------------------------------------------------------



% ----------------------------------------------------------------------------------------------------------------------
% Beginning Cell Extraction Actions
% ----------------------------------------------------------------------------------------------------------------------


function RunExtractedCells_Callback(hObject, eventdata, handles)
global Segmethods;
global cSegmentation;
global MFmethods
global cMethod
global nBasis;
global data
global Coeff
global ReconsDict
global X
global Coeff_single
global Dicts_single
global cNeurons
global ReconsDict_single


param_actual = Segmethods.param.all;
param_actual.Tsimilarity = Segmethods.ListSim{cSegmentation.SimilarSelected};
param_actual.type        = Segmethods.ListM{cSegmentation.MethodSelected};

if ~isempty(data)
    if ~isempty(ReconsDict)     
            disp('Extracting single cells ... '),
            if ~isempty(cNeurons) &&  ~param_actual.resetNeurons;
                 % ask if they want to update the current extracted cells                 
                 [cNeurons] = extractSingleCells(Coeff,ReconsDict,cNeurons,param_actual);
            else                
                 [cNeurons] = extractSingleCells(Coeff,ReconsDict,[],param_actual);
            end
            
            if ~isempty(cNeurons)                
                disp(['Extracted Single Cells is done.'])
                disp('Computing Temporal Evolution for the extracted single Cells')
                param.lambda = param_actual.lambda;
                param.pos    = getfield(MFmethods.Method{cMethod.MethodSelected}.params{cMethod.subMethodSelected}.all,'posD');
                [Coeff_single,Dicts_single,ReconsDict_single] = inferCoeffSingleCells(X,cNeurons,param);
                
                %% filter neurons that are not used for reconstructing
                [cNeurons,Coeff_single,Dicts_single,ReconsDict_single] = proneNeurons(cNeurons,Coeff_single,Dicts_single,ReconsDict_single,param_actual);
                
                disp('Temporal Evolution is computed')    
            else
                disp('No found single cells.')
            end
                
    else
         msgbox('Please run Matrix Factorization algorithm.','Error','error','modal');
    end
else
    msgbox('Please load a sequence before setting the number of basis functions.','Error','error','modal');
end

nBasis = length(cNeurons);
set(handles.DetectedCells,'String',num2str(nBasis));


function MethodExtractionList_Callback(hObject, eventdata, handles)
global cSegmentation;
num = get(hObject, 'Value');
cSegmentation.MethodSelected = num;
updateCellExtractionActions(handles);


function SimilarityFtList_Callback(hObject, eventdata, handles)
global cSegmentation;
num = get(hObject, 'Value');
cSegmentation.SimilarSelected = num;
updateCellExtractionActions(handles);

function DetectedCells_Callback(hObject, eventdata, handles)
global nBasis
set(hObject,'String',num2str(nBasis));



function SaveSingleExtractedCells_Callback(hObject, eventdata, handles)
ListVariables = get_Variablelist('segmentation');
[FileName,PathName] = uiputfile({'*.mat'},'File Selector');
if isequal(FileName,0)
    %disp('User pressed cancel')
else
    disp(['User selected File: ', [PathName FileName]])
    str = ['save('''  PathName FileName ''''];
    for i = 1:length(ListVariables),
        eval(['global ' ListVariables{i} ';'])
        str = [str ',''' ListVariables{i} ''''];
    end
    str = [str ');'];
    disp(['Saving Extracted Cells file'])
    disp(str)
    eval(str);
    disp('Saved file')
end




function LoadSingleExtractedCells_Callback(hObject, eventdata, handles)
global cNeurons
global Segmethods;
global cSegmentation;


param_actual = Segmethods.param.all;
param_actual.Tsimilarity = Segmethods.ListSim{cSegmentation.SimilarSelected};
param_actual.type        = Segmethods.ListM{cSegmentation.MethodSelected};

ListVariables = get_Variablelist('segmentation');
[FileName,PathName] = uigetfile({'*.mat'},'File Selector');
if isequal(FileName,0)
    %disp('User pressed cancel')
else
    disp(['User selected File: ', [PathName FileName]])
    S = load([PathName FileName]);
    names = fieldnames(S);
    %Extract data
    if ~isempty(names)        
        if ~isempty(cNeurons)
            selection = questdlg('Do you want to update the current extracted single components? No: Remove previous extraction. Yes: fuse load files and current extraction.','Close Request Function','Yes','No','Yes');
            switch selection
                case 'Yes'
                    if isfield(S,'cNeurons')
                        eval('tmpNeurons = S.cNeurons;')
                        cNeurons = identifyDistinctCells(cNeurons,tmpNeurons,param_actual);
                    end
                case 'No'
                    for i = 1:length(ListVariables),
                        if isfield(S,ListVariables{i})
                            eval([' global ' ListVariables{i}])
                            eval([ ListVariables{i} ' = S.' ListVariables{i} ';' ]);
                        else
                            disp([ ListVariables{i} ' is not updated. This variable does not exist in the matlab file.'])
                        end
                    end
            end
        else            
            for i = 1:length(ListVariables),
                if isfield(S,ListVariables{i})
                     eval([' global ' ListVariables{i} ';'])
                     eval([ ListVariables{i} ' = S.' ListVariables{i} ';']);
                else
                     disp([ ListVariables{i} ' is not updated. This variable does not exist in the matlab file.'])
                end
            end
        end       
        updateCellExtractionActions(handles);
        updateMatrixFactorizationActions(handles);
        updatePreprocessing(handles);
        updateDisplay(handles)
        disp(['Load mat file'])        
    else
        %disp('Input sequence is empty. Please select a new sequence (nrows x ncols x nFrames)');
        msgbox('Input file is empty. Please select a new file','Error','error','modal');
    end
end


function ParameterExtractionList_Callback(hObject, eventdata, handles)
global cSegmentation;
num = get(hObject, 'Value');
cSegmentation.paramSelected = num;
updateCellExtractionActions(handles)

function ParameterExtractionListVal_Callback(hObject, eventdata, handles)
global Segmethods;
global cSegmentation;


value  = str2double(get(hObject,'String'));
paramField = Segmethods.param.List{cSegmentation.paramSelected};
Segmethods.param.all = feval(str2func(['@set_parameters_Segmentation']),Segmethods.param.all,{paramField,value});
val = getfield(Segmethods.param.all,Segmethods.param.List{cSegmentation.paramSelected});
if val <0, tmp = '-';else  tmp = []; end
set(hObject,'String',[tmp num2str(abs(val))]);





function updateCellExtractionActions(handles)
global Segmethods;
global cSegmentation;
global nBasis

set(handles.MethodExtractionList,'String',Segmethods.ListM,'Value',cSegmentation.MethodSelected);
set(handles.SimilarityFtList,'String',Segmethods.ListSim,'Value',cSegmentation.SimilarSelected);
set(handles.ParameterExtractionList,'String',Segmethods.param.List,'Value',cSegmentation.paramSelected);

val = getfield(Segmethods.param.all,Segmethods.param.List{cSegmentation.paramSelected});
if val <0, tmp = '-';else  tmp = []; end
set(handles.ParameterExtractionListVal,'String',[tmp num2str(abs(val))]);
set(handles.DetectedCells,'String',num2str(nBasis));

% --- Executes during object creation, after setting all properties.
function ParameterExtractionListVal_CreateFcn(hObject, eventdata, handles)
global Segmethods;
global cSegmentation;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if ~exist('Segmethods.param','var'),
    initSegmentationMethods()
end

val = getfield(Segmethods.param.all,Segmethods.param.List{cSegmentation.paramSelected});
if val <0, tmp = '-';else  tmp = []; end
set(hObject,'String',[tmp num2str(abs(val))]);


function ParameterExtractionList_CreateFcn(hObject, eventdata, handles)
global Segmethods;
global cSegmentation;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if ~exist('Segmethods.param','var'),
    initSegmentationMethods()
end
set(hObject,'String',Segmethods.param.List,'Value',cSegmentation.paramSelected);


function DetectedCells_CreateFcn(hObject, eventdata, handles)
global Segmethods;
global cSegmentation;
global nGBasis
nGBasis = 0;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if ~exist('Segmethods.param','var'),
    initSegmentationMethods()
end
set(hObject,'String','0')

function SimilarityFtList_CreateFcn(hObject, eventdata, handles)
global Segmethods;
global cSegmentation;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
if ~exist('Segmethods.param','var'),
    initSegmentationMethods()
end
set(hObject,'String',Segmethods.ListSim,'Value',cSegmentation.SimilarSelected);


function MethodExtractionList_CreateFcn(hObject, eventdata, handles)
global Segmethods;
global cSegmentation;

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

if ~exist('Segmethods.param','var'),
    initSegmentationMethods()
end
set(hObject,'String',Segmethods.ListM,'Value',cSegmentation.MethodSelected);


% ----------------------------------------------------------------------------------------------------------------------
% Ending Cell Extraction Actions
% ----------------------------------------------------------------------------------------------------------------------

% ----------------------------------------------------------------------------------------------------------------------
% Beginning Display Actions
% ----------------------------------------------------------------------------------------------------------------------

function DisplayButton_display_Callback(hObject, eventdata, handles)
global data
global Coeff
global X
global Dict
global Coeff_single
global Dicts_single
global cNeurons
global displayFt
global displayParam
global GT
% displayFt.List = {'plotInputData','plotReconstructionMF',...
%     'plotReconstructionComponentExtraction','plotResultsMF',...
%     'plotResultsComponentExtraction','displaySingleCells','displayExploreCells'}


if ~isempty(data)
    if isempty(X)
        X = data;
    end
    switch lower(displayFt.List{displayFt.selected}),
        case 'plotinputdata'
            disp('Displaying input Data...')
            axes(handles.DisplaySequences);
            plotInputData(X,displayParam)
            disp('Displaying input data is done')
        case 'plotreconstructionmf'
            axes(handles.DisplaySequences);
            if ~isempty(Dict) && ~isempty(Coeff)
                disp('Displaying reconstructed data using Maxtrix Factorization result...')
                plotInputData(reshape(Dict*Coeff',size(X)),displayParam);
                disp('Displaying resulting data is done')
            else
                msgbox('Please run matrix factorization function before showing plotreconstructionMF.','Error','error','modal');
            end
        case 'plotreconstructioncomponentextraction'
            axes(handles.DisplaySequences);
            if ~isempty(Dicts_single) && ~isempty(Coeff_single)
                disp('Displaying reconstructed data using Extracting single components algorithm...')
               plotInputData(reshape(Dicts_single*Coeff_single',size(X)),displayParam);
                disp('Displaying resulting data is done')
            else
                msgbox('Please run Extraction of single components function before showing plotreconstructionCompenentExtraction.','Error','error','modal');
            end
        case 'plotresultsmf'
            if ~isempty(Dict) && ~isempty(Coeff)                
                disp('Displaying resuls from Matrix Factorization algorithm...')
                plotSequence(X,Coeff,Dict,displayParam);
                disp('Result displaying is done')
            else
                msgbox('Please run matrix factorization function before showing plotreconstructionMF.','Error','error','modal');
            end
        case 'plotresultscomponentextraction'   
            if ~isempty(Dicts_single) && ~isempty(Coeff_single)
                disp('Displaying resuls from Extracting single components algorithm...')
                plotSequence(X,Coeff_single,Dicts_single,displayParam);
                disp('Displaying resulting data is done')
            else
                msgbox('Please run Extraction of single components function before showing plotreconstructionCompenentExtraction.','Error','error','modal');
            end
        case 'displaysinglecells'
            if ~isempty(Dicts_single) && ~isempty(Coeff_single)
                disp('Exploring extracted cell ... done')
                axes(handles.DisplaySequences);
                displaySingleCellsIdentification(data,cNeurons);
                hold off
            else
                msgbox('Please run Extraction of single components function before showing plotreconstructionCompenentExtraction.','Error','error','modal');
            end
        case 'displayexplorecells'
            if ~isempty(Dicts_single) && ~isempty(Coeff_single)
                disp('Exploring extracted cells together with temporal evolution ... done')
                displayExploreCells(data,Coeff_single,cNeurons,displayParam);
            else
                msgbox('Please run Extraction of single components function before showing plotreconstructionCompenentExtraction.','Error','error','modal');
            end
        case 'displaytemporalevolutioncells'
            if ~isempty(Dicts_single) && ~isempty(Coeff_single)
                if ~isempty(GT)                    
                    [Ucells_GT_reord] = reorderGTMatrix(GT.extractedCells,Dicts_single,GT.Ucells,Coeff_single,0.5);
                    displayTemporalEvolutionCells(data,Ucells_GT_reord,Coeff_single,cNeurons)
                else
                    msgbox('This function is only used for comparing the results with the ground truth.','Error','error','modal');
                end
            else                
                msgbox('Please run Extraction of single components function before showing displayTemporalEvolutionCells.','Error','error','modal');
            end
        otherwise
            msgbox('Error on showing the results.','Error','error','modal');
    end
else
    msgbox('Please load a sequence before setting the number of basis functions.','Error','error','modal');
end

function updateDisplay(handles)
global displayParam
set(handles.EndFrameTag,'String',num2str(displayParam.endF));
set(handles.StepFrame,'String',num2str(displayParam.step));
set(handles.InitFrameTag,'String',num2str(displayParam.initF));

function DisplayFunctionsTag_Callback(hObject, eventdata, handles)
global displayFt
displayFt.selected = get(hObject,'Value'); 

function EndFrameTag_Callback(hObject, eventdata, handles)
global displayParam
global X

value = get(hObject,'String');
value = str2double(value);
value  = min(max(displayParam.initF,value),size(X,3));
displayParam.endF = value;
displayParam.step = min(max(displayParam.step,1),displayParam.endF-displayParam.initF+1);
updateDisplay(handles)

function ParametersDisplay_Callback(hObject, eventdata, handles)


function InitFrameTag_Callback(hObject, eventdata, handles)
global displayParam


value = get(hObject,'String');
value = str2double(value);
value  = max(min(displayParam.endF,value),1);
displayParam.initF = value;
displayParam.step = min(max(displayParam.step,1),displayParam.endF-displayParam.initF+1);
updateDisplay(handles)

function StepFrame_Callback(hObject, eventdata, handles)
global displayParam
value = get(hObject,'String');
value = str2double(value);
value  = min(max(value,1),displayParam.endF-displayParam.initF+1);
displayParam.step = value;
updateDisplay(handles)


% --- Executes during object creation, after setting all properties.
function InitFrameTag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','1')

function ParametersDisplay_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','Not implemented');


function EndFrameTag_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','1')

function DisplayFunctionsTag_CreateFcn(hObject, eventdata, handles)
global displayFt
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
displayFt.selected = 1; 
global displayParam
displayParam.initF = 1;
displayParam.endF  = 1;
displayParam.step  = 1;


% filesDir = dir([pwd '/display/*.m']);
% for i = 1:length(filesDir)
%     displayFt.List{i} = filesDir(i).name(1:end-2);
% end
displayFt.List = {'plotInputData','plotReconstructionMF',...
    'plotReconstructionComponentExtraction','plotResultsMF',...
    'plotResultsComponentExtraction','displaySingleCells','displayExploreCells','displayTemporalEvolutionCells'};
set(hObject,'String',displayFt.List,'Value',1);


% --- Executes during object creation, after setting all properties.
function StepFrame_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'String','1')



% ----------------------------------------------------------------------------------------------------------------------
% Ending Display Actions
% ----------------------------------------------------------------------------------------------------------------------



function initMethods()


global MFmethods;
global cMethod;
cMethod.MethodSelected = 1;
cMethod.subMethodSelected = 1;
cMethod.paramSelected = 1;

MFmethods.List = get_methods(pwd);
for i = 1:length(MFmethods.List)
    MFmethods.Method{i} = MF_methods(MFmethods.List{i});
end

function initSegmentationMethods()

global Segmethods;
global cSegmentation;
Segmethods.ListM = {'wavelet_denoising','image_smoothing','image_smoothing_powerlaw'};
Segmethods.ListSim = {'jointDistance','MahalanobisDistance','crosscorrelation'};
cSegmentation.MethodSelected = 1;
cSegmentation.SimilarSelected = 3;
cSegmentation.paramSelected = 1;
Segmethods.param = initial_parameter_Segmentation();


function initPreprocessingStep()
global paramgui
paramgui.prepro.List = {'none','wavelet transform','delta F','delta F + wavelet'};
paramgui.prepro.ListParam = {'BorderSize','quantile','halfWinMedian','tempSigma','spaSigma','wavScaleSpace','wavScaleTime','stepF0'};
paramgui.prepro.param = initial_parameter_preprocessing();
paramgui.prepro.methodSelected = 1;
paramgui.prepro.paramSelected = 1;


function initSequenceMethods()

global sequenceFts
sequenceFts.List = {'non','image resize','inverted intensities','wavelet denoising'};
sequenceFts.selected = 1;
sequenceFts.param = {' ','1','false','1'};

function [value] = fieldname_GUI(nameField)

value = 0;
switch lower(nameField)
    case 'data'
        value = 1;
    case 'coeff'
        value = 1;
    case 'x'
        value = 1;
    case 'dict'
        value = 1;
    case 'coeff_single'
        value = 1;
    case 'dicts_single'
        value = 1;
    case 'cneurons'
        value = 1;
    case 'displayft'
        value = 1;
    case 'segmethods'
        value = 1;
    case 'csegmentation'
        value = 1;
    case 'nbasis'
        value = 1;
    case 'mfmethods'
        value = 1;
    case 'cmethod'
        value = 1;
    case 'reconsdict'
        value = 1;
    case 'reconsdict_single'
        value = 1;
    case 'ngbasis'
        value = 1;
    case 'preproc_params'
        value =1;
    case 'articialseq'
        value =1;
    case 'cartificialseq'
        value =1;
    otherwise
        value = 0;
end


function List = get_Variablelist(subpart)

if nargin <1
    subpart = 'complete';
elseif nargin >1,
    error('Incorrect number of input variables at get_VariableList');
end

switch lower(subpart)
    case 'complete'
        List = {'data','data_orig','Coeff','X','Dict','Coeff_single','Dicts_single',...
            'cNeurons','displayFt','Segmethods','cSegmentation','nBasis','MFmethods',...
            'cMethod','nGBasis','ReconsDict_single','ReconsDict','preproc_params','sequenceFts',...
            'sequenceName','displayParam','paramgui'};
    case 'mf'
        List = {'Coeff','Dict','nBasis','MFmethods','cMethod','nGBasis'};
    case 'segmentation'
        List = {'Coeff_single','Dicts_single','cNeurons','Segmethods','cSegmentation','nGBasis','ReconsDict_single','paramgui'};
    otherwise 
        error('Incorrect input variable for extracting List of Variables')
end


% --- Executes on selection change in enhacementMenu.
function enhacementMenu_Callback(hObject, eventdata, handles)
% hObject    handle to enhacementMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns enhacementMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from enhacementMenu
global paramgui;
num = get(hObject, 'Value');
paramgui.prepro.methodSelected = num;
paramgui.prepro.param.method = num;
updateEnhancement(handles)

% --- Executes during object creation, after setting all properties.
function enhacementMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to enhacementMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global paramgui
if ~isfield(paramgui,'prepro')
    initialPreprocessingStep()
end
set(hObject,'String',paramgui.prepro.List,'Value',paramgui.prepro.methodSelected);



function valueEnhacement_Callback(hObject, eventdata, handles)
% hObject    handle to valueEnhacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of valueEnhacement as text
%        str2double(get(hObject,'String')) returns contents of valueEnhacement as a double
global paramgui
value = get(hObject,'String');
value = str2double(value);
paramgui.prepro.param = set_parameters_ArtificialSeq(paramgui.prepro.param,{paramgui.prepro.ListParam{paramgui.prepro.paramSelected},value});
updateEnhancement(handles)

% --- Executes during object creation, after setting all properties.
function valueEnhacement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to valueEnhacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global paramgui
if ~isfield(paramgui,'prepro')
    initPreprocessingStep()
end
val = getfield(paramgui.prepro.param,paramgui.prepro.ListParam{paramgui.prepro.paramSelected});
if val <0, tmp = '-';else  tmp = []; end
set(hObject,'String',[tmp num2str(abs(val))]);


% --- Executes on selection change in paramEnhacement.
function paramEnhacement_Callback(hObject, eventdata, handles)
% hObject    handle to paramEnhacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns paramEnhacement contents as cell array
%        contents{get(hObject,'Value')} returns selected item from paramEnhacement
global paramgui;
num = get(hObject, 'Value');
paramgui.prepro.paramSelected = num;
updateEnhancement(handles)

% --- Executes during object creation, after setting all properties.
function paramEnhacement_CreateFcn(hObject, eventdata, handles)
% hObject    handle to paramEnhacement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
global paramgui
if ~isfield(paramgui,'prepro')
    initPreprocessingStep()
end
    
set(hObject,'String',paramgui.prepro.ListParam,'Value',paramgui.prepro.paramSelected);

function updateEnhancement(handles)
global paramgui
val = getfield(paramgui.prepro.param,paramgui.prepro.ListParam{paramgui.prepro.paramSelected});
set(handles.valueEnhacement,'String',num2str(val))
set(handles.paramEnhacement,'Value',paramgui.prepro.paramSelected)
set(handles.enhacementMenu,'Value',paramgui.prepro.methodSelected)


% --- Executes on button press in runEnhancement.
function runEnhancement_Callback(hObject, eventdata, handles)
% hObject    handle to runEnhancement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global paramgui
global data

params = paramgui.prepro.param;

if ~isempty(data)
    data = paramgui.data(params.BorderSize+1:end-params.BorderSize,params.BorderSize+1:end-params.BorderSize,:);
    data_tmp = extractingNormalizedSeq(double(data),params);
    data = data_tmp.data_denoise;
else
    msgbox('Please load a sequence before setting the number of basis functions.','Error','error','modal');
end


