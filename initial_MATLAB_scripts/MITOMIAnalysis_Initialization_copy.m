function varargout = MITOMIAnalysis_Initialization(varargin)
% MITOMIANALYSIS_GUI MATLAB code for MITOMIAnalysis_GUI.fig
%      MITOMIANALYSIS_GUI, by itself, creates a new MITOMIANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = MITOMIANALYSIS_GUI returns the handle to a new MITOMIANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      MITOMIANALYSIS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MITOMIANALYSIS_GUI.M with the given input arguments.
%
%      MITOMIANALYSIS_GUI('Property','Value',...) creates a new MITOMIANALYSIS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MITOMIAnalysis_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MITOMIAnalysis_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MITOMIAnalysis_GUI

% Last Modified by GUIDE v2.5 29-Jan-2017 09:41:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MITOMIAnalysis_Initialization_OpeningFcn, ...
                   'gui_OutputFcn',  @MITOMIAnalysis_Initialization_OutputFcn, ...
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


% --- Executes just before MITOMIAnalysis_GUI is made visible.
function MITOMIAnalysis_Initialization_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for MITOMIAnalysis_GUI
handles.output = hObject;
set(handles.edit_output,'String',['MITOMIAnalysis_' datestr(now,30)]);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MITOMIAnalysis_GUI wait for user response (see UIRESUME)
uiwait(handles.figure_initialization);



% --- Outputs from this function are returned to the command line.
function [] = MITOMIAnalysis_Initialization_OutputFcn(hObject, eventdata, handles) 

% --- Executes on button press in pushbutton_navdir.
function pushbutton_navdir_Callback(hObject, eventdata, handles)

%Identify image directory, number of images, image names, and enable 
%'Continue' button
try
    dirname=uigetdir('','Select Directory for Analysis');
    cd(dirname);
    holderImages=[dir('*.tif');dir('*.jpg');dir('*.png')];
    nameImages={holderImages};
    numImages=size(holderImages);
    set(handles.text_numframe,'String',num2str(numImages(1)));
    set(handles.text_directory,'String',dirname);
    set(handles.text_directory,'UserData',{nameImages{1}.name;});
    set(handles.pushbutton_continue,'Enable','on','ForegroundColor','black');
    guidata(hObject,handles);
catch
    delete(handles.figure_initialization);
end

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)

%Close GUI, do not save any information
delete(handles.figure_initialization);


% --- Executes on button press in pushbutton_continue.
function pushbutton_continue_Callback(hObject, eventdata, handles)

%Pass data to global workspace for remaining MITOMIAnalysis functions
global Log
Log.dir=get(handles.text_directory,'String');
Log.type=get(get(handles.uibutton_analysisgroup,'SelectedObject'),'String');
Log.numFrames=str2double(get(handles.text_numframe,'String'));
Log.nameFrames=get(handles.text_directory,'UserData')';
Log.background=get(get(handles.uibutton_backgroundgroup,'SelectedObject'),'String');
Log.lastUpdated=20170129;
Log.nameOutput=get(handles.edit_output,'String');

%Close GUI
delete(handles.figure_initialization);



function edit_output_Callback(hObject, eventdata, handles)
% hObject    handle to edit_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_output as text
%        str2double(get(hObject,'String')) returns contents of edit_output as a double


% --- Executes during object creation, after setting all properties.
function edit_output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
