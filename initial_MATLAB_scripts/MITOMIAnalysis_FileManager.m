function varargout = MITOMIAnalysis_FileManager(varargin)
% MITOMIANALYSIS_FILEMANAGER MATLAB code for MITOMIAnalysis_FileManager.fig
%      MITOMIANALYSIS_FILEMANAGER, by itself, creates a new MITOMIANALYSIS_FILEMANAGER or raises the existing
%      singleton*.
%
%      H = MITOMIANALYSIS_FILEMANAGER returns the handle to a new MITOMIANALYSIS_FILEMANAGER or the handle to
%      the existing singleton*.
%
%      MITOMIANALYSIS_FILEMANAGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MITOMIANALYSIS_FILEMANAGER.M with the given input arguments.
%
%      MITOMIANALYSIS_FILEMANAGER('Property','Value',...) creates a new MITOMIANALYSIS_FILEMANAGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MITOMIAnalysis_FileManager_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MITOMIAnalysis_FileManager_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help MITOMIAnalysis_FileManager

% Last Modified by GUIDE v2.5 30-Jan-2017 12:37:04

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MITOMIAnalysis_FileManager_OpeningFcn, ...
                   'gui_OutputFcn',  @MITOMIAnalysis_FileManager_OutputFcn, ...
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


% --- Executes just before MITOMIAnalysis_FileManager is made visible.
function MITOMIAnalysis_FileManager_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MITOMIAnalysis_FileManager (see VARARGIN)

% Choose default command line output for MITOMIAnalysis_FileManager
handles.output = hObject;
global Log
Log.FileManager=[];
set(handles.listbox_files,'String',cellstr(Log.nameFrames))
set(handles.text_order,'String',sprintf('Equilbrium:\n\n1. Background\n    (if included)\n2. Surface Molecule\n    (traditionally protein)\n3. Solubilized Molecule \n    (traditionally DNA) \n\nDissociation:\n\n1. Background\n    (if included)\n2. Surface Molecule\n3. Solubilized Molecule #1\n :\nn. Solubilized Molecule #n \n    (in temporal sequence)'))
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes MITOMIAnalysis_FileManager wait for user response (see UIRESUME)
uiwait(handles.figure_filemanager);


% --- Outputs from this function are returned to the command line.
function [] = MITOMIAnalysis_FileManager_OutputFcn(hObject, eventdata, handles) 

% --- Executes on selection change in listbox_files.
function listbox_files_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listbox_files_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_moveup.
function pushbutton_moveup_Callback(hObject, eventdata, handles)

selection=get(handles.listbox_files,'Value');
if min(selection) > 1
    holderList=get(handles.listbox_files,'String');
    swapVal=holderList(min(selection)-1);
    holderList(selection-1)=holderList(selection);
    holderList(max(selection))=swapVal;
    set(handles.listbox_files,'String',holderList);
    set(handles.listbox_files,'Value',selection-1);
    guidata(hObject,handles);
end



% --- Executes on button press in pushbutton_movedown.
function pushbutton_movedown_Callback(hObject, eventdata, handles)

selection=get(handles.listbox_files,'Value');
if max(selection) < length(get(handles.listbox_files,'String'))
    holderList=get(handles.listbox_files,'String');
    swapVal=holderList(max(selection)+1);
    holderList(selection+1)=holderList(selection);
    holderList(min(selection))=swapVal;
    set(handles.listbox_files,'String',holderList);
    set(handles.listbox_files,'Value',selection+1);
    guidata(hObject,handles);
end



% --- Executes on button press in pushbutton_Continue.
function pushbutton_Continue_Callback(hObject, eventdata, handles)
global Log
Log.FileManager='Passed';
Log.nameFrames=get(handles.listbox_files,'String');
delete(handles.figure_filemanager)


% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)

delete(handles.figure_filemanager)
