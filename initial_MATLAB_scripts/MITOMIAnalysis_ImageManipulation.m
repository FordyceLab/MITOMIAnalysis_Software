function varargout = MITOMIAnalysis_ImageManipulation(varargin)
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

% Last Modified by GUIDE v2.5 03-Feb-2017 16:44:57

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MITOMIAnalysis_ImageManipulation_OpeningFcn, ...
                   'gui_OutputFcn',  @MITOMIAnalysis_ImageManipulation_OutputFcn, ...
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
function MITOMIAnalysis_ImageManipulation_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for MITOMIAnalysis_GUI
global Log
Log.Manipulation=[];
handles.output = hObject;
set(handles.uipanel_image,'UserData',mat2gray(varargin{1}));
set(handles.figure_manipulation,'UserData',varargin{2});
set(handles.text_directions,'String',sprintf('Directions:\nReorient the image such that chamber #1 is on the top left and chamber #2 is directly below that.\n\nAlso, set the dimensions of the array in the text boxes below'))
guidata(hObject, handles);

scrollAxes = axes('parent',handles.uipanel_image,'position',[0 0 1 1],'Units','pixels');
imshow(get(handles.uipanel_image,'UserData'),'parent',scrollAxes);

% UIWAIT makes MITOMIAnalysis_GUI wait for user response (see UIRESUME)
uiwait(handles.figure_manipulation);

% --- Outputs from this function are returned to the command line.
function MITOMIAnalysis_ImageManipulation_OutputFcn(hObject, eventdata, handles) 

function pushbutton_continue_Callback(hObject, eventdata, handles)

global Log

Log.ImageHolder=get(handles.figure_manipulation,'UserData');
Log.Rows=str2num(get(handles.edit_col,'String'));
Log.Cols=str2num(get(handles.edit_row,'String'));

try
    assert(Log.Rows>0 && rem(Log.Rows,1),'MITOMIAnalysis:ImageManipulation:rowValue','Row value declared was not an integer greater than 0');
    assert(Log.Cols>0 && rem(Log.Cols,1),'MITOIMAnalysis:ImageManipulation:colValue','Column value declared was not an integer greater than 0');
catch
    delete(handles.figure_manipulation)
end

Log.ImageManipulation='Passed';

delete(handles.figure_manipulation)


function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure_manipulation)

function slider_gamma_Callback(hObject, eventdata, handles)

display=imadjust(get(handles.uipanel_image,'UserData'),[],[],get(hObject,'Value'));
scrollAxes = axes('parent',handles.uipanel_image,'position',[0 0 1 1],'Units','pixels');
imshow(display,'parent',scrollAxes);

function pushbutton_lr_Callback(hObject, eventdata, handles)

    impreview=get(handles.uipanel_image,'UserData');
    images=get(handles.figure_manipulation,'UserData');
    impreview=fliplr(impreview);
    display=imadjust(impreview,[],[],get(handles.slider_gamma,'Value'));
    delete(handles.uipanel_image.Children)
    scrollAxes = axes('parent',handles.uipanel_image,'position',[0 0 1 1],'Units','pixels');
    imshow(display,'parent',scrollAxes);
    images= structfun(@fliplr, images, 'UniformOutput', false);
    set(handles.uipanel_image,'UserData',impreview);
    set(handles.figure_manipulation,'UserData',images)
    guidata(hObject,handles)

function pushbutton_tb_Callback(hObject, eventdata, handles)

    impreview=get(handles.uipanel_image,'UserData');
    images=get(handles.figure_manipulation,'UserData');
    impreview=flipud(impreview);   
    display=imadjust(impreview,[],[],get(handles.slider_gamma,'Value'));
    delete(handles.uipanel_image.Children)
    scrollAxes = axes('parent',handles.uipanel_image,'position',[0 0 1 1],'Units','pixels');
    imshow(display,'parent',scrollAxes);
    images= structfun(@flipud, images, 'UniformOutput', false);
    set(handles.uipanel_image,'UserData',impreview);
    set(handles.figure_manipulation,'UserData',images)
    guidata(hObject,handles)


function pushbutton_p90_Callback(hObject, eventdata, handles)

    impreview=get(handles.uipanel_image,'UserData');
    images=get(handles.figure_manipulation,'UserData');
    impreview=imrotate(impreview,90);   
    display=imadjust(impreview,[],[],get(handles.slider_gamma,'Value'));
    delete(handles.uipanel_image.Children)
    scrollAxes = axes('parent',handles.uipanel_image,'position',[0 0 1 1],'Units','pixels');
    imshow(display,'parent',scrollAxes);
    images= structfun(@(x) (imrotate(x,90)), images, 'UniformOutput', false);
    set(handles.uipanel_image,'UserData',impreview);
    set(handles.figure_manipulation,'UserData',images)
    guidata(hObject,handles)

function pushbutton_n90_Callback(hObject, eventdata, handles)

    impreview=get(handles.uipanel_image,'UserData');
    images=get(handles.figure_manipulation,'UserData');
    impreview=imrotate(impreview,-90);    
    display=imadjust(impreview,[],[],get(handles.slider_gamma,'Value'));
    delete(handles.uipanel_image.Children)
    scrollAxes = axes('parent',handles.uipanel_image,'position',[0 0 1 1],'Units','pixels');
    imshow(display,'parent',scrollAxes);
    images= structfun(@(x) (imrotate(x,-90)), images, 'UniformOutput', false);
    set(handles.uipanel_image,'UserData',impreview);
    set(handles.figure_manipulation,'UserData',images)
    guidata(hObject,handles)

function edit_row_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function edit_col_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function slider_gamma_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function edit_row_Callback(hObject,eventdata,handles)

function edit_col_Callback(hObject,eventdata,handles)
