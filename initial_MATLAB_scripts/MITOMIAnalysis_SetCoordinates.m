function varargout = MITOMIAnalysis_SetCoordinates(varargin)
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

% Last Modified by GUIDE v2.5 01-Feb-2017 18:00:22

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MITOMIAnalysis_SetCoordinates_OpeningFcn, ...
                   'gui_OutputFcn',  @MITOMIAnalysis_SetCoordinates_OutputFcn, ...
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
function MITOMIAnalysis_SetCoordinates_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for MITOMIAnalysis_GUI
handles.output = hObject;
set(handles.uipanel_scroll,'UserData',varargin{1});
set(handles.text_directions,'String',sprintf('1.Navigate to one of the four corners of the array \n\n2. Zoom into the spot to easily see the spot \n\n3. Adjust Gamma to make the background barely visible \n4. Lock the image \n5. Press the ''Set Coordinate'' button \n\n6. Click on three points on the circumference of the spot \n7. Unlock image and continue to the next corner of the array'))
guidata(hObject, handles);

global Log

scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
scrollImage = imshow(get(handles.uipanel_scroll,'UserData'),'parent',scrollAxes);
Log.ManipulationAPI = imscrollpanel(handles.uipanel_scroll,scrollImage); 
if Log.vertex==0
    Log.ImageManipulationSurface=[];
    Log.ImageManipulationBackground=[];
end

% UIWAIT makes MITOMIAnalysis_GUI wait for user response (see UIRESUME)
uiwait(handles.figure_manipulation);



% --- Outputs from this function are returned to the command line.
function MITOMIAnalysis_SetCoordinates_OutputFcn(hObject, eventdata, handles) 

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(hObject, eventdata, handles)
delete(handles.figure_manipulation)

% --- Executes on slider movement.
function slider_gamma_Callback(hObject, eventdata, handles)

global Log
Log.tempgamma=get(hObject,'Value');
API=iptgetapi(Log.ManipulationAPI);
display=imadjust(get(handles.uipanel_scroll,'UserData'),[],[],Log.tempgamma);
API.replaceImage(display,'PreserveView',true);


% --- Executes during object creation, after setting all properties.
function slider_gamma_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_zoom_Callback(hObject, eventdata, handles)

global Log
API=iptgetapi(Log.ManipulationAPI);
API.setMagnification(2^(get(hObject,'Value')-1));


% --- Executes during object creation, after setting all properties.
function slider_zoom_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton_lock.
function pushbutton_lock_Callback(hObject, eventdata, handles)

global Log

set(hObject,'UserData',get(hObject,'UserData')+1);
if rem(get(hObject,'UserData'),2)
    API=iptgetapi(Log.ManipulationAPI);
    dimImage=size(get(handles.uipanel_scroll,'UserData'));
    Log.currView=API.getVisibleImageRect();
    Log.magView=API.getMagnification();
    if round(Log.currView(1)+Log.currView(3))> dimImage(2)
        Log.currView(1)=floor(dimImage(2)-Log.currView(3));
    end
    if round(Log.currView(2)+Log.currView(4))> dimImage(1)
        Log.currView(2)=floor(dimImage(1)-Log.currView(4));
    end
    fullImage=get(handles.uipanel_scroll,'UserData');
    delete(handles.uipanel_scroll.Children)
    scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
    subImage=fullImage(round(Log.currView(2):round(Log.currView(2)+Log.currView(4))),round(Log.currView(1)):round(Log.currView(1)+Log.currView(3)));
    display=imadjust(subImage,[],[],Log.tempgamma);
    imshow(display,'parent',scrollAxes);
    set(handles.pushbutton_setcoor,{'Enable','ForegroundColor'},{'on','black'});
    set(handles.slider_gamma,'Enable','inactive');
    set(handles.slider_zoom,'Enable','inactive');
    set(handles.pushbutton_lock,'String','Unlock Image to Navigate')
else
    set(handles.pushbutton_setcoor,{'Enable','ForegroundColor','String'},{'inactive',[0.8 0.8 0.8],sprintf('Set Coordinate #%i',get(handles.pushbutton_setcoor,'UserData')+1)});
    delete(handles.uipanel_scroll.Children)
    scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
    display=imadjust(get(handles.uipanel_scroll,'UserData'),[],[],Log.tempgamma);
    scrollImage = imshow(display,'parent',scrollAxes);
    Log.ManipulationAPI = imscrollpanel(handles.uipanel_scroll,scrollImage);
    API=iptgetapi(Log.ManipulationAPI);
    API.setMagnification(Log.magView);
    API.setVisibleLocation(Log.currView(1:2));
    set(handles.slider_gamma,'Enable','on');
    set(handles.slider_zoom,'Enable','on');
    set(handles.pushbutton_lock,'String','Lock Image to Set Coordinates')

end
guidata(hObject,handles)

% --- Executes on button press in pushbutton_setcoor.
function pushbutton_setcoor_Callback(hObject, eventdata, handles)
global Log

set(hObject,{'UserData'},{get(hObject,'UserData')+1}); 
set(handles.pushbutton_lock,'Enable','inactive','ForegroundColor',[0.8 0.8 0.8]);
try
    sample=ginputax(handles.uipanel_scroll,3);
catch
    assert(false,'MITOMIAnalysis:ImageManipulation:axisError','User clicked off the image when setting coordinates')
    return
end
Log.vertex=Log.vertex+1;
Log.Coor{Log.vertex}=[sample(:,1)+Log.currView(1),sample(:,2)+Log.currView(2)];
Log.gamma(Log.vertex)=get(handles.slider_gamma,'Value');

if get(hObject,'UserData')>3
    if isempty(Log.ImageManipulationSurface)
        Log.ImageManipulationSurface='Passed';
    else
        Log.ImageManipulationBackground='Passed';
    end
    delete(handles.figure_manipulation)
    return
end
set(handles.pushbutton_lock,'Enable','on','ForegroundColor','black');
set(handles.pushbutton_setcoor,{'Enable','ForegroundColor','String'},{'inactive',[0.8 0.8 0.8],sprintf('Coordinate #%i Set !',get(hObject,'UserData'))});
guidata(hObject,handles);
