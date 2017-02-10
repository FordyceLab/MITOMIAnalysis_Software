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

% Last Modified by GUIDE v2.5 01-Feb-2017 17:58:34

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

%Set initial control text and pass image into object
handles.output = hObject;
set(handles.uipanel_scroll,'UserData',varargin{1});
set(handles.text_directions,'String',sprintf('1.Navigate to one of the four corners of the array \n\n2. Zoom into the spot to easily see the spot \n\n3. Adjust Gamma to make the background barely visible \n4. Lock the image \n5. Press the ''Set Coordinate'' button \n\n6. Click on three points on the circumference of the spot \n7. Unlock image and continue to the next corner of the array'))
guidata(hObject, handles);

global Log

%Generate initial image and grab controls
scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
scrollImage = imshow(get(handles.uipanel_scroll,'UserData'),'parent',scrollAxes);
Log.ManipulationAPI = imscrollpanel(handles.uipanel_scroll,scrollImage); 
Log.TempGamma=1;

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

%Grab image controls, current gamma and update image
Log.TempGamma=get(hObject,'Value');
API=iptgetapi(Log.ManipulationAPI);
display=imadjust(get(handles.uipanel_scroll,'UserData'),[],[],Log.TempGamma);
API.replaceImage(display,'PreserveView',true);


% --- Executes during object creation, after setting all properties.
function slider_gamma_CreateFcn(hObject, eventdata, handles)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_zoom_Callback(hObject, eventdata, handles)

global Log

%Grab image controls and update magnification
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

%Toggle lock state
set(hObject,'UserData',get(hObject,'UserData')+1);

if rem(get(hObject,'UserData'),2) %Lock GUI
    
    %Grab image controls and states
    API=iptgetapi(Log.ManipulationAPI);
    dimImage=size(get(handles.uipanel_scroll,'UserData'));
    Log.CurrView=API.getVisibleImageRect();
    Log.MagView=API.getMagnification();
    
    %Prevent image coordinates from exceeding image dimensions
    if round(Log.CurrView(1)+Log.CurrView(3))> dimImage(2)
        Log.CurrView(1)=floor(dimImage(2)-Log.CurrView(3));
    end
    
    if round(Log.CurrView(2)+Log.CurrView(4))> dimImage(1)
        Log.CurrView(2)=floor(dimImage(1)-Log.CurrView(4));
    end
    
    %Replace full image with subimage
    fullImage=get(handles.uipanel_scroll,'UserData');
    delete(handles.uipanel_scroll.Children)
    scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
    subImage=fullImage(round(Log.CurrView(2):round(Log.CurrView(2)+Log.CurrView(4))),round(Log.CurrView(1)):round(Log.CurrView(1)+Log.CurrView(3)));
    display=imadjust(subImage,[],[],Log.TempGamma);
    imshow(display,'parent',scrollAxes);
    
    %Update GUI controls
    set(handles.pushbutton_setcoor,{'Enable','ForegroundColor'},{'on','black'});
    set(handles.slider_gamma,'Enable','inactive');
    set(handles.slider_zoom,'Enable','inactive');
    set(handles.pushbutton_lock,'String','Unlock Image to Navigate')
    
else %Unlock GUI
    
    %Replace subimage with full image
    delete(handles.uipanel_scroll.Children)
    scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
    display=imadjust(get(handles.uipanel_scroll,'UserData'),[],[],Log.TempGamma);
    scrollImage = imshow(display,'parent',scrollAxes);
    
    %Grab new image controls and reset states
    Log.ManipulationAPI = imscrollpanel(handles.uipanel_scroll,scrollImage);
    API=iptgetapi(Log.ManipulationAPI);
    API.setMagnification(Log.MagView);
    API.setVisibleLocation(Log.CurrView(1:2));
    
    %Update GUI controls
    set(handles.pushbutton_setcoor,{'Enable','ForegroundColor','String'},{'inactive',[0.8 0.8 0.8],sprintf('Set Coordinate #%i',get(handles.pushbutton_setcoor,'UserData')+1)});
    set(handles.slider_gamma,'Enable','on');
    set(handles.slider_zoom,'Enable','on');
    set(handles.pushbutton_lock,'String','Lock Image to Set Coordinates')

end
guidata(hObject,handles)

% --- Executes on button press in pushbutton_setcoor.
function pushbutton_setcoor_Callback(hObject, eventdata, handles)

global Log

%Lock GUI operations
set(hObject,{'UserData'},{get(hObject,'UserData')+1}); 
set(handles.pushbutton_lock,'Enable','inactive','ForegroundColor',[0.8 0.8 0.8]);

try
    
    %Collect points from circumference
    [sample]=ginputax(handles.uipanel_scroll,3);
    
    %Find center of circle and approximate radius
    p1=[sample(1,:) 1]; 
    p2=[sample(2,:) 1]; 
    p3=[sample(3,:) 1];
    t = p2-p1; 
    u = p3-p1; 
    v = p3-p2;
    w = cross(t,u);
    t2 = sum(t.^2); 
    u2 = sum(u.^2); 
    w2 = sum(w.^2);
    c = p1+(t2*sum(u.*v)*u-u2*sum(t.*v)*t)/(2*w2);
    RadiusSample = (sqrt(t2*u2*sum(v.^2)/w2))/2;

catch
    assert(false,'MITOMIAnalysis:ImageManipulation:axisError','User clicked off the image when setting coordinates')
    return %%%CONFIRM THIS RETURN CLOSES GUI ON ERROR%%%
end

%Pass to Log
Log.Vertex=Log.Vertex+1;
Log.RadiusSample(Log.Vertex)=RadiusSample;
Log.CoorList(Log.Vertex,:)=[c(1,1)+Log.CurrView(1),c(1,2)+Log.CurrView(2)];
Log.Gamma(Log.Vertex)=Log.TempGamma;

if get(hObject,'UserData')>3 %User has selected 4 vertices
    
    if isempty(Log.SetCoordinatesSurface)  %User's first pass on func
        
        %Interpolate array for button positions
        vertices=round(sortrows(Log.CoorList(1:4,:),2));
        TopRowX=sort(linspace(vertices(1,1),vertices(2,1),Log.Cols));
        BotRowX=sort(linspace(vertices(3,1),vertices(4,1),Log.Cols));
        TopRowY=interp1(vertices(1:2,1),vertices(1:2,2),TopRowX);
        BotRowY=interp1(vertices(3:4,1),vertices(3:4,2),BotRowX);
        
        for j=1:Log.Cols
            ColYVal=sort(linspace(TopRowY(j),BotRowY(j),Log.Rows));
            ColXVal=interp1([TopRowY(j) BotRowY(j)],[TopRowX(j) BotRowX(j)],ColYVal);
            Log.CoorButtons=round(vertcat(Log.CoorButtons,[ColXVal' ColYVal']));
        end
        
        Log.ApproxButtonRadius=round(mean(Log.RadiusSample(1:4)));
        Log.ApproxGammaSurface=mean(Log.Gamma(1:4));
        Log.SetCoordinatesSurface='Passed';
        
    else %User has completed function once before
        
        %Interpolate array for background positions
        vertices=round(sortrows(Log.CoorList(5:8,:),2));
        TopRowX=sort(linspace(vertices(1,1),vertices(2,1),Log.Cols));
        BotRowX=sort(linspace(vertices(3,1),vertices(4,1),Log.Cols));
        TopRowY=interp1(vertices(1:2,1),vertices(1:2,2),TopRowX);
        BotRowY=interp1(vertices(3:4,1),vertices(3:4,2),BotRowX);
        for j=1:Log.Cols
            ColYVal=sort(linspace(TopRowY(j),BotRowY(j),Log.Rows));
            ColXVal=interp1([TopRowY(j) BotRowY(j)],[TopRowX(j) BotRowX(j)],ColYVal);
            Log.CoorBackground=round(vertcat(Log.CoorBackground,[ColXVal' ColYVal']));
        end
        
        Log.ApproxBackgroundRadius=round(mean(Log.RadiusSample(5:8)));
        Log.SetCoordinatesBackground='Passed';
    end
    
    %Close GUI
    delete(handles.figure_manipulation)
    return
    
end

%Reset GUI for user
set(handles.pushbutton_lock,'Enable','on','ForegroundColor','black');
set(handles.pushbutton_setcoor,{'Enable','ForegroundColor','String'},{'inactive',[0.8 0.8 0.8],sprintf('Coordinate #%i Set !',get(hObject,'UserData'))});
guidata(hObject,handles);
