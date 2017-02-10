function varargout = MITOMIAnalysis_UserEdit(varargin)
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

% Last Modified by GUIDE v2.5 10-Feb-2017 13:45:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MITOMIAnalysis_UserEdit_OpeningFcn, ...
                   'gui_OutputFcn',  @MITOMIAnalysis_UserEdit_OutputFcn, ...
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
function MITOMIAnalysis_UserEdit_OpeningFcn(hObject, eventdata, handles, varargin)

%Set initial control text and pass images and data into objects
Image=varargin{1};
Data=varargin{2};
set(handles.uipanel_scroll,'UserData',Image);
set(handles.figure_manipulation,'UserData',Data);
set(handles.text_directions,'String',sprintf('1.Navigate to one of the four corners of the array \n\n2. Zoom into the spot to easily see the spot \n\n3. Adjust Gamma to make the background barely visible \n4. Lock the image \n5. Press the ''Set Coordinate'' button \n\n6. Click on three points on the circumference of the spot \n7. Unlock image and continue to the next corner of the array'))
set(handles.pushbutton_update,'UserData',0);
set(handles.uipanel_direction,'UserData',0);
set(handles.pushbutton_continue,'UserData',cell(1));
guidata(hObject, handles);

global Log

%Generate initial image and grab controls
AutofindOnImage(Image.Surface,Data.AutofindButtons,Data.ButtonsXCoor,Data.ButtonsYCoor,Data.ButtonsRadius,hObject,handles)
scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
scrollImage = imshow(get(handles.pushbutton_continue,'UserData'),'parent',scrollAxes);
Log.ManipulationAPI = imscrollpanel(handles.uipanel_scroll,scrollImage); 


% UIWAIT makes MITOMIAnalysis_GUI wait for user response (see UIRESUME)
uiwait(handles.figure_manipulation);

% --- Executes on button press in pushbutton_continue.
function pushbutton_continue_Callback(hObject, ~, handles)

% --- Executes on button press in pushbutton_correct.
function pushbutton_correct_Callback(hObject, ~, handles)

global Log

%Recall index
index=get(handles.pushbutton_update,'UserData')+1;

%Place point on image
Log.RelocPoints{index}=impoint(gca,[]);
setColor(Log.RelocPoints{index},'m');

% Construct boundary constraint function and enforce
fcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(Log.RelocPoints{index},fcn);
set(handles.pushbutton_erase,{'Enable','ForegroundColor'},{'on','black'});
set(handles.pushbutton_update,{'Enable','ForegroundColor','UserData'},{'on','black',index});
guidata(hObject,handles);
uiwait(handles.figure_manipulation)

% --- Executes on button press in pushbutton_update.
function pushbutton_update_Callback(hObject, ~, handles)

global Log
Data=get(handles.figure_manipulation,'UserData');
%Find objects that have not been deleted
Log.RelocValid=cell2mat(cellfun(@isvalid,Log.RelocPoints,'UniformOutput',0)');

%Generate list of coordinates for remaining objects
Log.RelocCoor=cell2mat(cellfun(@getPosition,Log.RelocPoints(Log.RelocValid),'UniformOutput',0)');

%Iterate through positions and identify closest feature
MoveFeat=zeros(length(Log.RelocCoor(:,1)),1);
for i = 1:length(Log.RelocCoor(:,1))
    [~,Nearest]=sortrows((Data.ButtonsXCoor-Log.RelocCoor(i,1)).^2+(Data.ButtonsYCoor-Log.RelocCoor(i,2)).^2);
    MoveFeat(i)=Nearest(1);
end

%Update feature locations
Data.ButtonsXCoor(MoveFeat)=Log.RelocCoor(:,1);
Data.ButtonsYCoor(MoveFeat)=Log.RelocCoor(:,2);
Data.AutofindButtons(MoveFeat)=false;

%Reset coordinate tracking variables
iterateDelete=find(Log.RelocValid==1);
for i = 1:length(iterateDelete)
    delete(Log.RelocPoints{iterateDelete(i)});
end
Log.RelocPoints=[];
Log.RelocValid=[];
Log.RelocCoor=[];
set(hObject,'UserData',0);
set(handles.figure_manipulation,'UserData',Data);

%Generate new image with updated coordinates
image=get(handles.uipanel_scroll,'UserData');
AutofindOnImage(image.Surface,Data.AutofindButtons,Data.ButtonsXCoor,Data.ButtonsYCoor,Data.ButtonsRadius,hObject,handles);
API=iptgetapi(Log.ManipulationAPI);
API.replaceImage(get(handles.pushbutton_continue,'UserData'),'PreserveView',1);
guidata(hObject);
uiwait(handles.figure_manipulation);

% --- Executes on button press in pushbutton_erase.
function pushbutton_erase_Callback(hObject, ~, handles)

global Log

erasePoint=impoint(gca,[]);
erasePointPosition=getPosition(erasePoint);

%Find objects that have not been deleted
Log.RelocValid=cell2mat(cellfun(@isvalid,Log.RelocPoints,'UniformOutput',0)');

%Generate list of coordinates for remaining objects
Log.RelocCoor=cell2mat(cellfun(@getPosition,Log.RelocPoints(Log.RelocValid),'UniformOutput',0)');

[~,Nearest]=sortrows((Log.RelocCoor(:,1)-erasePointPosition(1,1)).^2+(Log.RelocCoor(:,2)-erasePointPosition(1,2)).^2);
originalIndex=0;
sumValid=0;
while sumValid~=Nearest(1)
    originalIndex=originalIndex+1;
    sumValid=sumValid+Log.RelocValid(originalIndex);
end
delete(Log.RelocPoints{originalIndex});
delete(erasePoint);

if sum(Log.RelocValid)==0
    set(hObject,{'Enable','ForegroundColor'},{'inactive',[0.8 0.8 0.8]});
end
guidata(hObject,handles)
uiwait(handles.figure_manipulation);

function []=AutofindOnImage(Image,Autofind,XCoor,YCoor,Radius,hObject,handles)
AutoFeatures=find(Autofind==true);
AutoLength=length(AutoFeatures);
MissFeatures=find(Autofind==false);
MissLength=length(MissFeatures);
CellAuto=cell(AutoLength,1);
CellFull=cell(AutoLength+MissLength,1);
iax=cellfun('isempty',CellAuto);
CellFull(iax)={'green'};
imx=cellfun('isempty',CellFull);
CellFull(imx)={'blue'};
DispImage=imadjust(Image,[],[],get(handles.slider_gamma,'Value'));
ImWithAutoMiss=insertShape(DispImage,'circle',[[XCoor(AutoFeatures),YCoor(AutoFeatures),Radius(AutoFeatures)];[XCoor(MissFeatures),YCoor(MissFeatures),Radius(MissFeatures)]],'Color',CellFull,'LineWidth',3);
set(handles.pushbutton_continue,'UserData',ImWithAutoMiss);
guidata(hObject, handles);


% --- Executes on button press in pushbutton_flag.
function pushbutton_flag_Callback(hObject, ~, handles)


% --- Executes during object creation, after setting all properties.
function slider_zoom_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes during object creation, after setting all properties.
function slider_gamma_CreateFcn(hObject, ~, ~)

if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function slider_zoom_Callback(hObject, ~, ~)

global Log

%Grab image controls and update magnification
API=iptgetapi(Log.ManipulationAPI);
API.setMagnification(2^(get(hObject,'Value')-1));

% --- Executes on slider movement.
function slider_gamma_Callback(hObject, ~, handles)

global Log

%Grab image controls, current gamma and update image
Log.TempGamma=get(hObject,'Value');
API=iptgetapi(Log.ManipulationAPI);
display=imadjust(get(handles.pushbutton_continue,'UserData'),[],[],Log.TempGamma);
API.replaceImage(display,'PreserveView',true);

% --- Executes on button press in pushbutton_cancel.
function pushbutton_cancel_Callback(~, ~, handles)
uiresume(handles.figure_manipulation)

% --- Outputs from this function are returned to the command line.
function MITOMIAnalysis_UserEdit_OutputFcn(~, ~, handles)

delete(handles.figure_manipulation)
