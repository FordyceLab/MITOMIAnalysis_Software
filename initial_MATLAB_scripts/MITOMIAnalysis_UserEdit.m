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

% Last Modified by GUIDE v2.5 10-Feb-2017 17:36:52

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
function MITOMIAnalysis_UserEdit_OpeningFcn(hObject, ~, handles, varargin)

global Log

%Set initial control text and pass images and data into objects
Log=varargin{1};
Image=varargin{2};
Data=varargin{3};
set(handles.uipanel_scroll,'UserData',Image);
set(handles.figure_manipulation,'UserData',Data);
set(handles.text_directions,'String',sprintf('1.Navigate to one of the four corners of the array \n\n2. Zoom into the spot to easily see the spot \n\n3. Adjust Gamma to make the background barely visible \n4. Lock the image \n5. Press the ''Set Coordinate'' button \n\n6. Click on three points on the circumference of the spot \n7. Unlock image and continue to the next corner of the array'))
set(handles.pushbutton_update,'UserData',0); %index for button corrections
set(handles.pushbutton_update2,'UserData',0); %index for BG corrections
set(handles.text_directions,'UserData',1); %track correction phase
set(handles.pushbutton_continue,'UserData',0); %curr displayed image holder
guidata(hObject, handles);


%Generate initial image and grab controls
AutofindOnImage(Image.Surface,Data.AutofindButtons,Data.ButtonsXCoor,Data.ButtonsYCoor,Data.ButtonsRadius,Data.Flag,Data.Remove,hObject,handles)
scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
scrollImage = imshow(get(handles.pushbutton_continue,'UserData'),'parent',scrollAxes);
Log.ManipulationAPI = imscrollpanel(handles.uipanel_scroll,scrollImage); 


% UIWAIT makes MITOMIAnalysis_GUI wait for user response (see UIRESUME)
uiwait(handles.figure_manipulation);

% --- Executes on button press in pushbutton_continue.
function pushbutton_continue_Callback(hObject, ~, handles)

global Log
phase=get(handles.text_directions,'UserData');
Image=get(handles.uipanel_scroll,'UserData');
Data=get(handles.figure_manipulation,'UserData');

switch phase
    case 1
        %Clear uncommitted updates
        Log.RelocPoints=[];
        Log.RelocValid=[];
        Log.RelocCoor=[];
        Log.SurfaceFeatureEvaluation='Passed';
        
        %Generate image for next stage
        RemoveFromImage(uint16(squeeze(Image.Captured(:,:,1))),Data.ButtonsXCoor,Data.ButtonsYCoor,Data.ButtonsRadius,Data.Flag,Data.Remove,hObject,handles)
        scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
        scrollImage = imshow(get(handles.pushbutton_continue,'UserData'),'parent',scrollAxes);
        Log.ManipulationAPI = imscrollpanel(handles.uipanel_scroll,scrollImage); 
        
        %Update uibutton controls
        set(handles.text_directions,'UserData',2);
        set(handles.pushbutton_correct,{'Enable','Visible'},{'inactive','off'});
        set(handles.pushbutton_erase,{'Enable','Visible'},{'inactive','off'});
        set(handles.pushbutton_update,{'Enable','Visible'},{'inactive','off'});
        set(handles.pushbutton_remove,{'Enable','Visible'},{'on','on'});
        set(handles.pushbutton_undoremove,{'Enable','Visible','ForegroundColor'},{'on','on','black'});
        set(handles.pushbutton_flag,{'Enable','Visible'},{'on','on'});
        set(handles.pushbutton_undoflag,{'Enable','Visible','ForegroundColor'},{'on','on','black'});
        
    case 2

        Log.CapturedFeatureEvaluation='Passed';
        %Generate image for next stage if user defined
        if sum(sum(Image.Background))~=0
            AutofindOnImage(Image.Background,Data.AutofindChambers,Data.ChamberXCoor,Data.ChamberYCoor,Data.ChamberRadius,Data.Flag,Data.Remove,hObject,handles)
            scrollAxes = axes('parent',handles.uipanel_scroll,'position',[0 0 1 1],'Units','pixels');
            scrollImage = imshow(get(handles.pushbutton_continue,'UserData'),'parent',scrollAxes);
            Log.ManipulationAPI = imscrollpanel(handles.uipanel_scroll,scrollImage); 


            %Update uibutton controls
            set(handles.text_directions,'UserData',3);
            set(handles.pushbutton_remove,{'Enable','Visible'},{'inactive','off'});
            set(handles.pushbutton_undoremove,{'Enable','Visible'},{'inactive','off'});
            set(handles.pushbutton_flag,{'Enable','Visible'},{'inactive','off'});
            set(handles.pushbutton_undoflag,{'Enable','Visible'},{'inactive','off'});
            set(handles.pushbutton_correct2,{'Enable','Visible'},{'on','on'});
            set(handles.pushbutton_erase2,{'Enable','Visible'},{'on','on'});
            set(handles.pushbutton_flag,{'Enable','Visible'},{'inactive','off'});
            set(handles.pushbutton_update2,{'Enable','Visible'},{'on','on'});
        else
            Log.RelocPoints=[];
            Log.RelocValid=[];
            Log.RelocCoor=[];
            Log.ManipulationAPI='Deleted';
            [Data]=CompileData(Image,Data);
            fprintMITOMI(Image,Data);
            uiresume(handles.figure_manipulation)         
        end
        
    case 3
        Log.RelocPoints=[];
        Log.RelocValid=[];
        Log.RelocCoor=[];
        Log.ManipulationAPI='Deleted';
        Log.BackgroundFeatureEvaluation='Passed';
        [Data]=CompileData(Image,Data);
        fprintMITOMI(Image,Data);
        uiresume(handles.figure_manipulation);
        
end

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

set(handles.pushbutton_erase,{'Enable','ForegroundColor'},{'inactive',[0.8,0.8,0.8]});
set(hObject,{'Enable','ForegroundColor'},{'inactive',[0.8 0.8 0.8]});
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
Data.ButtonsXCoor(MoveFeat)=round(Log.RelocCoor(:,1));
Data.ButtonsYCoor(MoveFeat)=round(Log.RelocCoor(:,2));
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
Image=get(handles.uipanel_scroll,'UserData');
AutofindOnImage(Image.Surface,Data.AutofindButtons,Data.ButtonsXCoor,Data.ButtonsYCoor,Data.ButtonsRadius,Data.Flag,Data.Remove,hObject,handles);
API=iptgetapi(Log.ManipulationAPI);
API.replaceImage(get(handles.pushbutton_continue,'UserData'),'PreserveView',1);

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
uiresume()

% --- Outputs from this function are returned to the command line.
function MITOMIAnalysis_UserEdit_OutputFcn(~, ~, handles)
delete(handles.figure_manipulation)

% --- Executes on button press in pushbutton_correct2.
function pushbutton_correct2_Callback(hObject, ~, handles)

global Log

%Recall index
index=get(handles.pushbutton_update2,'UserData')+1;

%Place point on image
Log.RelocPoints{index}=impoint(gca,[]);
setColor(Log.RelocPoints{index},'m');

% Construct boundary constraint function and enforce
fcn = makeConstrainToRectFcn('impoint',get(gca,'XLim'),get(gca,'YLim'));
setPositionConstraintFcn(Log.RelocPoints{index},fcn);
set(handles.pushbutton_erase2,{'Enable','ForegroundColor'},{'on','black'});
set(handles.pushbutton_update2,{'Enable','ForegroundColor','UserData'},{'on','black',index});

guidata(hObject,handles);
uiwait(handles.figure_manipulation)

% --- Executes on button press in pushbutton_update2.
function pushbutton_update2_Callback(hObject, ~, handles)
        
global Log

set(handles.pushbutton_erase2,{'Enable','ForegroundColor'},{'inactive',[0.8,0.8,0.8]});
set(hObject,{'Enable','ForegroundColor'},{'inactive',[0.8 0.8 0.8]});
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
Data.ChamberXCoor(MoveFeat)=round(Log.RelocCoor(:,1));
Data.ChamberYCoor(MoveFeat)=round(Log.RelocCoor(:,2));
Data.AutofindChambers(MoveFeat)=false;

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
Image=get(handles.uipanel_scroll,'UserData');
AutofindOnImage(Image.Background,Data.AutofindChambers,Data.ChamberXCoor,Data.ChamberYCoor,Data.ChamberRadius,Data.Flag,Data.Remove,hObject,handles);

API=iptgetapi(Log.ManipulationAPI);
API.replaceImage(get(handles.pushbutton_continue,'UserData'),'PreserveView',1);
guidata(hObject);
uiwait(handles.figure_manipulation);

% --- Executes on button press in pushbutton_erase2.
function pushbutton_erase2_Callback(hObject, ~, handles)

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
uiwait(handles.figure_manipulation)

function []=AutofindOnImage(Image,Autofind,XCoor,YCoor,Radius,Flag,Remove,hObject,handles)

%Check status of each feature
AutoFeatures=find(Autofind==true & Remove==false & Flag==false);
AutoLength=length(AutoFeatures);
MissFeatures=find(Autofind==false & Remove==false & Flag==false);
MissLength=length(MissFeatures);
FlagFeatures=find(Flag==true & Remove==false);
FlagLength=length(FlagFeatures);

%Setup cell array for coloring features by status
ColorAuto=cell(0,1);
ColorMiss=cell(0,1);
ColorFlag=cell(0,1);

if ~isempty(AutoLength)
    ColorAuto(1:AutoLength)={'green'};
end
if ~isempty(MissLength)
    ColorMiss(1:MissLength)={'blue'};
end
if ~isempty(FlagLength)
    ColorFlag(1:FlagLength)={'magenta'};
end

ColorFeatures=[ColorAuto;ColorMiss;ColorFlag];

DispImage=imadjust(Image,[],[],get(handles.slider_gamma,'Value'));
ImWithAutoMiss=insertShape(DispImage,'circle',[[XCoor(AutoFeatures),YCoor(AutoFeatures),Radius(AutoFeatures)];[XCoor(MissFeatures),YCoor(MissFeatures),Radius(MissFeatures)];[XCoor(FlagFeatures),YCoor(FlagFeatures),Radius(FlagFeatures)]],'Color',ColorFeatures,'LineWidth',3);
set(handles.pushbutton_continue,'UserData',ImWithAutoMiss);
guidata(hObject, handles);

function []=RemoveFromImage(Image,XCoor,YCoor,Radius,Flag,Remove,hObject,handles)
Features=find(Remove==false & Flag==false);
FeaturesLength=length(Features);
FlagFeatures=find(Flag==true & Remove==false);
FlagLength=length(FlagFeatures);

ColorFeat=cell(0,1);
ColorFlag=cell(0,1);

if ~isempty(FeaturesLength)
    ColorFeat(1:FeaturesLength)={'red'};
end
if ~isempty(FlagLength)
    ColorFlag(1:FlagLength)={'magenta'};
end

ColorFeatures=[ColorFeat;ColorFlag];

DispImage=imadjust(Image,[],[],get(handles.slider_gamma,'Value'));
ImWithFlag=insertShape(DispImage,'circle',[[XCoor(Features),YCoor(Features),Radius(Features)];[XCoor(FlagFeatures),YCoor(FlagFeatures),Radius(FlagFeatures)]],'Color',ColorFeatures,'LineWidth',3);
set(handles.pushbutton_continue,'UserData',ImWithFlag);
guidata(hObject,handles);

% --- Executes on button press in pushbutton_remove.
function pushbutton_remove_Callback(hObject, ~, handles)

global Log

Data=get(handles.figure_manipulation,'UserData');
[p1a,p1b]=ginputax(handles.uipanel_scroll,1);
rbbox; 
axes_handle = gca;
p2=get(axes_handle,'CurrentPoint');

%Sort vertices in preparation of data selection
if ( p1a < p2(1,1) )
   lowX = p1a; highX = p2(1,1);
else
   lowX = p2(1,1); highX = p1a;
end

if ( p1b < p2(1,2) )
   lowY = p1b; highY = p2(1,2);
else
   lowY = p2(1,2); highY = p1b;
end

RemovedFeatures=((Data.ButtonsXCoor > lowX) & (Data.ButtonsXCoor < highX) & (Data.ButtonsYCoor > lowY) & (Data.ButtonsYCoor < highY));
Data.Remove(RemovedFeatures)=true;
set(handles.figure_manipulation,'UserData',Data);
guidata(hObject,Data);

API=iptgetapi(Log.ManipulationAPI);
Image=get(handles.uipanel_scroll,'UserData');
RemoveFromImage(uint16(squeeze(Image.Captured(:,:,1))),Data.ButtonsXCoor,Data.ButtonsYCoor,Data.ButtonsRadius,Data.Flag,Data.Remove,hObject,handles)

API=iptgetapi(Log.ManipulationAPI);
API.replaceImage(get(handles.pushbutton_continue,'UserData'),'PreserveView',1);

% --- Executes on button press in pushbutton_undoremove.
function pushbutton_undoremove_Callback(hObject, ~, handles)

global Log

Data=get(handles.figure_manipulation,'UserData');
[p1a,p1b]=ginputax(handles.uipanel_scroll,1);
rbbox; 
axes_handle = gca;
p2=get(axes_handle,'CurrentPoint');

%Sort vertices in preparation of data selection
if ( p1a < p2(1,1) )
   lowX = p1a; highX = p2(1,1);
else
   lowX = p2(1,1); highX = p1a;
end

if ( p1b < p2(1,2) )
   lowY = p1b; highY = p2(1,2);
else
   lowY = p2(1,2); highY = p1b;
end

RestoredFeatures=((Data.ButtonsXCoor > lowX) & (Data.ButtonsXCoor < highX) & (Data.ButtonsYCoor > lowY) & (Data.ButtonsYCoor < highY));
Data.Remove(RestoredFeatures)=false;
set(handles.figure_manipulation,'UserData',Data);
guidata(hObject,Data);

Image=get(handles.uipanel_scroll,'UserData');
RemoveFromImage(uint16(squeeze(Image.Captured(:,:,1))),Data.ButtonsXCoor,Data.ButtonsYCoor,Data.ButtonsRadius,Data.Flag,Data.Remove,hObject,handles)

API=iptgetapi(Log.ManipulationAPI);
API.replaceImage(get(handles.pushbutton_continue,'UserData'),'PreserveView',1);

% --- Executes on button press in pushbutton_flag.
function pushbutton_flag_Callback(hObject, ~, handles)

global Log

Data=get(handles.figure_manipulation,'UserData');
[p1a,p1b]=ginputax(handles.uipanel_scroll,1);
rbbox; 
axes_handle = gca;
p2=get(axes_handle,'CurrentPoint');

%Sort vertices in preparation of data selection
if ( p1a < p2(1,1) )
   lowX = p1a; highX = p2(1,1);
else
   lowX = p2(1,1); highX = p1a;
end

if ( p1b < p2(1,2) )
   lowY = p1b; highY = p2(1,2);
else
   lowY = p2(1,2); highY = p1b;
end

FlaggedFeatures=((Data.ButtonsXCoor > lowX) & (Data.ButtonsXCoor < highX) & (Data.ButtonsYCoor > lowY) & (Data.ButtonsYCoor < highY));
Data.Flag(FlaggedFeatures)=true;
set(handles.figure_manipulation,'UserData',Data);
guidata(hObject,Data);

Image=get(handles.uipanel_scroll,'UserData');
RemoveFromImage(uint16(squeeze(Image.Captured(:,:,1))),Data.ButtonsXCoor,Data.ButtonsYCoor,Data.ButtonsRadius,Data.Flag,Data.Remove,hObject,handles)
 
API=iptgetapi(Log.ManipulationAPI);
API.replaceImage(get(handles.pushbutton_continue,'UserData'),'PreserveView',1);

% --- Executes on button press in pushbutton_undoflag.
function pushbutton_undoflag_Callback(hObject, ~, handles)

global Log

Data=get(handles.figure_manipulation,'UserData');
[p1a,p1b]=ginputax(handles.uipanel_scroll,1);
rbbox; 
axes_handle = gca;
p2=get(axes_handle,'CurrentPoint');

%Sort vertices in preparation of data selection
if ( p1a < p2(1,1) )
   lowX = p1a; highX = p2(1,1);
else
   lowX = p2(1,1); highX = p1a;
end

if ( p1b < p2(1,2) )
   lowY = p1b; highY = p2(1,2);
else
   lowY = p2(1,2); highY = p1b;
end

UnflaggedFeatures=((Data.ButtonsXCoor > lowX) & (Data.ButtonsXCoor < highX) & (Data.ButtonsYCoor > lowY) & (Data.ButtonsYCoor < highY));
Data.Flag(UnflaggedFeatures)=false;
set(handles.figure_manipulation,'UserData',Data);
guidata(hObject,Data);

Image=get(handles.uipanel_scroll,'UserData');
RemoveFromImage(uint16(squeeze(Image.Captured(:,:,1))),Data.ButtonsXCoor,Data.ButtonsYCoor,Data.ButtonsRadius,Data.Flag,Data.Remove,hObject,handles)
 
API=iptgetapi(Log.ManipulationAPI);
API.replaceImage(get(handles.pushbutton_continue,'UserData'),'PreserveView',1);

function [Data]=CompileData(Image,Data)

global Log

    index=0;
    Log.NumWells=Log.Rows*Log.Cols;
    Log.NumSamples=sum(~Data.Remove);
    WAIT=waitbar(0,'Extracting data from features...','Name','Data Extraction Percent Complete: ');

    %Masks are made such that button is always centered
    %Chamber is then inserted into mask relative to button coordinates
    
    if sum(sum(Image.Background))~=0 %If a background image was used
                             
        window=4*Log.ApproxBackgroundRadius+1;
        [MaskX,MaskY]=meshgrid(1:window,1:window);
        
        for i=1:Log.NumWells
        
            index=index+double(~Data.Remove(i));
            
            if Data.Remove(i)==0
               Data.Index(i)=index;
            end
            
            %Generate unique button masks in case variable radius is used
            ButtonMask=uint16(sqrt((MaskX-(2*Log.ApproxBackgroundRadius+1)).^2+(MaskY-(2*Log.ApproxBackgroundRadius+1)).^2)<=Data.ButtonsRadius(i));
            ChamberBGMask =       uint16(sqrt((MaskX-(Data.ChamberXCoor(i)-Data.ButtonsXCoor(i)+2*Log.ApproxBackgroundRadius+1)).^2+(MaskY-(Data.ChamberYCoor(i)-Data.ButtonsYCoor(i)+2*Log.ApproxBackgroundRadius+1)).^2)>=Log.ApproxBackgroundRadius*1.1 & sqrt((MaskX-(Data.ChamberXCoor(i)-Data.ButtonsXCoor(i)+2*Log.ApproxBackgroundRadius+1)).^2+(MaskY-(Data.ChamberYCoor(i)-Data.ButtonsYCoor(i)+2*Log.ApproxBackgroundRadius+1)).^2)<=Log.ApproxBackgroundRadius*1.3 );
            ChamberNoButtonMask = uint16(sqrt((MaskX-(Data.ChamberXCoor(i)-Data.ButtonsXCoor(i)+2*Log.ApproxBackgroundRadius+1)).^2+(MaskY-(Data.ChamberYCoor(i)-Data.ButtonsYCoor(i)+2*Log.ApproxBackgroundRadius+1)).^2)<=Log.ApproxBackgroundRadius &~ sqrt((MaskX-(2*Log.ApproxBackgroundRadius+1)).^2+(MaskY-(2*Log.ApproxBackgroundRadius+1)).^2)<=Log.ApproxBackgroundRadius*1.2 );
            Data.ChamberAreaFG(i)=sum(sum(ChamberNoButtonMask));
            Data.ChamberAreaBG(i)=sum(sum(ChamberBGMask));    
            Data.ButtonsAreaFG(i)=sum(sum(ButtonMask));
            Data.ButtonsAreaBG(i)=sum(sum(ChamberNoButtonMask));
            
            %Collect data from solubilized/background chambers
            ImageSol=Image.Background((Data.ButtonsYCoor(i)-2*Log.ApproxBackgroundRadius):(Data.ButtonsYCoor(i)+2*Log.ApproxBackgroundRadius),(Data.ButtonsXCoor(i)-2*Log.ApproxBackgroundRadius):(Data.ButtonsXCoor(i)+2*Log.ApproxBackgroundRadius));
            DNAChamber = double(ImageSol.*ChamberNoButtonMask);
            DNAChamberBG = double(ImageSol.*ChamberBGMask);

            Data.SolubilizedMedianFG(i)=median(DNAChamber(DNAChamber(:)>0));
            Data.SolubilizedMeanFG(i)=mean(DNAChamber(DNAChamber(:)>0));
            Data.SolubilizedSTDFG(i)=std(DNAChamber(DNAChamber(:)>0));
            Data.SolubilizedTotalFG(i)=sum(DNAChamber(DNAChamber(:)>0));
            Data.SolubilizedMedianBG(i)=median(DNAChamberBG(DNAChamberBG(:)>0));
            Data.SolubilizedMeanBG(i)=mean(DNAChamberBG(DNAChamberBG(:)>0));
            Data.SolubilizedSTDBG(i)=std(DNAChamberBG(DNAChamberBG(:)>0));
            Data.SolubilizedTotalBG(i)=sum(DNAChamberBG(DNAChamberBG(:)>0))*Data.ChamberAreaFG(i)./Data.ChamberAreaBG(i);
            Data.SolubilizedFractionSaturatedFG(i)=length(find(DNAChamber==65535))./length(find(DNAChamber(:)>0));
            Data.SolubilizedFractionSaturatedBG(i)=length(find(DNAChamberBG==65535))./length(find(DNAChamberBG(:)>0));           

            %Collect data from surface immobilized molecule
            ImageSur=Image.Surface((Data.ButtonsYCoor(i)-2*Log.ApproxBackgroundRadius):(Data.ButtonsYCoor(i)+2*Log.ApproxBackgroundRadius),(Data.ButtonsXCoor(i)-2*Log.ApproxBackgroundRadius):(Data.ButtonsXCoor(i)+2*Log.ApproxBackgroundRadius));        
            SurfaceButton=double(ImageSur.*ButtonMask);
            SurfaceBG=double(ImageSur.*ChamberNoButtonMask);

            Data.SurfaceMedianFG(i)=median(SurfaceButton(SurfaceButton(:)>0));
            Data.SurfaceAverageFG(i)=mean(SurfaceButton(SurfaceButton(:)>0));
            Data.SurfaceSTDFG(i)=std(SurfaceButton(SurfaceButton(:)>0));
            Data.SurfaceTotalFG(i)=sum(SurfaceButton(SurfaceButton(:)>0));
            Data.SurfaceMedianBG(i)=median(SurfaceBG(SurfaceBG(:)>0));
            Data.SurfaceAverageBG(i)=mean(SurfaceBG(SurfaceBG(:)>0));
            Data.SurfaceSTDBG(i)=std(SurfaceBG(SurfaceBG(:)>0));
            Data.SurfaceTotalBG(i)=sum(SurfaceBG(SurfaceBG(:)>0))*Data.ButtonsAreaFG(i)./Data.ButtonsAreaBG(i);
            Data.SurfaceFractionSaturatedFG(i)=length(find(SurfaceButton==65535))./length(find(SurfaceButton(:)>0));
            Data.SurfaceFractionSaturatedBG(i)=length(find(SurfaceBG==65535))./length(find(SurfaceBG(:)>0));

            %Collect data from captured molecule images
            for FrameCap=1:Log.CapturedFrames

                ImageCap=uint16(squeeze(Image.Captured((Data.ButtonsYCoor(i)-2*Log.ApproxBackgroundRadius):(Data.ButtonsYCoor(i)+2*Log.ApproxBackgroundRadius),(Data.ButtonsXCoor(i)-2*Log.ApproxBackgroundRadius):(Data.ButtonsXCoor(i)+2*Log.ApproxBackgroundRadius),FrameCap)));
                CapturedButton=double(ImageCap.*ButtonMask);
                CapturedBG=double(ImageCap.*ChamberNoButtonMask);

                Data.CapturedMedianFG(i,FrameCap)=median(CapturedButton(CapturedButton(:)>0));
                Data.CapturedAverageFG(i,FrameCap)=mean(CapturedButton(CapturedButton(:)>0));
                Data.CapturedSTDFG(i,FrameCap)=std(CapturedButton(CapturedButton(:)>0));
                Data.CapturedTotalFG(i,FrameCap)=sum(CapturedButton(CapturedButton(:)>0));
                Data.CapturedMedianBG(i,FrameCap)=median(CapturedBG(CapturedBG(:)>0));
                Data.CapturedAverageBG(i,FrameCap)=mean(CapturedBG(CapturedBG(:)>0));
                Data.CapturedSTDBG(i,FrameCap)=std(CapturedBG(CapturedBG(:)>0));
                Data.CapturedTotalBG(i,FrameCap)=sum(CapturedBG(CapturedBG(:)>0))*Data.ButtonsAreaFG(i)./Data.ButtonsAreaBG(i); 
                Data.CapturedFractionSaturatedFG(i,FrameCap)=length(find(CapturedButton==65535))./length(find(CapturedButton(:)>0));
                Data.CapturedFractionSaturatedBG(i,FrameCap)=length(find(CapturedBG==65535))./length(find(CapturedBG(:)>0));

            end

            waitbar(i/Log.NumWells,WAIT,sprintf('%6.3f',i/Log.NumWells*100));

        end

        delete(WAIT)
        
    else %a background image was not used
        
        window=4*Log.ApproxButtonRadius+1;
        [MaskX,MaskY]=meshgrid(1:window,1:window);
        
        for j=1:Log.NumWells
            
            index=index+double(~Data.Remove(j));
            
            if Data.Remove(j)==0
               Data.Index(j)=index;
            end

            %Generate unique button masks in case variable radius is used          
            ButtonMask=uint16(sqrt((MaskX-(2*Log.ApproxButtonRadius+1)).^2+(MaskY-(2*Log.ApproxButtonRadius+1)).^2)<=Data.ButtonsRadius(j));
            ChamberNoButtonMask = uint16(sqrt((MaskX-(2*Log.ApproxButtonRadius+1)).^2+(MaskY-(2*Log.ApproxButtonRadius+1)).^2)>=Data.ButtonsRadius(j)*1.3 & sqrt((MaskX-(2*Log.ApproxButtonRadius+1)).^2+(MaskY-(2*Log.ApproxButtonRadius+1)).^2)<=Data.ButtonsRadius(j)*1.9 );            
            Data.ButtonsAreaFG(j)=sum(sum(ButtonMask));
            Data.ButtonsAreaBG(j)=sum(sum(ChamberNoButtonMask));

            %Collect data from surface immobilized molecules            
            ImageSur=Image.Surface((Data.ButtonsYCoor(j)-2*Log.ApproxButtonRadius):(Data.ButtonsYCoor(j)+2*Log.ApproxButtonRadius),(Data.ButtonsXCoor(j)-2*Log.ApproxButtonRadius):(Data.ButtonsXCoor(j)+2*Log.ApproxButtonRadius));
            SurfaceButton=double(ImageSur.*ButtonMask);
            SurfaceBG=double(ImageSur.*ChamberNoButtonMask);

            Data.SurfaceMedianFG(j)=median(SurfaceButton(SurfaceButton(:)>0));
            Data.SurfaceAverageFG(j)=mean(SurfaceButton(SurfaceButton(:)>0));
            Data.SurfaceSTDFG(j)=std(SurfaceButton(SurfaceButton(:)>0));
            Data.SurfaceTotalFG(j)=sum(SurfaceButton(SurfaceButton(:)>0));
            Data.SurfaceMedianBG(j)=median(SurfaceBG(SurfaceBG(:)>0));
            Data.SurfaceAverageBG(j)=mean(SurfaceBG(SurfaceBG(:)>0));
            Data.SurfaceSTDBG(j)=std(SurfaceBG(SurfaceBG(:)>0));
            Data.SurfaceTotalBG(j)=sum(SurfaceBG(SurfaceBG(:)>0))*Data.ButtonsAreaFG(j)./Data.ButtonsAreaBG(j);
            Data.SurfaceFractionSaturatedFG(j)=length(find(SurfaceButton==65535))./length(find(SurfaceButton(:)>0));
            Data.SurfaceFractionSaturatedBG(j)=length(find(SurfaceBG==65535))./length(find(SurfaceBG(:)>0));

            %Collect data from captured molecule images
            for FrameCap=1:Log.CapturedFrames

                ImageCap=uint16(squeeze(Image.Captured((Data.ButtonsYCoor(j)-2*Log.ApproxButtonRadius):(Data.ButtonsYCoor(j)+2*Log.ApproxButtonRadius),(Data.ButtonsXCoor(j)-2*Log.ApproxButtonRadius):(Data.ButtonsXCoor(j)+2*Log.ApproxButtonRadius),FrameCap)));
                CapturedButton=double(ImageCap.*ButtonMask);
                CapturedBG=double(ImageCap.*ChamberNoButtonMask);

                Data.CapturedMedianFG(j,FrameCap)=median(CapturedButton(CapturedButton(:)>0));
                Data.CapturedAverageFG(j,FrameCap)=mean(CapturedButton(CapturedButton(:)>0));
                Data.CapturedSTDFG(j,FrameCap)=std(CapturedButton(CapturedButton(:)>0));
                Data.CapturedTotalFG(j,FrameCap)=sum(CapturedButton(CapturedButton(:)>0));
                Data.CapturedMedianBG(j,FrameCap)=median(CapturedBG(CapturedBG(:)>0));
                Data.CapturedAverageBG(j,FrameCap)=mean(CapturedBG(CapturedBG(:)>0));
                Data.CapturedSTDBG(j,FrameCap)=std(CapturedBG(CapturedBG(:)>0));
                Data.CapturedTotalBG(j,FrameCap)=sum(CapturedBG(CapturedBG(:)>0))*Data.ButtonsAreaFG(j)./Data.ButtonsAreaBG(j); 
                Data.CapturedFractionSaturatedFG(j,FrameCap)=length(find(CapturedButton==65535))./length(find(CapturedButton(:)>0));
                Data.CapturedFractionSaturatedBG(j,FrameCap)=length(find(CapturedBG==65535))./length(find(CapturedBG(:)>0));

            end

            waitbar(j/Log.NumWells,WAIT,sprintf('%6.3f',j/Log.NumWells*100));

        end

        delete(WAIT)

    end         

function []=fprintMITOMI(Image,Data)

global Log

if isempty(Log.InitiatedFileSave)
    savename=['Completed_' Log.NameOutput];
    savemat=[savename '.mat'];
else
    savename=['Reviewed_' Log.NameOutput '_on_' datestr(now,30)];
    savemat=[ savename '.mat'];
end
    Log.InitiatedFileSave='Passed';
    saveMSGBOX=msgbox('Saving large files. Please be patient','Saving files');
    save(savemat,'Log','Image','Data','-v7.3')

    %create string header for dissociation data
    HeaderFormat={'Index','ColIndex','RowIndex','Removed','Flagged','ButtonXCoor','ButtonYCoor','ButtonRadius','ButtonAreaFG','ButtonAreaBG','ButtonAutoFind','BNDMedFG','BNDAvgFG','BNDStdFG','BNDSumFG','BNDSatFG','BNDMedBG','BNDAvgBG','BNDStdBG','BNDSumBG','BNDSatBG'};

    for z=1:Log.CapturedFrames
        HeaderFormat(21+z)={['CAPMedFG' num2str(z)]};
        HeaderFormat(21+z+Log.CapturedFrames)={['CAPAvgFG' num2str(z)]};
        HeaderFormat(21+z+Log.CapturedFrames*2)={['CAPStdFG' num2str(z)]};
        HeaderFormat(21+z+Log.CapturedFrames*3)={['CAPSumFG' num2str(z)]};
        HeaderFormat(21+z+Log.CapturedFrames*4)={['CAPSatFG' num2str(z)]};
        HeaderFormat(21+z+Log.CapturedFrames*5)={['CAPMedBG' num2str(z)]};
        HeaderFormat(21+z+Log.CapturedFrames*6)={['CAPAvgBG' num2str(z)]};
        HeaderFormat(21+z+Log.CapturedFrames*7)={['CAPStdBG' num2str(z)]};
        HeaderFormat(21+z+Log.CapturedFrames*8)={['CAPSumBG' num2str(z)]};
        HeaderFormat(21+z+Log.CapturedFrames*9)={['CAPSatBG' num2str(z)]};
    end

    HeaderFormat(end+1:end+16)={'ChamberXCoor','ChamberYCoor','ChamberRadius','ChamberAreaFG','ChamberAreaBG','ChamberAutoFind','SOLMedFG','SOLAvgFG','SOLStdFG','SOLSumFG','SOLSatFG','SOLMedBG','SOLAvgBG','SOLStdBG','SOLSumBG','SOLSatBG'};

    %Save txt file with headers
    DataText=fopen([savename '.txt'],'w');
    IntermediateFormat=struct2cell(Data);
    DataFormat=cat(2,IntermediateFormat{:});
    fprintf(DataText,'%s\t',HeaderFormat{:} );
    fprintf(DataText,'\r\n');

    for Z=1:Log.NumWells;
        fprintf(DataText,'%f\t',DataFormat(Z,:)');
        fprintf(DataText,'\r\n');
    end
    
    fclose(DataText);
    
    %Save txt file without headers
    savetxt=[savename '_NoHeaders.txt'];
    save(savetxt,'DataFormat','-ascii','-tabs')
    
    if ishandle(saveMSGBOX)
        close(saveMSGBOX);
    end
    

