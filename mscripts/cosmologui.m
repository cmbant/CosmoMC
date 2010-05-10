function varargout = cosmologui(varargin)
% COSMOLOGUI M-file for cosmologui.fig
%      COSMOLOGUI, by itself, creates a new COSMOLOGUI or raises the existing
%      singleton*.
%
%      H = cosmologui returns the handle to a new cosmologui or the handle to
%      the existing singleton*.
%
%      cosmologui('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in cosmologui.M with the given input arguments.
%
%      cosmologui('Property','Value',...) creates a new cosmologui or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cosmologui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cosmologui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cosmologui

% Last Modified by GUIDE v2.5 16-Oct-2004 10:45:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cosmologui_OpeningFcn, ...
                   'gui_OutputFcn',  @cosmologui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

% doesn't work if docked
set(0,'DefaultFigureWindowStyle','normal') 


if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% -----------------------------------------------------------------
% -----------------------------------------------------------------
% --- Executes just before cosmologui is made visible.
function cosmologui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cosmologui (see VARARGIN)

% Choose default command line output for cosmologui
handles.output = hObject;

movegui([10 -10]);

% choices associated with making the plot
handles.xcol=1; % Line number
handles.ycol=1; % Probability
handles.zcol=1; % Blank
handles.linewidth = 2; % default line width
handles.linecolor = [1 0 0]; % default line color (red)
handles.linestyle = '-'; % default line style (solid)
handles.clevels = [0.68 0.95]; % default confidence levels
handles.res = 0.1; % default smoothing width, as fraction of plot width
handles.nsampperbox = 5; % default no. samples per bin
handles.plottype = 1; % default: do contours and shading
handles.fontsize = 14;
handles.markersize=2; % for scatter plot
handles.colormap=1; % jet
handles.isheld = 0; % don't hold plots by default
handles.cuts=[]; % no cuts at first. Don't save these either.

% load in choices from last time
fid=fopen('cosmologui_prefs.mat');
if (fid>0)
    fclose(fid);
    htmp=load('cosmologui_prefs');
    % how to do this better!?
    handles.xcol=htmp.xcol;
    handles.ycol=htmp.ycol;
    handles.zcol=htmp.zcol;
    handles.linewidth=htmp.linewidth;
    handles.linecolor=htmp.linecolor;
    handles.linestyle=htmp.linestyle;
    handles.clevels=htmp.clevels;
    handles.res=htmp.res; %resolution
    handles.nsampperbox=htmp.nsampperbox;
    handles.plottype=htmp.plottype;
    handles.fontsize=htmp.fontsize;
    handles.markersize=htmp.markersize;
    handles.colormap=htmp.colormap;
    handles.isheld=htmp.isheld;
end

% set the menus to reflect the choices
set(handles.PlotTypePopup,'Value',handles.plottype);
set(handles.ResEdit,'Value',handles.res); % resolution
set(handles.HoldRadio,'Value',handles.isheld); % whether hold is on or not

if ~isempty(varargin)
    % samples may be given as arguments, so no need to read them in
    handles.samples=varargin{1}; % messy.. copying around big arrays?
    handles.names=varargin{2};
    handles.isampled=varargin{3};
    handles.cuts=varargin{4};
    handles.figno=varargin{5};

    handles=doxyzpulldown(handles);
    
    % open new figure, if necessary
    if (~isfield(handles,'figno'))
        scrsz = get(0,'ScreenSize');
        handles.figno=figure('Position',[10 10 scrsz(3)/2 scrsz(4)/2]);
    end
    plotsamps(handles);
    
else
    
    % load in some samples using the getsamps gui
     handles=Load_Callback(hObject, eventdata, handles);
% !!!!?? what if it comes back empty?     
end

set(handles.ColormapPopup,'Value',handles.colormap);
ColormapPopup_Callback(hObject, eventdata, handles);
set(handles.ResEdit,'String',handles.res); % resolution

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cosmologui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = cosmologui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
%varargout{1} = handles.output;
if (isempty(handles))
    varargout=[];
else
  varargout = struct2cell(handles);
end


% --------------------------------------------------------------------
function handles=Load_Callback(hObject, eventdata, handles)
% hObject    handle to Open (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load in some samples using the getsamps gui
[handles.samples,handles.names,handles.isampled,handles.cuts,h,handles.colors]=...
    getsamps('Position',[20 20 150 22]);
if ((isempty(handles.samples))|(handles.samples(1)<0))
    % go to quit button
    Exit_Callback;
end

% shut down the getsamps gui
if (ishandle(h))
    status=close(h);
    if (status==0) fprintf(1,'Error closing figure\n'); end
else
    fprintf(1,'Error loading samples\n');
    return;
end
if (status==0) fprintf(1,'Error closing figure\n'); end

handles=doxyzpulldown(handles);

% open new figure, if necessary
if (~isfield(handles,'figno'))
    scrsz = get(0,'ScreenSize');
    handles.figno=figure('Position',[10 10 scrsz(3)/2 scrsz(4)/2]);
    handles.logoaxes=axes;
    h=text(0.5,1.05,'Created using CosmoloGUI','color',[0 0.5 0.1]);
    set(handles.logoaxes,'visible','off','layer','top')
    v=axes;
    
end

plotsamps(handles);

guidata(hObject,handles); % store the changes


% --- Executes during object creation, after setting all properties.
function XlabelPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to XlabelPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in XlabelPopupmenu.
function XlabelPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to XlabelPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns XlabelPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from XlabelPopupmenu

handles.xcol = handles.xvals(get(handles.XlabelPopupmenu, 'Value'));
figure(handles.figno); clflogo % no point in overplotting if y axis is now different
plotsamps(handles);

% change the pull-down menus so eg. can't plot x against x
handles=doxyzpulldown(handles);
guidata(hObject,handles); % store the changes

% --- Executes during object creation, after setting all properties.
function YlabelPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to YlabelPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in YlabelPopupmenu.
function YlabelPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to YlabelPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns YlabelPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from YlabelPopupmenu

handles.ycol = handles.yvals(get(handles.YlabelPopupmenu, 'Value'));
figure(handles.figno); clflogo % no point in overplotting if y axis is now different
plotsamps(handles);

% change the pull-down menus so eg. can't plot x against x
handles=doxyzpulldown(handles);
guidata(hObject,handles); % store the changes


% --- Executes during object creation, after setting all properties.
function ZlabelPopupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZlabelPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in ZlabelPopupmenu.
function ZlabelPopupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to ZlabelPopupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ZlabelPopupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ZlabelPopupmenu

handles.zcol = handles.zvals(get(handles.ZlabelPopupmenu, 'Value'));

%if ((handles.zcol>1)&(handles.xcol==1)) % this is a bit naff..
%        handles.xcol=handles.isampled(2); 
%        fprintf(1,'*** NB. cant have X=Line number and Z not blank ***\n')
%        fprintf(1,'*** Switched X to next available column ***\n')
%end
%if ((handles.zcol>1)&(handles.ycol==1)) % this is a bit naff..
%        handles.ycol=handles.isampled(3); 
%        fprintf(1,'*** NB. cant have Y=Line number and Z not blank ***\n')
%        fprintf(1,'*** Switched Y to next available column ***\n')
%end

figure(handles.figno); clf % won't want to overplot if x axis has changed
plotsamps(handles);

% change the pull-down menus so eg. can't plot x against x
handles=doxyzpulldown(handles);

guidata(hObject,handles); % store the changes

% --- Executes during object creation, after setting all properties.
function PlotTypePopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PlotTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in PlotTypePopup.
function PlotTypePopup_Callback(hObject, eventdata, handles)
% hObject    handle to PlotTypePopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns PlotTypePopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PlotTypePopup

handles.plottype=get(hObject,'Value');
figure(handles.figno);
plotsamps(handles);
guidata(hObject,handles); % store the changes


% --- Executes during object creation, after setting all properties.
function ColormapPopup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColormapPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in ColormapPopup.
function ColormapPopup_Callback(hObject, eventdata, handles)
% hObject    handle to ColormapPopup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ColormapPopup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ColormapPopup

handles.colormap = get(handles.ColormapPopup, 'Value');
guidata(hObject,handles); % store the changes

figure(handles.figno);
switch handles.colormap
    case 1
        colormap(jet)
    case 2
        red
    case 3
        blue
    case 4
        green
    case 5
        bw
    case 6
        colormap(hsv)
    case 7
        yellowblue
end

% --- Executes on button press in InvertButton.
function InvertButton_Callback(hObject, eventdata, handles)
% hObject    handle to InvertButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(handles.figno);
invert

% --- Executes on button press in LightenButton.
function LightenButton_Callback(hObject, eventdata, handles)
% hObject    handle to LightenButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(handles.figno);
lighten


% --- Executes on button press in ClearButton.
function ClearButton_Callback(hObject, eventdata, handles)
% hObject    handle to ClearButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(handles.figno); 
clflogo

% --- Executes during object creation, after setting all properties.
function ResEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ResEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

%set(handles.ResEdit,'Value',handles.res); % resolution


% --- Executes on changing Text in Resolution field
function ResEdit_Callback(hObject, eventdata, handles)
% hObject    handle to ResEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ResEdit as text
%        str2double(get(hObject,'String')) returns contents of ResEdit as a double

handles.res=str2double(get(hObject,'String'))
guidata(hObject,handles); % store the changes

plotsamps(handles);


% --- Executes on button press in UpdateButton.
function UpdateButton_Callback(hObject, eventdata, handles)
% hObject    handle to UpdateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

plotsamps(handles);

% -------------------------------------------------------------------
% --- Called by functions above, to sort out X,Y,Z pulldown menus
% after reading in new data, or changing the value xcol, ycol, zcol
function handles=doxyzpulldown(handles);

% check that we aren't trying to plot x against x
if (  ((handles.xcol==handles.ycol)&(handles.xcol~=1))    | ...
        ((handles.ycol==handles.zcol)&(handles.ycol~=1))    | ...
        ((handles.zcol==handles.xcol)&(handles.zcol~=1)) )
    % (if they are both 1 then this is a special case)
    handles.xcol=1; % reset everything to 1 
    handles.ycol=1; 
    handles.zcol=1; 
end    
    
% make the pulldown menus match the names
% and remove ycol and zcol from the list (unless 1)
handles.xvals=[];
for i=1:length(handles.isampled)
    is=handles.isampled(i);
    if ( ((is~=handles.ycol)|(is==1)) & ((is~=handles.zcol)|(is==1)) )
        handles.xvals=[handles.xvals is];
    end
end
% sort out a special case: if zcol isn't blank, can't have xcol=line no.
if (handles.zcol>1) 
    handles.xvals=handles.xvals(2:length(handles.xvals)); 
    tmpcell=handles.names(handles.xvals);
else
    tmpcell=handles.names(handles.xvals);
    tmpcell(1)={'Line number'};
end
set(handles.XlabelPopupmenu,'String',tmpcell);

% sim Y
handles.yvals=[];
for i=1:length(handles.isampled)
    is=handles.isampled(i);
    if ( ((is~=handles.xcol)|(is==1)) & ((is~=handles.zcol)|(is==1)) )
        handles.yvals=[handles.yvals is];
    end
end
% sort out a special case: if zcol isn't blank, can't have xcol=line no.
if (handles.zcol>1) 
    handles.yvals=handles.yvals(2:length(handles.yvals)); 
    tmpcell=handles.names(handles.yvals);
else
    tmpcell=handles.names(handles.yvals);
    tmpcell(1)={'Probability'};
end
set(handles.YlabelPopupmenu,'String',tmpcell);

% sim Z
handles.zvals=[];
for i=1:length(handles.isampled)
    is=handles.isampled(i);
    if ( ((is~=handles.xcol)|(is==1)) & ((is~=handles.ycol)|(is==1)) )
        handles.zvals=[handles.zvals is];
    end
end
% insert an extra option: probability up z axis:
if ((handles.xcol==1)|(handles.ycol==1))
  handles.zvals=1;
  tmpcell=handles.names(handles.zvals);
  tmpcell(1)={'- - - - - - - - - - -'};
else
  handles.zvals=[1 2 handles.zvals(2:length(handles.zvals))];
  tmpcell=handles.names(handles.zvals);
  tmpcell(1)={'- - - - - - - - - - -'};
  tmpcell(2)={'Probability'};
end
set(handles.ZlabelPopupmenu,'String',tmpcell);

% check that selected xyzcols were actually sampled from
ix=find(handles.xcol==handles.xvals); % hopefully length(ix)==1 !
iy=find(handles.ycol==handles.yvals);
iz=find(handles.zcol==handles.zvals);

if (isempty(ix)|isempty(iy)|isempty(iz))
    % one of the old columns was not sampled from in the new file
    handles.xcol=1; % reset everything to 1 
    handles.ycol=1; 
    handles.zcol=1; 
    set(handles.XlabelPopupmenu,'Value',1);
    set(handles.YlabelPopupmenu,'Value',1);
    set(handles.ZlabelPopupmenu,'Value',1);
else 
    % set up the popupmenus
    set(handles.XlabelPopupmenu,'Value',ix(1));
    set(handles.YlabelPopupmenu,'Value',iy(1));
    set(handles.ZlabelPopupmenu,'Value',iz(1));
end

% if the z axis is not blank, change the plottype menu
if (handles.zcol==1) % blank
    tmpcell={'Contours and Shading'; 'Contours'; 'Shading'; ...
            'Scatter (weights ignored)'; 'Histogram'};
    set(handles.PlotTypePopup,'String',tmpcell);
    if (handles.plottype>length(tmpcell))
        handles.plottype=1; % reset to default
    end
    set(handles.PlotTypePopup,'Value',handles.plottype);
end
if (handles.zcol==2) % probability
    % can't really do just contours in 3d, but otherwise same
    % as handles.zcol==1;
    tmpcell={'Scatter (weights ignored)';'Surface';'Histogram'};
    set(handles.PlotTypePopup,'String',tmpcell);
    if (handles.plottype>length(tmpcell))
        handles.plottype=1; % reset to default
    end
    set(handles.PlotTypePopup,'Value',handles.plottype);
end
if (handles.zcol>2) % probability
    % can't really do just contours in 3d, but otherwise same
    % as handles.zcol==1;
    tmpcell={'Scatter (weights ignored)'; '68 % isosurface'; '95 % isosurface'};
    set(handles.PlotTypePopup,'String',tmpcell);
    if (handles.plottype>length(tmpcell))
        handles.plottype=1; % reset to default
    end
    set(handles.PlotTypePopup,'Value',handles.plottype);
end

% --------------------------------------------------------------------
function Exit_Callback(hObject, eventdata, handles)
% hObject    handle to Exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%pos_size = get(handles.figure1,'Position')
%user_response = modaldlg('Title','Confirm Close');
% couldn't get the above modaldlg to work yet.. see matlab help for 
% example of how it should work..

% save preferences for next time
% first remove the large data entry from handles, to save space
%rmfield(handles,'samples')
%save cosmologui_prefs handles
if (isfield(handles,'xcol'))
  xcol=handles.xcol;
  ycol=handles.ycol;
  zcol=handles.zcol;
  linewidth=handles.linewidth;
  linecolor=handles.linecolor;
  linestyle=handles.linestyle;
  clevels=handles.clevels;
  res=handles.res;
  nsampperbox=handles.nsampperbox;
  plottype=handles.plottype;
  fontsize=handles.fontsize;
  markersize=handles.markersize;
  colormap=handles.colormap;
  isheld=handles.isheld;
  save cosmologui_prefs xcol ycol zcol ...
      linewidth linecolor linestyle clevels res ...
      nsampperbox plottype fontsize markersize colormap isheld
  end
      
status=close(gcf);
if (status==0) fprintf('Error closing figure\n'); end


% --- Executes on button press in RotateButton.
function RotateButton_Callback(hObject, eventdata, handles)
% hObject    handle to RotateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

figure(handles.figno)

if (handles.zcol>1)
  fprintf(1,'for i=1:5:360 view([(-20+i) 30]); pause(0.03); end\n')
  for i=1:5:360 
      figure(handles.figno)
      view([(-20+i) 30]); 
      pause(0.05); 
  end
end

% --- Executes on key press over cosmologui with no controls selected.
function cosmologui_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to cosmologui (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --- Executes on button press in SaveMovieButton.
function SaveMovieButton_Callback(hObject, eventdata, handles)
% hObject    handle to SaveMovieButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.zcol>1)
    
    figure(handles.figno)
    
    % make an avi file see help avifile
    set(gcf,'Color',[1 1 1]); % want a white background for ppt movie
    axis vis3d;
    for i=1:round(360/5)
        figure(handles.figno)
        view([(-20+i*5) 30]); 
        %F(i)=getframe(gcf,[1 1 500 380]);
        F(i)=getframe(gcf);
        pause(0.01); 
    end
    %movie(F)
    ! rm cosmologui.avi
    fprintf(1,'Compressing and saving movie to cosmologui.avi ...')
    if (ispc)
        movie2avi(F,'cosmologui.avi','compression','cinepak','fps',5,'quality',75)
    else
        movie2avi(F,'cosmologui.avi','fps',5); % unix/linux versions doesn't compress
    end    
    
    % set back to the way it was
    set(gcf,'Color',get(0,'defaultFigureColor'));
    
end


% --- Executes on button press in HoldRadio.
function HoldRadio_Callback(hObject, eventdata, handles)
% hObject    handle to HoldRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of HoldRadio

handles.isheld=get(handles.HoldRadio,'Value');
guidata(hObject,handles); % store the changes


% --------------------------------------------------------------------
function EditMenu_Callback(hObject, eventdata, handles)
% hObject    handle to CutsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function CutsMenu_Callback(hObject, eventdata, handles)
% hObject    handle to CutsMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% choose some cuts using the cutsamps gui
[handles.cuts h]=cutsamps_devel(handles.samples,handles.names,handles.cuts,...
    'Position',[20 20 100 21]);

% shut down the cutsamps gui
if (ishandle(h))
    status=close(h);
    if (status==0) fprintf(1,'Error closing figure\n'); end
else
    fprintf(1,'Error choosing cuts\n');
    return;
end
if (status==0) fprintf(1,'Error closing figure\n'); end

plotsamps(handles);

guidata(hObject,handles); % store the changes


% --- Executes on button press in PushPriors.
function PushPriors_Callback(hObject, eventdata, handles)
% hObject    handle to PushPriors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% choose some cuts using the cutsamps gui
[handles.cuts h]=cutsamps_devel(handles.samples,handles.names,handles.cuts,...
    'Position',[20 20 100 21]);

% shut down the cutsamps gui
if (ishandle(h))
    status=close(h);
    if (status==0) fprintf(1,'Error closing figure\n'); end
else
    fprintf(1,'Error choosing cuts\n');
    return;
end
if (status==0) fprintf(1,'Error closing figure\n'); end

plotsamps(handles);

guidata(hObject,handles); % store the changes


% --- Executes on button press in PrefsButton.
function PrefsButton_Callback(hObject, eventdata, handles)
% hObject    handle to PrefsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


