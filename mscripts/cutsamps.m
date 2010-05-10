function varargout = cutsamps(varargin)
% CUTSAMPS M-file for cutsamps.fig
%      CUTSAMPS, by itself, creates a new CUTSAMPS or raises the existing
%      singleton*.
%
%      H = CUTSAMPS returns the handle to a new CUTSAMPS or the handle to
%      the existing singleton*.
%
%      CUTSAMPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CUTSAMPS.M with the given input arguments.
%
%      CUTSAMPS('Property','Value',...) creates a new CUTSAMPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cutsamps_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cutsamps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cutsamps

% Last Modified by GUIDE v2.5 04-Oct-2004 16:09:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cutsamps_OpeningFcn, ...
                   'gui_OutputFcn',  @cutsamps_OutputFcn, ...
                   'gui_LayoutFcn',  @cutsamps_LayoutFcn, ...
                   'gui_Callback',   []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before cutsamps is made visible.
function cutsamps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cutsamps (see VARARGIN)

% sort out the input args
if ~isempty(varargin)
    handles.samples=varargin{1}; % messy.. copying around big arrays?
    handles.names=varargin{2};
    handles.cuts=varargin{3};
else    
  % complain
  fprintf(1,'ERROR: cutsamps was called with no inputs\n')
  handles.output={[],hObject};
  guidata(hObject, handles);
  return
end

% find out which parameters have any range worth cutting on
sz=size(handles.samples);
handles.irange=[];
for i=1:sz(1)
    if (max(handles.samples(i,:))>min(handles.samples(i,:)))
        handles.irange=[handles.irange i];
    end
end

% set up the various options
handles.priortypes={'Top Hat'}; % a bit boring at the mo'

%% fill up the various menus
% pulldown list of parameter names 
tmpcell=handles.names(handles.irange);
set(handles.ParameterMenu,'String',tmpcell);
set(handles.PriortypeMenu,'String',handles.priortypes);
set(handles.ParameterMenu,'Value',1);
handles=dopriormenus(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes cutsamps wait for user response (see UIRESUME)
% uiwait(handles.figure1);
uiwait; % don't let it assign the output arguments until we've got them


% --- Outputs from this function are returned to the command line.
function varargout = cutsamps_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.cuts;
varargout{2} = hObject;


% --- Executes during object creation, after setting all properties.
function ParameterMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ParameterMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in ParameterMenu.
function ParameterMenu_Callback(hObject, eventdata, handles)
% hObject    handle to ParameterMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns ParameterMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ParameterMenu

handles=dopriormenus(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function Par1Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Par1Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function Par1Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Par1Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.Par1ValText,'String',num2str(get(handles.Par1Slider,'Value')))
set(handles.Par2Slider,'Min',get(handles.Par1Slider,'Value'));
set(handles.Par2MinText,'String',num2str(get(handles.Par2Slider,'Min')));
itmp=handles.irange(get(handles.ParameterMenu,'Value'));
handles.cuts(2,itmp)=get(handles.Par1Slider,'Value');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function Par2Slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Par2Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background, change
%       'usewhitebg' to 0 to use default.  See ISPC and COMPUTER.
usewhitebg = 1;
if usewhitebg
    set(hObject,'BackgroundColor',[.9 .9 .9]);
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function Par2Slider_Callback(hObject, eventdata, handles)
% hObject    handle to Par2Slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.Par2ValText,'String',num2str(get(handles.Par2Slider,'Value')))
set(handles.Par1Slider,'Max',get(handles.Par2Slider,'Value'));
set(handles.Par1MaxText,'String',num2str(get(handles.Par1Slider,'Max')));
itmp=handles.irange(get(handles.ParameterMenu,'Value'));
handles.cuts(3,itmp)=get(handles.Par2Slider,'Value');

guidata(hObject, handles);


% --- Executes on button press in ListButton.
function ListButton_Callback(hObject, eventdata, handles)
% hObject    handle to ListButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in RemoveAllButton.
function RemoveAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to RemoveAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% go back to defaults
sz=size(handles.samples);
for i=1:sz(1)
    handles.cuts(1,i)=1;
    handles.cuts(2,i)= min(handles.samples(i,:));
    handles.cuts(3,i)= max(handles.samples(i,:));
end

handles=dopriormenus(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function PriortypeMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to PriortypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on selection change in PriortypeMenu.
function PriortypeMenu_Callback(hObject, eventdata, handles)
% hObject    handle to PriortypeMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns PriortypeMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from PriortypeMenu


% --- Executes on button press in CancelButton.
function CancelButton_Callback(hObject, eventdata, handles)
% hObject    handle to CancelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf)

% --- a function called by the others
function handles=dopriormenus(handles);

set(handles.PriortypeMenu,'Value',...
    handles.cuts(1,handles.irange(get(handles.ParameterMenu,'Value'))));

if (get(handles.PriortypeMenu,'Value')==1)
  handles.priorparnames={'Minimum','Maximum'};
  itmp=handles.irange(get(handles.ParameterMenu,'Value'));
  mintmp=min(handles.samples(itmp,:));
  maxtmp=max(handles.samples(itmp,:));
  set(handles.Par1Slider,'Min',mintmp);
  set(handles.Par1Slider,'Max',maxtmp);
  set(handles.Par2Slider,'Min',mintmp);
  set(handles.Par2Slider,'Max',maxtmp);
  set(handles.Par1Slider,'Value',handles.cuts(2,itmp));
  set(handles.Par2Slider,'Value',handles.cuts(3,itmp));
  % parameter 1
  set(handles.Par1Slider,'Min',mintmp);
  set(handles.Par1Slider,'Max',get(handles.Par2Slider,'Value'));
  % parameter 2
  set(handles.Par2Slider,'Min',get(handles.Par1Slider,'Value'));
  set(handles.Par2Slider,'Max',maxtmp);
end
set(handles.Par1NameText,'String',handles.priorparnames{1});
set(handles.Par2NameText,'String',handles.priorparnames{2});
set(handles.Par1MinText,'String',num2str(get(handles.Par1Slider,'Min')));
set(handles.Par1MaxText,'String',num2str(get(handles.Par1Slider,'Max')));
set(handles.Par1ValText,'String',num2str(get(handles.Par1Slider,'Value')));
set(handles.Par2MinText,'String',num2str(get(handles.Par2Slider,'Min')));
set(handles.Par2MaxText,'String',num2str(get(handles.Par2Slider,'Max')));
set(handles.Par2ValText,'String',num2str(get(handles.Par2Slider,'Value')));


% --- Executes on button press in DoneButton.
function DoneButton_Callback(hObject, eventdata, handles)
% hObject    handle to DoneButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

uiresume; % OK, can now assign the output arguments now we've finished


% --- Creates and returns a handle to the GUI figure. 
function h1 = cutsamps_LayoutFcn(policy)
% policy - create a new figure or use a singleton. 'new' or 'reuse'.

persistent hsingleton;
if strcmpi(policy, 'reuse') & ishandle(hsingleton)
    h1 = hsingleton;
    return;
end

appdata = [];
appdata.GUIDEOptions = struct(...
    'active_h', [], ...
    'taginfo', struct(...
    'figure', 2, ...
    'text', 12, ...
    'popupmenu', 3, ...
    'slider', 4, ...
    'pushbutton', 5), ...
    'override', 0, ...
    'release', 13, ...
    'resize', 'none', ...
    'accessibility', 'callback', ...
    'mfile', 1, ...
    'callbacks', 1, ...
    'singleton', 1, ...
    'syscolorfig', 1, ...
    'blocking', 0, ...
    'lastSavedFile', 'D:\home\sarah\work\joint\cosmomc\mscripts\cutsamps.m');
appdata.lastValidTag = 'cutsamps_devel';
appdata.GUIDELayoutEditor = [];

h1 = figure(...
'Units','characters',...
'PaperUnits',get(0,'defaultfigurePaperUnits'),...
'Color',[0.87843137254902 0.874509803921569 0.890196078431373],...
'Colormap',[0 0 0.5625;0 0 0.625;0 0 0.6875;0 0 0.75;0 0 0.8125;0 0 0.875;0 0 0.9375;0 0 1;0 0.0625 1;0 0.125 1;0 0.1875 1;0 0.25 1;0 0.3125 1;0 0.375 1;0 0.4375 1;0 0.5 1;0 0.5625 1;0 0.625 1;0 0.6875 1;0 0.75 1;0 0.8125 1;0 0.875 1;0 0.9375 1;0 1 1;0.0625 1 1;0.125 1 0.9375;0.1875 1 0.875;0.25 1 0.8125;0.3125 1 0.75;0.375 1 0.6875;0.4375 1 0.625;0.5 1 0.5625;0.5625 1 0.5;0.625 1 0.4375;0.6875 1 0.375;0.75 1 0.3125;0.8125 1 0.25;0.875 1 0.1875;0.9375 1 0.125;1 1 0.0625;1 1 0;1 0.9375 0;1 0.875 0;1 0.8125 0;1 0.75 0;1 0.6875 0;1 0.625 0;1 0.5625 0;1 0.5 0;1 0.4375 0;1 0.375 0;1 0.3125 0;1 0.25 0;1 0.1875 0;1 0.125 0;1 0.0625 0;1 0 0;0.9375 0 0;0.875 0 0;0.8125 0 0;0.75 0 0;0.6875 0 0;0.625 0 0;0.5625 0 0],...
'IntegerHandle','off',...
'InvertHardcopy',get(0,'defaultfigureInvertHardcopy'),...
'MenuBar','none',...
'Name','cutsamps_devel',...
'NumberTitle','off',...
'PaperPosition',get(0,'defaultfigurePaperPosition'),...
'PaperSize',[20.98404194812 29.67743169791],...
'PaperType',get(0,'defaultfigurePaperType'),...
'Position',[103.8 42.4615384615385 105.6 19],...
'Renderer',get(0,'defaultfigureRenderer'),...
'RendererMode','manual',...
'Resize','off',...
'HandleVisibility','callback',...
'Tag','cutsamps_devel',...
'UserData',[],...
'Behavior',get(0,'defaultfigureBehavior'),...
'Visible','on',...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'TitleText';

h2 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',14,...
'ForegroundColor',[0 0.501960784313725 0],...
'ListboxTop',0,...
'Position',[10.2 16.3076923076923 47.8 1.76923076923077],...
'String','Add or Modify Cuts or Priors',...
'Style','text',...
'Tag','TitleText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'ParameterMenu';

h3 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','cutsamps(''ParameterMenu_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[35.8 13.3846153846154 48.2 1.61538461538462],...
'String','List of pars should appear here',...
'Style','popupmenu',...
'TooltipString','Choose the parameter to cut on',...
'Value',1,...
'CreateFcn', {@local_CreateFcn, 'cutsamps(''ParameterMenu_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','ParameterMenu',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'Par1NameText';

h4 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ForegroundColor',[0 0.501960784313725 0],...
'ListboxTop',0,...
'Position',[3.8 9.15384615384616 12.2 1.23076923076923],...
'String','Par1Name',...
'Style','text',...
'Tag','Par1NameText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'Par2NameText';

h5 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ForegroundColor',[0 0.501960784313725 0],...
'ListboxTop',0,...
'Position',[3.8 6.07692307692308 12.2 1.23076923076923],...
'String','Par2Name',...
'Style','text',...
'Tag','Par2NameText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'Par1Slider';

h6 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback','cutsamps(''Par1Slider_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[35.8 8.76923076923077 48.2 1.69230769230769],...
'String',{  '' },...
'Style','slider',...
'CreateFcn', {@local_CreateFcn, 'cutsamps(''Par1Slider_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','Par1Slider',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'Par2Slider';

h7 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[0.9 0.9 0.9],...
'Callback','cutsamps(''Par2Slider_Callback'',gcbo,[],guidata(gcbo))',...
'ListboxTop',0,...
'Position',[35.8 5.6923076923077 48.2 1.61538461538462],...
'String',{  '' },...
'Style','slider',...
'CreateFcn', {@local_CreateFcn, 'cutsamps(''Par2Slider_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','Par2Slider',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'Par1MinText';

h8 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[35.8 7.61538461538462 16.2 1.23076923076923],...
'String','Par1Min',...
'Style','text',...
'Tag','Par1MinText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'Par1MaxText';

h9 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[67.8 7.61538461538462 16.4 1.23076923076923],...
'String','Par1Max',...
'Style','text',...
'Tag','Par1MaxText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'Par2ValText';

h10 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[91.8 6.07692307692308 10.2 1.23076923076923],...
'String','Par2Val',...
'Style','text',...
'Tag','Par2ValText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'Par2MinText';

h11 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[35.8 4.53846153846154 9.6 1.23076923076923],...
'String','Par2Min',...
'Style','text',...
'Tag','Par2MinText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'Par2MaxText';

h12 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[73.8 4.53846153846154 10.4 1.23076923076923],...
'String','Par2Max',...
'Style','text',...
'Tag','Par2MaxText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'Par1ValText';

h13 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[89.8 9.15384615384616 14.2 1.23076923076923],...
'String','Par1Val',...
'Style','text',...
'Tag','Par1ValText',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'ListButton';

h14 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','cutsamps(''ListButton_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[5.8 0.923076923076927 13.2 1.76923076923077],...
'String','List All',...
'TooltipString','List all the cuts and priors in the main matlab window.',...
'Tag','ListButton',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'text10';

h15 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ForegroundColor',[0 0.501960784313725 0],...
'ListboxTop',0,...
'Position',[3.8 13.7692307692308 30.2 1.23076923076923],...
'String','Modify cuts for parameter:',...
'Style','text',...
'TooltipString','Choose the parameter to change the cuts for.',...
'Tag','text10',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'RemoveAllButton';

h16 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','cutsamps(''RemoveAllButton_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[21.8 0.846153846153847 14.6 1.84615384615385],...
'String','Clear All',...
'TooltipString','Remove all cuts and priors.',...
'Tag','RemoveAllButton',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'text11';

h17 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'FontSize',10,...
'ForegroundColor',[0 0.501960784313725 0],...
'ListboxTop',0,...
'Position',[3.8 11.8461538461538 11.8 1.23076923076923],...
'String','Prior Type',...
'Style','text',...
'TooltipString','Decides the type of cut or prior. Top hat is the same as a simple cut.',...
'Tag','text11',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'PriortypeMenu';

h18 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'BackgroundColor',[1 1 1],...
'Callback','cutsamps(''PriortypeMenu_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[35.8 11.0769230769231 48.2 1.61538461538462],...
'String','Types of priors should appear here',...
'Style','popupmenu',...
'Value',1,...
'CreateFcn', {@local_CreateFcn, 'cutsamps(''PriortypeMenu_CreateFcn'',gcbo,[],guidata(gcbo))', appdata} ,...
'Tag','PriortypeMenu',...
'Behavior',get(0,'defaultuicontrolBehavior'));

appdata = [];
appdata.lastValidTag = 'CancelButton';

h19 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','cutsamps(''CancelButton_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[87.8 1 13 1.76923076923077],...
'String','Cancel',...
'TooltipString','Cancel modifying cuts or priors',...
'Tag','CancelButton',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );

appdata = [];
appdata.lastValidTag = 'DoneButton';

h20 = uicontrol(...
'Parent',h1,...
'Units','characters',...
'Callback','cutsamps(''DoneButton_Callback'',gcbo,[],guidata(gcbo))',...
'FontSize',10,...
'ListboxTop',0,...
'Position',[60.8 1 13.2 1.76923076923077],...
'String','Done',...
'Tag','DoneButton',...
'Behavior',get(0,'defaultuicontrolBehavior'),...
'CreateFcn', {@local_CreateFcn, '', appdata} );


hsingleton = h1;


% --- Set application data first then calling the CreateFcn. 
function local_CreateFcn(hObject, eventdata, createfcn, appdata)

if ~isempty(appdata)
   names = fieldnames(appdata);
   for i=1:length(names)
       name = char(names(i));
       setappdata(hObject, name, getfield(appdata,name));
   end
end

if ~isempty(createfcn)
   eval(createfcn);
end


% --- Handles default GUIDE GUI creation and callback dispatch
function varargout = gui_mainfcn(gui_State, varargin)


%   GUI_MAINFCN provides these command line APIs for dealing with GUIs
%
%      CUTSAMPS, by itself, creates a new CUTSAMPS or raises the existing
%      singleton*.
%
%      H = CUTSAMPS returns the handle to a new CUTSAMPS or the handle to
%      the existing singleton*.
%
%      CUTSAMPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CUTSAMPS.M with the given input arguments.
%
%      CUTSAMPS('Property','Value',...) creates a new CUTSAMPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before untitled_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to untitled_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".

%   Copyright 1984-2004 The MathWorks, Inc.
%   $Revision: 1.4.6.8 $ $Date: 2004/04/15 00:06:57 $

gui_StateFields =  {'gui_Name'
                    'gui_Singleton'
                    'gui_OpeningFcn'
                    'gui_OutputFcn'
                    'gui_LayoutFcn'
                    'gui_Callback'};
gui_Mfile = '';
for i=1:length(gui_StateFields)
    if ~isfield(gui_State, gui_StateFields{i})
        error('Could not find field %s in the gui_State struct in GUI M-file %s', gui_StateFields{i}, gui_Mfile);        
    elseif isequal(gui_StateFields{i}, 'gui_Name')
        gui_Mfile = [gui_State.(gui_StateFields{i}), '.m'];
    end
end

numargin = length(varargin);

if numargin == 0
    % CUTSAMPS
    % create the GUI
    gui_Create = 1;
elseif isequal(ishandle(varargin{1}), 1) && ispc && iscom(varargin{1}) && isequal(varargin{1},gcbo)
    % CUTSAMPS(ACTIVEX,...)    
    vin{1} = gui_State.gui_Name;
    vin{2} = [get(varargin{1}.Peer, 'Tag'), '_', varargin{end}];
    vin{3} = varargin{1};
    vin{4} = varargin{end-1};
    vin{5} = guidata(varargin{1}.Peer);
    feval(vin{:});
    return;
elseif ischar(varargin{1}) && numargin>1 && isequal(ishandle(varargin{2}), 1)
    % CUTSAMPS('CALLBACK',hObject,eventData,handles,...)
    gui_Create = 0;
else
    % CUTSAMPS(...)
    % create the GUI and hand varargin to the openingfcn
    gui_Create = 1;
end

if gui_Create == 0
    varargin{1} = gui_State.gui_Callback;
    if nargout
        [varargout{1:nargout}] = feval(varargin{:});
    else
        feval(varargin{:});
    end
else
    if gui_State.gui_Singleton
        gui_SingletonOpt = 'reuse';
    else
        gui_SingletonOpt = 'new';
    end
    
    % Open fig file with stored settings.  Note: This executes all component
    % specific CreateFunctions with an empty HANDLES structure.
    
    % Do feval on layout code in m-file if it exists
    if ~isempty(gui_State.gui_LayoutFcn)
        gui_hFigure = feval(gui_State.gui_LayoutFcn, gui_SingletonOpt);
        % openfig (called by local_openfig below) does this for guis without
        % the LayoutFcn. Be sure to do it here so guis show up on screen.
        movegui(gui_hFigure,'onscreen')
    else
        gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        % If the figure has InGUIInitialization it was not completely created
        % on the last pass.  Delete this handle and try again.
        if isappdata(gui_hFigure, 'InGUIInitialization')
            delete(gui_hFigure);
            gui_hFigure = local_openfig(gui_State.gui_Name, gui_SingletonOpt);            
        end
    end
    
    % Set flag to indicate starting GUI initialization
    setappdata(gui_hFigure,'InGUIInitialization',1);

    % Fetch GUIDE Application options
    gui_Options = getappdata(gui_hFigure,'GUIDEOptions');
    
    if ~isappdata(gui_hFigure,'GUIOnScreen')
        % Adjust background color
        if gui_Options.syscolorfig 
            set(gui_hFigure,'Color', get(0,'DefaultUicontrolBackgroundColor'));
        end

        % Generate HANDLES structure and store with GUIDATA
        guidata(gui_hFigure, guihandles(gui_hFigure));
    end
    
    % If user specified 'Visible','off' in p/v pairs, don't make the figure
    % visible.
    gui_MakeVisible = 1;
    for ind=1:2:length(varargin)
        if length(varargin) == ind
            break;
        end
        len1 = min(length('visible'),length(varargin{ind}));
        len2 = min(length('off'),length(varargin{ind+1}));
        if ischar(varargin{ind}) && ischar(varargin{ind+1}) && ...
                strncmpi(varargin{ind},'visible',len1) && len2 > 1
            if strncmpi(varargin{ind+1},'off',len2)
                gui_MakeVisible = 0;
            elseif strncmpi(varargin{ind+1},'on',len2)
                gui_MakeVisible = 1;
            end
        end
    end
    
    % Check for figure param value pairs
    for index=1:2:length(varargin)
        if length(varargin) == index
            break;
        end
        try set(gui_hFigure, varargin{index}, varargin{index+1}), catch break, end
    end

    % If handle visibility is set to 'callback', turn it on until finished
    % with OpeningFcn
    gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
    if strcmp(gui_HandleVisibility, 'callback')
        set(gui_hFigure,'HandleVisibility', 'on');
    end
    
    feval(gui_State.gui_OpeningFcn, gui_hFigure, [], guidata(gui_hFigure), varargin{:});
    
    if ishandle(gui_hFigure)
        % Update handle visibility
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
        
        % Make figure visible
        if gui_MakeVisible
            set(gui_hFigure, 'Visible', 'on')
            if gui_Options.singleton 
                setappdata(gui_hFigure,'GUIOnScreen', 1);
            end
        end

        % Done with GUI initialization
        rmappdata(gui_hFigure,'InGUIInitialization');
    end
    
    % If handle visibility is set to 'callback', turn it on until finished with
    % OutputFcn
    if ishandle(gui_hFigure)
        gui_HandleVisibility = get(gui_hFigure,'HandleVisibility');
        if strcmp(gui_HandleVisibility, 'callback')
            set(gui_hFigure,'HandleVisibility', 'on');
        end
        gui_Handles = guidata(gui_hFigure);
    else
        gui_Handles = [];
    end
    
    if nargout
        [varargout{1:nargout}] = feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    else
        feval(gui_State.gui_OutputFcn, gui_hFigure, [], gui_Handles);
    end
    
    if ishandle(gui_hFigure)
        set(gui_hFigure,'HandleVisibility', gui_HandleVisibility);
    end
end    

function gui_hFigure = local_openfig(name, singleton)

% openfig with three arguments was new from R13. Try to call that first, if
% failed, try the old openfig.
try 
    gui_hFigure = openfig(name, singleton, 'auto');
catch
    % OPENFIG did not accept 3rd input argument until R13,
    % toggle default figure visible to prevent the figure
    % from showing up too soon.
    gui_OldDefaultVisible = get(0,'defaultFigureVisible');
    set(0,'defaultFigureVisible','off');
    gui_hFigure = openfig(name, singleton);
    set(0,'defaultFigureVisible',gui_OldDefaultVisible);
end

