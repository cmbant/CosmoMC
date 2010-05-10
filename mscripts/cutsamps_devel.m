function varargout = cutsamps_devel(varargin)
% CUTSAMPS_DEVEL M-file for cutsamps_devel.fig
%      CUTSAMPS_DEVEL, by itself, creates a new CUTSAMPS_DEVEL or raises the existing
%      singleton*.
%
%      H = CUTSAMPS_DEVEL returns the handle to a new CUTSAMPS_DEVEL or the handle to
%      the existing singleton*.
%
%      CUTSAMPS_DEVEL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CUTSAMPS_DEVEL.M with the given input arguments.
%
%      CUTSAMPS_DEVEL('Property','Value',...) creates a new CUTSAMPS_DEVEL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before cutsamps_devel_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to cutsamps_devel_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help cutsamps_devel

% Last Modified by GUIDE v2.5 08-Jul-2004 20:14:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @cutsamps_devel_OpeningFcn, ...
                   'gui_OutputFcn',  @cutsamps_devel_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
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


% --- Executes just before cutsamps_devel is made visible.
function cutsamps_devel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to cutsamps_devel (see VARARGIN)

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

% UIWAIT makes cutsamps_devel wait for user response (see UIRESUME)
% uiwait(handles.figure1);
uiwait; % don't let it assign the output arguments until we've got them


% --- Outputs from this function are returned to the command line.
function varargout = cutsamps_devel_OutputFcn(hObject, eventdata, handles)
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
