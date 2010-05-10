function varargout = getsamps(varargin)
% GETSAMPS M-file for getsamps.fig
%      GETSAMPS, by itself, creates a new GETSAMPS or raises the existing
%      singleton*.
%
%      [SAMPLES, NAMES, ISAMPLED] = GETSAMPS returns the burn-removed samples info
%      SAMPLES, with column names NAMES as specified when running GETSAMPS
%
%      SAMPLES is just the burn-removed contents of the samples files.
%      NAMES is a cell array containing the column names
%      ISAMPLED is an array of the column numbers worth plotting
%
%      GETSAMPS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GETSAMPS.M with the given input arguments.
%
%      GETSAMPS('Property','Value',...) creates a new GETSAMPS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before getsamps_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to getsamps_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help getsamps

% Last Modified by GUIDE v2.5 16-Oct-2004 11:06:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @getsamps_OpeningFcn, ...
                   'gui_OutputFcn',  @getsamps_OutputFcn, ...
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


% --- Executes just before getsamps is made visible.
function getsamps_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to getsamps (see VARARGIN)

% set default handles associated with reading in samples file
handles.sampsdir= '../chains/'; % default samples directory
handles.namesdir= './'; % default names directory
handles.namesfile= 'CMB.names'; % default names filename
handles.nburn= 100; % default number to burn
handles.nthin=10; % default thinning factor
samplesfiles=''; % ??
handles.colors = [0 0 1]; % default colors

fid=fopen('getsamps_prefs.mat');
if (fid>0)
    fclose(fid);
    htmp=load('getsamps_prefs'); % sampsdir namesdir namesfile nburn nthin filescell
    % NB. loading the above into handles wipes the handles struct !
    % is there a way to add elements to a struct from another struct?
    % ie. handles= concatstruct(handles,htmp) or something? 
    % maybe could write one?
    handles.sampsdir=htmp.sampsdir;
    handles.namesdir=htmp.namesdir;
    handles.namesfile=htmp.namesfile;
    handles.nburn=htmp.nburn;
    handles.nthin=htmp.nthin;
    samplesfiles=htmp.samplesfiles;
    usefiles=htmp.usefiles;
end
    
set(handles.NamesButton,'String',[handles.namesfile]);
set(handles.NoBurnEdit,'String',num2str(handles.nburn));
set(handles.NoThinEdit,'String',num2str(handles.nthin));
set(handles.SamplesFilesListbox,'String',samplesfiles);
if (exist('usefiles'))
    set(handles.SamplesFilesListbox,'Value',usefiles);
end

% set default output values 
handles.samples = -1; % the contents of the samples file (ie. nothing yet..)
handles.names = 'Not known yet'; % the names of the columns (dummy)
handles.isampled= -1; % the columns worth plotting
handles.cuts= []; % default cuts

% Update handles structure
guidata(hObject, handles);

%% bring up the file browser immediately - will definitely want it!
%AddSamplesFilesButton_Callback(hObject, eventdata, handles)
% can't do this here, it seems, as this confuses it..

% UIWAIT makes getsamps wait for user response (see UIRESUME)
% uiwait(handles.getsamps);
uiwait; % don't let it assign the output arguments until we've got them


% --- Outputs from this function are returned to the command line.
function varargout = getsamps_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if (isempty(handles))
    varargout{1} = -1;
    varargout{2} = 'Not known yet';
    varargout{3} = -1;
    varargout{4}=-1;
    varargout{5}=-1;
else    
    varargout{1} = handles.samples;
    varargout{2} = handles.names;
    varargout{3} = handles.isampled;
    varargout{4} = handles.cuts;
    if sum(sum(handles.samples~=-1))
        varargout{5} = gcf;
    else
        varargout{5}=-1;
    end
    varargout{6} = handles.colors;
end

    
% --- Executes on button press in AddSamplesFilesButton.
function AddSamplesFilesButton_Callback(hObject, eventdata, handles)
% hObject    handle to AddSamplesFilesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile([handles.sampsdir,'*.txt']);
if ~isequal(file, 0)
    %open(file);
    handles.sampsdir=path; % probs will want to use this next time too
    oldfiles=get(handles.SamplesFilesListbox,'String');
    newfiles=strvcat(oldfiles,[path,file]);
    tmp=size(newfiles);
    nfiles=tmp(1);
    if (nfiles>get(handles.SamplesFilesListbox,'Max'))
        set(handles.SamplesFilesListbox,'Max',nfiles+10); 
        % why not just increasing as reqd.
%        fprintf(1,'Error: need to increase Max for Listbox SamplesFileList')
    end
    set(handles.SamplesFilesListbox,'String',newfiles);
    % highlight the new entry as well as previously highlighted entries
    set(handles.SamplesFilesListbox,'Value',...
        [get(handles.SamplesFilesListbox,'Value') nfiles]);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function SamplesFilesListbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SamplesFilesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on selection change in SamplesFilesListbox.
function SamplesFilesListbox_Callback(hObject, eventdata, handles)
% hObject    handle to SamplesFilesListbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns SamplesFilesListbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from SamplesFilesListbox


% --- Executes on button press in NamesButton.
function NamesButton_Callback(hObject, eventdata, handles)
% hObject    handle to NamesButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file,path] = uigetfile([handles.namesdir,'*.names']);
if ~isequal(file, 0)
    handles.namesdir=path;
    handles.namesfile=file;
    set(handles.NamesButton,'String',[handles.namesfile]);
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function NoBurnEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NoBurnEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

% --- Executes on changing number in NoBurnEdit box
function NoBurnEdit_Callback(hObject, eventdata, handles)
% do nothing

% --- Executes on button press in LoadButton.
function LoadButton_Callback(hObject, eventdata, handles)
% hObject    handle to LoadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.nburn=str2double(get(handles.NoBurnEdit,'String'));
handles.nthin=str2double(get(handles.NoThinEdit,'String'));
handles.filescell=cellstr(get(handles.SamplesFilesListbox,'String'));

% get rid of repetitions in 
% get(handles.SamplesFilesListbox,'Value')
% why does this happen? Surely we would never want it.
vals=get(handles.SamplesFilesListbox,'Value');
iin=1;
for i=2:length(vals)
    if ~(vals(i)==vals(i-1))
        iin=[iin i];
    end
end

goodfile=0;
while (goodfile~=1)
[handles.samples,handles.names]=readsampsf( ...
    handles.filescell(vals(iin)), ...
    [handles.namesdir,handles.namesfile],...
    handles.nburn,handles.nthin);
    if isempty(handles.samples)
        fprintf('ERROR: No samples were stored\n')
        return;
    end
    if (length(handles.samples)>10)
        goodfile=1;
    else
        fprintf(1,'ERROR: File contains less than 10 lines!\n')
        return
%        AddSamplesfile_Callback(hObject,eventdata,handles);
    end    
end

% find out which columns are sampled from and which are fixed
handles.isampled=find(min(handles.samples,[],2)-max(handles.samples,[],2));
% we don't want to include the second column
% (also will change meaning of request for first column)
tmp=find(handles.isampled>=3);
handles.isampled=[1 handles.isampled(tmp)'];
%handles.isampled=handles.isampled([1 3:length(handles.isampled)]);
% don't include any columns that don't have names
tmp=size(handles.names);
handles.isampled=handles.isampled(find(handles.isampled<=tmp(1)));

handles=setdefaultcuts(handles);
guidata(hObject,handles); % store the changes

% save the inputs to a file for next-time
% annoying.. I can't put bits of structures to save..
% how can I do the below (and corresponding read in..) better?
sampsdir=handles.sampsdir;
namesdir=handles.namesdir;
namesfile=handles.namesfile;
nburn=handles.nburn;
nthin=handles.nthin;
samplesfiles=get(handles.SamplesFilesListbox,'String');
usefiles=get(handles.SamplesFilesListbox,'Value');
save getsamps_prefs sampsdir namesdir namesfile nburn nthin samplesfiles usefiles

uiresume; % OK, can now assign the output arguments now we've finished


% --- Executes on key press over getsamps with no controls selected.
function getsamps_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to getsamps (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if a button is pressed at any time within this window, just make this
% equivalent to ...

% hmm, what if this window isn't called getsamps ??

cchar=get(gcf,'CurrentCharacter');
fprintf(1,['currentchar=',cchar,'\n'])
if (cchar=='q')
    close(gcf)
elseif (cchar=='a')
    AddSamplesFilesButton_Callback(hObject, eventdata, handles)
elseif (cchar=='l')
    LoadButton_Callback(hObject, eventdata, handles)
end

% whatever character is pressed, load the data
% (would ideally like this to happen only if `return' is pressed..
%  but how to check for this character?)
%LoadButton_Callback(hObject, eventdata, handles)


% --- Executes on button press in RemovePush.
function RemovePush_Callback(hObject, eventdata, handles)
% hObject    handle to RemovePush (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

samplesfiles=get(handles.SamplesFilesListbox,'String');
tmp=size(samplesfiles);
nfiles=tmp(1);

% remove all the highlighted files from the list
removefiles=get(handles.SamplesFilesListbox,'Value');
usefiles=[];
for i=1:nfiles
    if (~sum(i==removefiles))
        usefiles=[usefiles i]
    end
end

set(handles.SamplesFilesListbox,'String',samplesfiles(usefiles,:))
set(handles.SamplesFilesListbox,'Value',1); % ? anything more sensible?


% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

close(gcf)


% -------------------------------------------------------------------
% --- Called by functions above, to set cuts to default values
% esp. after reading in new data
function handles=setdefaultcuts(handles);

% Cuts are stored in handles.cuts
% format is
% handles.cuts(1,i) is the type of cut for parameter handles.isampled(i)
%                   cut type = 1 is just a top hat range
% handles.cuts(2,i) is the first parameter for this cut type 
%                   for cut type = 1 this is the minimum of the range
% handles.cuts(3,i) is the second parameter for this cut type
%                   for cut type = 1 this is the maximum of the range
% handles.cuts(j,i) is the j-1 th parameter for this cut type..
%
% In this subroutine just set to the default, ie. no cuts
% which is equivalent to handles.cuts(1,i)=1
% handles.cuts(2,i) = minimum value of this parameter
% handles.cuts(3,i) = maximum value of this parameter
% Set this up now, to make it easier to modify them later.

sz=size(handles.samples);
for i=1:sz(1)
    handles.cuts(1,i)=1;
    handles.cuts(2,i)= min(handles.samples(i,:));
    handles.cuts(3,i)= max(handles.samples(i,:));
end




% --- Executes during object creation, after setting all properties.
function NoThinEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NoThinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end



function NoThinEdit_Callback(hObject, eventdata, handles)
% hObject    handle to NoThinEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NoThinEdit as text
%        str2double(get(hObject,'String')) returns contents of NoThinEdit as a double



function ColorText_Callback(hObject, eventdata, handles)
% hObject    handle to ColorText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ColorText as text
%        str2double(get(hObject,'String')) returns contents of ColorText as a double


% --- Executes during object creation, after setting all properties.
function ColorText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ColorText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end


% --- Executes on slider movement.
function RSlider_Callback(hObject, eventdata, handles)
% hObject    handle to RSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.RColor,'String',num2str(get(handles.RSlider,'Value')))
handles.colors(1)=get(handles.RSlider,'Value');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function RSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to RSlider (see GCBO)
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
function GSlider_Callback(hObject, eventdata, handles)
% hObject    handle to GSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.GColor,'String',num2str(get(handles.GSlider,'Value')))
handles.colors(2)=get(handles.GSlider,'Value');

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function GSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GSlider (see GCBO)
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
function BSlider_Callback(hObject, eventdata, handles)
% hObject    handle to BSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

set(handles.BColor,'String',num2str(get(handles.BSlider,'Value')))
handles.colors(3)=get(handles.BSlider,'Value');

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function BSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to BSlider (see GCBO)
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


