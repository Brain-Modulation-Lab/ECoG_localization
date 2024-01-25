function varargout = dbs_fidlocalizer(varargin)
% DBSFIGURE MATLAB code for dbsfigure.fig
%      DBSFIGURE, by itself, creates a new DBSFIGURE or raises the existing
%      singleton*.
%
%      H = DBSFIGURE returns the handle to a new DBSFIGURE or the handle to
%      the existing singleton*.
%
%      DBSFIGURE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DBSFIGURE.M with the given input arguments.
%
%      DBSFIGURE('Property','Value',...) creates a new DBSFIGURE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dbs_fidlocalizer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dbs_fidlocalizer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dbsfigure

% Last Modified by GUIDE v2.5 11-Jan-2024 14:16:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @dbs_fidlocalizer_OpeningFcn, ...
                   'gui_OutputFcn',  @dbs_fidlocalizer_OutputFcn, ...
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


% --- Executes just before dbsfigure is made visible.
function dbs_fidlocalizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dbsfigure (see VARARGIN)

% Choose default command line output for dbsfigure
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

dbsroot=dbs_getroot;

if ~isdeployed
    addpath(genpath(dbsroot));
end

set(handles.dbsfigure,'name','Fiducial Localizer');

% add recent patients...
load([dbsroot,'dbs_recentpatients.mat']);
if iscell(fullrpts)
    options.uipatdir = fullrpts(1);
end
dbs_load_pts(handles,options.uipatdir);
clear fullrpts

% add recent fsfolders...
load([dbsroot,'dbs_recentfsfolders.mat']);
if iscell(fullrfs)
    options.uifsdir = fullrfs(1);
end
dbs_load_fs(handles,options.uifsdir);
clear fullrfs

%%%%%%%%%%%%%%%%%%%
% set(0,'gca',handles.logoaxes);
set(0,'CurrentFigure',handles.dbsfigure);
im=imread([dbsroot,'icons',filesep,'logo_brainmodulationlab_small.png']);
axes(handles.axes1);
image(im);
axis off;
axis equal;

% %%%%%%%%%%%%
% prefs = dbs_prefs('');
% options.filename = fullfile(options.uipatdir,prefs.fidlocalizer.firstimage); %Default postopmri.nii first
% dbs_fidshowimage(handles,options)

%%%%%%%%%%%%
% set(handles.dbsfigure,'WindowButtonDownFcn',@MouseClick)
% set(handles.dbsfigure,'pointer','cross')
% UIWAIT makes dbsfigure wait for user response (see UIRESUME)
% uiwait(handles.dbsfigure);

% --- Outputs from this function are returned to the command line.
function varargout = dbs_fidlocalizer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes during object creation, after setting all properties.
function fidlocalizer_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fidlocalizer (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in patdir_choosebox.
function patdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dbs_getpatients(handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over patdir_choosebox.
function patdir_choosebox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to patdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dbs_storeui(handles)

% --- Executes on button press in fsdir_choosebox.
function fsdir_choosebox_Callback(hObject, eventdata, handles)
% hObject    handle to fsdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dbs_getfsfolder(handles)

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over fsdir_choosebox.
function fsdir_choosebox_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to fsdir_choosebox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dbs_storeui(handles)

% --- Executes on button press in openfsdir.
function openfsdir_Callback(hObject, eventdata, handles)
% hObject    handle to openfsdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prefs = dbs_prefs('');
outfolder=get(handles.fsdir_choosebox,'String');

if strcmp(outfolder,'No Patient Selected')
    msgbox('Please set the working directory first!', 'Error','error');
    return;
end

if ismac
    system(['open ', outfolder]);
elseif isunix
    system(['xdg-open ', outfolder]);
elseif ispc
    system(['explorer ', outfolder]);
end

if strcmp(prefs.patdir.showpath,'on')
    disp(outfolder)
end
if strcmp(prefs.patdir.cd,'on')
    cd(outfolder);
end

% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over openfsdir.
function openfsdir_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to openfsdir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in openimage.
function openimage_Callback(hObject, eventdata, handles)
% hObject    handle to openimage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% prefs = dbs_prefs('');
options.uipatdir=getappdata(handles.dbsfigure,'uipatdir');
options.uifsdir=getappdata(handles.dbsfigure,'uifsdir');
[filename,path,ext] = dbs_uigetfile(options.uipatdir,'Choose NifTi'); %Default postopmri.nii first
% if ~strcmp(filename(1),'r')
%     error('Please choose a coregistered CT')
% end
options.filename = char(fullfile(path,[filename,ext]));
dbs_fidshowimage(handles,options)
set(handles.dbsfigure,'KeyPressFcn',@KeyPress)
set(handles.dbsfigure,'KeyReleaseFcn', @KeyRelease)
% set(gcf,'WindowButtonDownFcn', @KeyRelease)
% set(gcf,'WindowButtonUpFcn', @KeyRelease)

% --- Executes on button press in saveEls.
function saveEls_Callback(hObject, eventdata, handles)
% hObject    handle to saveEls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in saveFids.
function saveFids_Callback(hObject, eventdata, handles)
% hObject    handle to saveFids (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in viewEls.
function viewEls_Callback(hObject, eventdata, handles)
% hObject    handle to viewEls (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when dbsfigure is resized.
function dbsfigure_SizeChangedFcn(hObject, eventdata, handles)
% hObject    handle to dbsfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% -=< Key Press callback function >=-
function KeyPress(hObject, eventdata, handles)
%     KeyStat = eventdata.Key;
    modifiers = get(gcf,'currentModifier'); 
    if ~isempty(find(ismember({'shift'},modifiers)));
        handles=guidata(hObject);
        %         disp('shift')
        set(handles.dbsfigure,'pointer','cross')
        set(handles.dbsfigure,'WindowButtonDownFcn',@MouseClick)
    end

    
% -=< Mouse button released callback function >=-
function KeyRelease (hObject,eventdata)
        handles=guidata(hObject);
        %         disp('up')
        set(handles.dbsfigure,'pointer','arrow')
        set(handles.dbsfigure,'WindowButtonDownFcn','')

        
% --- Executes on mouse click over figure background.
function MouseClick(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MouseStat = get(gcbf, 'SelectionType');
if sum(strcmp(MouseStat,{'normal','extend'}))     %   Left CLICK
    drawElec
end

handles=guidata(hObject);

% % a=get(hObject,'CurrentPoint');
% % b=get(handles.axes1,'Position');
% % % c=get(handles.axes2,'Position');
% % % d=get(handles.axes3,'Position');
% % % % if handles.View==1
% %     if a(1)>b(1) && a(1)<b(1)+b(3) && a(2)>b(2) && a(2)<b(2)+b(4)
% % 
% %         xy=get(handles.axes1,'CurrentPoint'); %current point x and y
% % 
% %         %check for click within image
% %         if xy(1,2)>5 && xy(1,2)<=length(handles.data.mrData(:,1,1))-5 && ...
% %             xy(1,1)>5 && xy(1,1)<=length(handles.data.mrData(1,:,1))-5 % I do not know for sure whether X and Y are correct here
% %             % if already an electrode in that position remove
% %             if handles.data.elecMap(round(xy(1,2)),round(xy(1,1)),round(get(handles.sl_Z,'Value')))>0;
% %                 removeElec(hObject, eventdata, handles,xy);
% %             else
% %                 drawElec(hObject, eventdata, handles,xy);
% %             end
% %         end
% %     else
% %         %set(handles.display_feedback,'String','mis');
% %     end
% % elseif handles.View==2
% %     if a(1)>b(1) && a(1)<b(1)+b(3) && a(2)>b(2) && a(2)<b(2)+b(4)
% %        xy=get(handles.axes1,'CurrentPoint'); %current point x and y
% %         %check for click within image
% %         if xy(1,2)>5 && xy(1,2)<=length(handles.data.mrData(:,1,1))-5 && ...
% %             xy(1,1)>5 && xy(1,1)<=length(handles.data.mrData(1,1,:))-5 % I do not know for sure whether X and Y are correct here
% %             % if already an electrode in that position remove
% %             if handles.data.elecMap(round(xy(1,2)),round(get(handles.sl_Z,'Value')),round(xy(1,1)))>0;
% %                 %drawElec(hObject, eventdata, handles,xy);
% %                 removeElec(hObject, eventdata, handles,xy);
% %             else
% %                 drawElec(hObject, eventdata, handles,xy);
% %             end
% %         end
% %     else
% %         %set(handles.display_feedback,'String','mis');
% %     end
% % elseif handles.View==3
% %     if a(1)>b(1) && a(1)<b(1)+b(3) && a(2)>b(2) && a(2)<b(2)+b(4)
% %        xy=get(handles.axes1,'CurrentPoint'); %current point x and y
% %         %check for click within image
% %         if xy(1,2)>5 && xy(1,2)<=length(handles.data.mrData(1,:,1))-5 && ...
% %             xy(1,1)>5 && xy(1,1)<=length(handles.data.mrData(1,1,:))-5 % I do not know for sure whether X and Y are correct here
% %             % if already an electrode in that position remove
% %             if handles.data.elecMap(round(get(handles.sl_Z,'Value')),round(xy(1,2)),round(xy(1,1)))>0;
% %                 removeElec(hObject, eventdata, handles,xy);
% %             else
% %                 drawElec(hObject, eventdata, handles,xy);
% %             end
% %         end
% %     else
% %         %set(handles.display_feedback,'String','mis');
% %     end
% % end

% --- Function that draws electrode
function drawElec(hObject, eventdata, handles,xy) 
% set stuff
% radiusElect=str2double(get(handles.max_size,'String')); %5 in voxels, not mm
% radiusElect=1;
% 
% %     if handles.View==1
%         current.XYZ=[round(xy(1,2)) round(xy(1,1)) round(get(handles.ImgData.S,'Value'))];
% %     elseif handles.View==2
% %         current.XYZ=[round(xy(1,2)) round(get(handles.sl_Z,'Value')) round(xy(1,1))];
% %     elseif handles.View==3
% %         current.XYZ=[round(get(handles.sl_Z,'Value')) round(xy(1,2)) round(xy(1,1))];
% %     end
%     
% 
%     % make an empty matrix for electrode locations
%     data.elecMap=zeros(size(handles.data.mrData));
%     % put ones on a box of radiusElect aournd where there was clicked
%     data.elecMap(...
%         current.XYZ(1)-radiusElect:current.XYZ(1)+radiusElect,...
%         current.XYZ(2)-radiusElect:current.XYZ(2)+radiusElect,...
%         current.XYZ(3)-radiusElect:current.XYZ(3)+radiusElect)=1;
%     % put a 2 at the point were there was clicked
%     data.elecMap(...
%         current.XYZ(1),current.XYZ(2),current.XYZ(3))=2;
%     
%     % set data
%     handles.data.elecMap(handles.data.elecMap==1|data.elecMap==1)=1;
%     handles.data.elecMap(data.elecMap==2)=2;
% %    set(handles.display_feedback,'String','found electrode');
%     guidata(hObject,handles);
%     UpdateAxes(hObject,eventdata,handles);
%     
% function removeElec(hObject, eventdata, handles,xy)
% %set(handles.display_feedback,'String','removing electrode');
% % function that removes electrode
% if handles.View==1
% current.XYZ=[round(xy(1,2)) round(xy(1,1)) round(get(handles.sl_Z,'Value'))];
% elseif handles.View==2
% current.XYZ=[round(xy(1,2)) round(get(handles.sl_Z,'Value')) round(xy(1,1))];
% elseif handles.View==3
% current.XYZ=[round(get(handles.sl_Z,'Value')) round(xy(1,2)) round(xy(1,1))];
% end
% 
% data.elecMap=handles.data.elecMap;
% data.elecMap(current.XYZ(1),current.XYZ(2),current.XYZ(3))=3;
% data.elecMapV=data.elecMap(:);
% 
% %label blobs
% data.elecMap_bwlabel=bwlabeln(data.elecMap);
% data.elecMapV_bwlabel=data.elecMap_bwlabel(:);
% 
% % get label where the mouse clicked:
% data.clickpointbool=(data.elecMapV==3);
% data.clickpointbwlabel=data.elecMapV_bwlabel(data.clickpointbool);
% data.correctlabelbool=(data.elecMapV_bwlabel==data.clickpointbwlabel);
% 
% % find blob with correct label
% data.current_electV=data.correctlabelbool';
% data.elecMapV(data.current_electV)=0;
% data.elecMap(:)=data.elecMapV;
% 
% handles.data.elecMap=data.elecMap;
% %set(handles.display_feedback,'String','removed electrode');
% guidata(hObject,handles);
% UpdateAxes(hObject,eventdata,handles);
disp('Sorry not ready yet, come back soon!!')
