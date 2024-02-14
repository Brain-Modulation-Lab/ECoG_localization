function varargout = DBS_Elec_Localizer(varargin)
%SDK012215

% Edit the above text to modify the response to help DBS_Elec_Localizer

% Last Modified by GUIDE v2.5 04-Jan-2024 15:32:35

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DBS_Elec_Localizer_OpeningFcn, ...
                   'gui_OutputFcn',  @DBS_Elec_Localizer_OutputFcn, ...
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


% --- Executes just before DBS_Elec_Localizer is made visible.
function DBS_Elec_Localizer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DBS_Elec_Localizer (see VARARGIN)

% Choose default command line output for DBS_Elec_Localizer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes DBS_Elec_Localizer wait for user response (see UIRESUME)
% uiwait(handles.figure1);
axes(handles.ax1);
pos=get(gca,'Position');
set(handles.ax1,'Color','none')
hold on

axes(handles.ax2);
set(gca,'Position',pos);
axis off;
axes(handles.ax3);
set(gca,'Position',pos);
axis off;
axes(handles.ax2);

set(handles.ax1,'projection','perspective')%CHANGES FROM ORTHOGRAPHIC TO PERSPECTIVE PROJECTION


match = regexp(dbs_getrecentptfolder, '(?<=sub-)[a-zA-Z0-9-]+' ,'match');
assert(length(match)>=1); match = match{1}; 
handles.SubjectID = match; 

handles.Cortex.Hp=gobjects(1);
handles.Fl.HI=gobjects(1);
handles.El.He=gobjects(1);
handles.CortLocalizer.LH1=[];
handles.CortLocalizer.LH2=[];
handles.CortLocalizer.CEH=[];
handles.FidLocalizer.LH1=[];
handles.FidLocalizer.LH2=[];
handles.FidLocalizer.SEH=[];
handles.FidLocalizer.MEH=[];
handles.misc.LockAx=0;
handles.Skull=struct;
handles.SurfFlip = 0;


handles.CortLocalizer.Status1=0;
handles.CortLocalizer.Status2=0;
handles.CortLocalizer.ElecLoc={};
handles.CortLocalizer.ElecLoc2={};
handles.FidLocalizer.Status1=0;
handles.FidLocalizer.Status2=0;
handles.FidLocalizer.ElecLoc={};
handles.FidLocalizer.ElecLoc2={};
handles.CortLocalizer.numelec=0;
handles.Nudge.mode=0;
handles.Nudge.ScaleMax=20;

handles.FluoroLocalizer.Status1=0;
handles.FluoroLocalizer.emitter2nose = 650; % default distance from  mm from nose to screen
handles.FluoroLocalizer.emitter2screen = 780; % fixed distance from  mm from emitter to screen
handles.FluoroLocalizer.emitterSide = 'left'; % left or right side? 
handles.FluoroLocalizer.width = 244; % fixed width of fluoro in mm
handles.FluoroLocalizer.fov_radius = 114; % fixed radius of fluoro in mm
handles.FluoroLocalizer.Landmarks = table();
handles.FluoroLocalizer.LandmarksFileName = ['./sub-' handles.SubjectID '_fluoro-landmarks.tsv'];
handles.FluoroLocalizer.LandmarksNames = {'pintip_front', 'pintip_occ', 'dbslead_bottom', 'dbslead_top'};
% handles.FluoroLocalizer.OptimizationSolutionFileName = ['./sub-' handles.SubjectID '_optimized-camera-solution.tsv'];
handles.FluoroLocalizer.OptimizationSolutionMATFileName = ['./sub-' handles.SubjectID '_eleclocalizer-camera-optimization-solution.mat'];
handles.FluoroLocalizer.OptimizationSolutionElectrodesFileName = ['./sub-' handles.SubjectID '_eleclocalizer-all-landmarks-final.tsv'];


set(handles.ed_Cmd,'KeyReleaseFcn',@CmdKeyRelease)

guidata(hObject,handles)


% --- Outputs from this function are returned to the command line.
function varargout = DBS_Elec_Localizer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function mn_File_Callback(hObject, eventdata, handles)
% hObject    handle to mn_File (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.Camera.cth=cameratoolbar;
handles.Camera.cp=get(handles.ax1,'CameraPosition');
handles.Camera.ct=get(handles.ax1,'CameraTarget');
handles.Origin.trans=[0,0,0];



% --------------------------------------------------------------------
function mn_LoadElecLocs_Callback(hObject, eventdata, handles)
% hObject    handle to mn_LoadElecLocs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[El_fname,El_pname,~] = uigetfile('*.mat');
load(fullfile(El_pname,El_fname));
handles.El.elecmatrix=elecmatrix;
delete(handles.El.He);
axes(handles.ax1);

elecmatrix1 = elecmatrix(elecmatrix(:,1)<0,:);
elecmatrix2 = elecmatrix(elecmatrix(:,1)>0,:);

% for i=1:size(elecmatrix,1)
%     handles.El.He(i)=plot3(elecmatrix(i,1),elecmatrix(i,2),elecmatrix(i,3),'.','MarkerSize',25);
% end
if ~isempty(elecmatrix1)
    [~,si] = sort(elecmatrix1(:,3));
    elecmatrix1 = elecmatrix1(si,:);
%     for i=1:3
%         elecmatrix1(:,i) = smooth(elecmatrix1(:,i),3) ;
%     end
    handles.El.He(1) = plot3(elecmatrix1(:,1),elecmatrix1(:,2),elecmatrix1(:,3),'linewidth',3,'color','r');
end
if ~isempty(elecmatrix2)
     [~,si] = sort(elecmatrix2(:,3));
    elecmatrix2 = elecmatrix2(si,:);
%     for i=1:3
%         elecmatrix2(:,i) = smooth(elecmatrix2(:,i),3) ;
%     end
    handles.El.He(2) = plot3(elecmatrix2(:,1),elecmatrix2(:,2),elecmatrix2(:,3),'linewidth',3,'color','b');
end
guidata(hObject,handles)

% --------------------------------------------------------------------
function handles = mn_LoadFluoro_Callback(hObject, eventdata, handles)
% hObject    handle to mn_LoadFluoro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[Fl_fname,Fl_pname] = dbs_getfluoroname;
%[folders, files] = dbs_subdir(ptroot);
%strip_dir = folders{~cellfun(@isempty,strfind(folders,'Strip'))};
fsroot = dbs_getrecentfsfolder;
ptroot = dbs_getrecentptfolder;
try
    %handles.Fl.I=imread(fullfile(Fl_pname,Fl_fname));
    fnames = strtrim(string(ls(ptroot))); 
    fname = fnames(contains(fnames, 'fluoro', 'IgnoreCase', true)); 
    handles.Fl.I=imread(fullfile(ptroot,fname));
catch
    [Fl_fname,Fl_pname,ext] = dbs_uigetfile(ptroot,'Choose Fluoro');%{'*.tif'});
    handles.Fl.I=imread(fullfile(Fl_pname,[Fl_fname,ext]));
end

if size(handles.Fl.I,3)>3
    handles.Fl.I=handles.Fl.I(:,:,1:3);
end


% [Fl_fname,Fl_pname] = dbs_getfluoroname;
% try
%     handles.Fl.I=imread(fullfile(Fl_pname,Fl_fname));
% catch
%     [Fl_fname,Fl_pname,~] = uigetfile({'*.tif'});
%     handles.Fl.I=imread(fullfile(Fl_pname,Fl_fname));
% end
% axes(handles.ax2)
% if size(handles.Fl.I,3)>3
%     handles.Fl.I=handles.Fl.I(:,:,1:3);
% end

axes(handles.ax2)
handles.Fl.HI=imshow(handles.Fl.I);

set(handles.Fl.HI, 'visible', 0); 
set(handles.rb_fluoro_visible,'Value', 1); 
rb_fluoro_visible_Callback(hObject, eventdata, handles);



%% Load fluoro landmarks
prompt = {'Enter distance from fluoro emitter to nasion:'};
dlgtitle = '';
fieldsize = [1 45];
definput = {'670'};
answer = inputdlg(prompt,dlgtitle,fieldsize,definput); 
handles.FluoroLocalizer.emitter2nose = str2num(answer{1}); 

answer = questdlg('Fluoro emitter/camera side?', '', ...
    'left', "right", 'left');
handles.FluoroLocalizer.emitterSide = string(answer); 

flip = double(handles.FluoroLocalizer.emitterSide=="left")*2-1; 
handles.Camera.cp = [-1*flip*handles.FluoroLocalizer.emitter2nose, 0, 0];
handles.Camera.ct = [0, 0, 0];
handles.Camera.cva = 20;
handles.Camera.uv = [0, 0, 1]; 
bt_ResetCamera_Callback(hObject, eventdata, handles);


opts = detectImportOptions(handles.FluoroLocalizer.LandmarksFileName, 'FileType', 'text');
opts = setvartype(opts, {'hemi', 'name'}, "string");
lm = readtable(handles.FluoroLocalizer.LandmarksFileName, opts);
lm.pos = [zeros(height(lm), 1) -lm.pos_1 lm.pos_2 ];
lm.coordframe = repmat("fluoro", [height(lm), 1]);
lm(:, {'pos_1', 'pos_2'}) = [];

% % if fluoro was marked at a different distance, 
% lm = handles.FluoroLocalizer.Landmarks;
% landmarks_fluoro_depth = mean(lm.pos(lm.coordframe=="fluoro", 1)); 
% fluoro_depth = handles.FluoroLocalizer.emitter2screen - handles.FluoroLocalizer.emitter2nose; 
% % if abs(fluoro_depth - landmarks_fluoro_depth) > 5 
% %     warning("Depth of fluoro landmarks don't match distance. Removing fluoro landmarks from landmarks table.")
% %     lm(lm.coordframe=="fluoro", :) = []
% % end

% project landmarks according to offset distance and fluoro width
scale = eye(3)*(handles.FluoroLocalizer.width/2);
% flip = double(handles.FluoroLocalizer.emitterSide=="left")*2-1; 
translate = [1; 0; 0]*(handles.FluoroLocalizer.emitter2screen - handles.FluoroLocalizer.emitter2nose);
R = [scale translate; 
    0 0 0  1        ];
pos = [lm.pos'; ones([1 width(lm.pos')])];
pos = (R*pos)';
lm.pos = pos(:, 1:3);

handles.FluoroLocalizer.Landmarks = [handles.FluoroLocalizer.Landmarks; lm];


guidata(hObject, handles); 


% --- Executes on button press in bt_CortStack.
function bt_CortStack_Callback(hObject, eventdata, handles)
% hObject    handle to bt_CortStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bt_CortStack
% set(handles.bt_FluoroStack,'value',0)
% set(handles.bt_SkullStack,'value',0)
axes(handles.ax1);

% --- Executes on button press in bt_FluoroStack.
function bt_FluoroStack_Callback(hObject, eventdata, handles)
% hObject    handle to bt_FluoroStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hint: get(hObject,'Value') returns toggle state of bt_FluoroStack
% set(handles.bt_CortStack,'value',0)
% set(handles.bt_SkullStack,'value',0)
% set(handles.bt_FluoroStack,'value',0)
% set(handles.bt_CortStack,'value',0)
set(handles.bt_FluoroStack,'value',0)
set(handles.bt_CortStack,'value',0)
axes(handles.ax2);

guidata(hObject,handles)



% --- Executes on button press in bt_SkullStack.
function bt_SkullStack_Callback(hObject, eventdata, handles)
% hObject    handle to bt_SkullStack (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of bt_SkullStack
set(handles.bt_FluoroStack,'value',0)
set(handles.bt_CortStack,'value',0)
axes(handles.ax3);


% --- Executes on slider movement.
function sl_Cort_Callback(hObject, eventdata, handles)
% hObject    handle to sl_Cort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.Cortex.Hp,'FaceAlpha',get(hObject,'Value'))


% --- Executes during object creation, after setting all properties.
function sl_Cort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_Cort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function sl_Fluoro_Callback(hObject, eventdata, handles)
% hObject    handle to sl_Fluoro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.Fl.HI,'AlphaData',get(hObject,'Value'))


% --- Executes during object creation, after setting all properties.
function sl_Fluoro_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_Fluoro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sl_Skull_Callback(hObject, eventdata, handles)
% hObject    handle to sl_Skull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.Skull.Hp,'FaceAlpha', get(hObject,'Value'))

% --- Executes during object creation, after setting all properties.
function sl_Skull_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_Skull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function DrawElecLine(hObject, eventdata, handles)
handles=guidata(hObject);

% Witek debug
fprintf('in DrawElecLine...');

if strcmp(get(gcf,'SelectionType'), 'normal') %add contact on left click
    pos=get(gca,'CurrentPoint');
    handles.CortLocalizer.ElecLoc{end+1}=pos;
    dthr = 1;
    vert=handles.Cortex.hull; %Use cortical hull
    
    Rside = get(handles.rb_elecSideR,'value');
    Lside = get(handles.rb_elecSideL,'value');
    
    if Rside == 1
        vert = vert(vert(:,1)>0,:); %CHOOSE LEFT OR RIGHT BRAIN
    elseif Lside == 1
        vert = vert(vert(:,1)<0,:);
    end
    
    lm = handles.FluoroLocalizer.Landmarks;
    lm = lm(startsWith(lm.name, "ecog"), :);
    lm = sortrows(lm, 'name');
    cp = get(handles.ax1, 'CameraPosition'); 

    for ielec = 1:height(lm)
        pos = lm.sol_pos(ielec, :);
        pos1 = [linspace(pos(1), cp(1), 5000);
                linspace(pos(2), cp(2), 5000);
                linspace(pos(3), cp(3), 5000);]';

        %     pos1=[linspace(pos(1,1),pos(2,1),5000)', linspace(pos(1,2),pos(2,2),5000)', linspace(pos(1,3),pos(2,3),5000)'];
        dst = single(pdist2(vert,pos1));
        mi = (min(dst,[],2)<dthr); %Find vertices within dthr of our line
        vert2 = vert(mi,:);
        [~,mi2] = max(abs(vert2(:,1))); % find the more lateral intersection
        epos = vert2(mi2,:);


        % %---Orthogonal Distance Method
        % campos=get(gca,'CameraPosition');
        % mag=@(A) sqrt( sum(A.^2) );
        % vert0=vert;
        % tmp=bsxfun(@minus,pos,campos);
        % for i=1:2
        %     d(i)=mag(tmp(i,:));
        % end
        % [~,i1]=min(d);
        % [~,i2]=max(d);
        % O=pos(i1,:);
        % B=pos(i2,:)-O;%Move B so O is origin
        %
        % vert=bsxfun(@minus,vert,O);%Move verts so O is origin
        %
        % fc = @(A,B) sqrt( mag(A)^2 - ( dot(A,B) / mag(B) )^2   );
        %
        % A=num2cell(vert,2);
        % B1=num2cell( repmat(B,size(A,1),1), 2);
        %
        % c=cellfun(fc,A,B1);
        %
        % odthr=2; %MUST CHANGE THRESHOLD TO ACTUAL DISTANCE RATHER THAN PIXELS;
        % ci=find(c<odthr);
        %
        %
        % tmp=num2cell(vert(ci,:),2);
        % tmp=cellfun(mag,tmp);
        % [~,ci2]=min(tmp);
        %
        % tmp=vert0(ci(ci2),:);
        % for i=1:3
        %    tmp1(:,i)=linspace(pos(1,i),pos(2,i),1000);
        % end
        %
        % [~,mi]=min(cellfun(mag,num2cell(bsxfun(@minus,tmp1,tmp),2)));
        % epos=tmp1(mi,:);

        colormap jet
        col=colormap;
        colormap gray
        coldr=round(length(col)/8);
        col=col(3:17:end-3,:);

        handles.CortLocalizer.numelec = handles.CortLocalizer.numelec + 1;
        handles.CortLocalizer.ElecLocCort{handles.CortLocalizer.numelec}=epos;
        coli=repmat(1:length(col),1,10);
        % Witek edit: draw all electrodes in blue instead of cycling through
        % colors
        %     handles.CortLocalizer.CEH(handles.CortLocalizer.numelec)=plot3(epos(1),epos(2),epos(3),'.','MarkerSize',20,'Color',col( coli( handles.CortLocalizer.numelec ),:) );
        handles.CortLocalizer.CEH(handles.CortLocalizer.numelec)=plot3(handles.ax1, epos(1),epos(2),epos(3),'.','MarkerSize',20,'Color','b' );
        drawnow

        nlm = [];
        nlm.coordframe = "recon";
        nlm.sol_pos = epos;
        nlm.pos = [nan nan nan];
        nlm.hemi = lm.hemi(ielec);
        nlm.name = lm.name(ielec);
        handles.FluoroLocalizer.Landmarks = [handles.FluoroLocalizer.Landmarks; struct2table(nlm)];
        % Witek debug
        fprintf('... contact %d drawn.\n', length(handles.CortLocalizer.ElecLoc))

    end
else %delete last contact on right click
    if ~isempty(handles.CortLocalizer.ElecLoc)
        delete(handles.CortLocalizer.CEH(handles.CortLocalizer.numelec));
        handles.CortLocalizer.CEH(handles.CortLocalizer.numelec) = [];
        % Witek debug
        fprintf('... contact %d deleted.\n', length(handles.CortLocalizer.ElecLoc))
        handles.CortLocalizer.ElecLoc = handles.CortLocalizer.ElecLoc(1:end-1);
        handles.CortLocalizer.ElecLocCort = handles.CortLocalizer.ElecLocCort(1:end-1);
        handles.CortLocalizer.numelec = handles.CortLocalizer.numelec - 1;
    else
        % Witek debug
        fprintf('... no contacts to delete.\n')
    end
end

guidata(hObject,handles)
fprintf('Finished projecting electrodes \n')


% --- Executes on button press in bt_ResetCamera.
function bt_ResetCamera_Callback(hObject, eventdata, handles)
% hObject    handle to bt_ResetCamera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.ax1,'CameraPosition',handles.Camera.cp,...
'CameraTarget',handles.Camera.ct,...
'CameraViewAngle',handles.Camera.cva,...
'CameraUpVector',handles.Camera.uv)
guidata(hObject,handles)


% --- Executes on button press in bt_SaveCamera.
function handles = bt_SaveCamera_Callback(hObject, eventdata, handles)
% hObject    handle to bt_SaveCamera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.Camera.cp=get(handles.ax1,'CameraPosition');
handles.Camera.ct=get(handles.ax1,'CameraTarget');
handles.Camera.cva=get(handles.ax1,'CameraViewAngle');
handles.Camera.uv=get(handles.ax1,'CameraUpVector');

cam = handles.Camera; 
guidata(hObject,handles)


function [cp, ct, cuv, cva] = get_current_camera_props(hObject, eventdata, handles)
cp=get(handles.ax1,'CameraPosition');
ct=get(handles.ax1,'CameraTarget');
cuv=get(handles.ax1,'CameraUpVector');
cva=get(handles.ax1,'CameraViewAngle');

% --- Executes on button press in bt_CoRegCort1.
function bt_CoRegCort1_Callback(hObject, eventdata, handles)
% hObject    handle to bt_CoRegCort1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

DrawElecLine(hObject, eventdata, handles)

% handles.CortLocalizer.ax1ch=[handles.ax1;get(handles.ax1,'children')];
% if handles.CortLocalizer.Status1==0
%     axes(handles.ax1)
%     set(handles.CortLocalizer.ax1ch,'ButtonDownFcn',@DrawElecLine)
%     handles.CortLocalizer.Status1=1;
%     set(handles.tx_LocMode,'string','Cort Localizer Mode ON','ForegroundColor',[1,0,0])
%     set(handles.figure1,'pointer','crosshair')
% elseif handles.CortLocalizer.Status1==1
%     set(handles.CortLocalizer.ax1ch,'ButtonDownFcn',[]);
%     handles.CortLocalizer.Status1=0;
%     set(handles.tx_LocMode,'string','Cort Localizer Mode OFF','ForegroundColor',[0,0,0])
%     set(handles.figure1,'pointer','arrow')
% end
% guidata(hObject,handles)


% --- Executes on button press in bt_CoRegCort2.
function bt_CoRegCort2_Callback(hObject, eventdata, handles)
% hObject    handle to bt_CoRegCort2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.CortLocalizer.ax1ch=[handles.ax1;get(handles.ax1,'children')];
if handles.CortLocalizer.Status2==0
    axes(handles.ax1)
    set(handles.CortLocalizer.ax1ch,'ButtonDownFcn',@SelectElecLine)
    handles.CortLocalizer.Status2=1;
    set(handles.tx_LocMode,'string','Cort Localizer Mode ON','ForegroundColor',[1,0,0])
elseif handles.CortLocalizer.Status2==1
    set(handles.CortLocalizer.ax1ch,'ButtonDownFcn',[]);
    handles.CortLocalizer.Status2=0;
    set(handles.tx_LocMode,'string','Cort Localizer Mode OFF','ForegroundColor',[0,0,0])
end
guidata(hObject,handles)

function SelectElecLine(hObject, eventdata, handles)
handles=guidata(hObject);
handles.CortLocalizer.co=gco;
handles.CortLocalizer.LH1_i=find(ismember(handles.CortLocalizer.LH1,handles.CortLocalizer.co));
if size(handles.CortLocalizer.LH1_i,2)==1
    set(handles.CortLocalizer.ax1ch,'ButtonDownFcn',@IntersectElecLine)
    guidata(hObject,handles)
else
    wh=warndlg(sprintf('Try again\nClick Electrode Line'));
    uiwait(wh);
end

function IntersectElecLine(hObject, eventdata, handles)
handles=guidata(hObject);

a=1;

co=handles.CortLocalizer.co;
col=get(co,'Color');
LH1_i=handles.CortLocalizer.LH1_i;
pos=get(gca,'CurrentPoint');
try
    if ~isempty(handles.CortLocalizer.ElecLoc2{LH1_i})
        delete(handles.CortLocalizer.LH2(LH1_i))
    end
catch
end
handles.CortLocalizer.ElecLoc2{LH1_i}=pos;
handles.CortLocalizer.LH2(LH1_i)=plot3(pos(:,1),pos(:,2),pos(:,3),'LineWidth',6,'Marker','o','MarkerSize',10,'Color',col);
guidata(hObject,handles)

A=handles.CortLocalizer.ElecLoc{LH1_i};
B=handles.CortLocalizer.ElecLoc2{LH1_i};
A1=A(1,:);A2=A(2,:);
B1=B(1,:);B2=B(2,:);

nA = dot(cross(B2-B1,A1-B1),cross(A2-A1,B2-B1));
nB = dot(cross(A2-A1,A1-B1),cross(A2-A1,B2-B1));
d = dot(cross(A2-A1,B2-B1),cross(A2-A1,B2-B1));
A0 = A1 + (nA/d)*(A2-A1);
B0 = B1 + (nB/d)*(B2-B1);

handles.CortLocalizer.ElecLocCort{LH1_i}=A0;
handles.CortLocalizer.CEH(LH1_i)=plot3(A0(1),A0(2),A0(3),'.','MarkerSize',50,'Color',col);

handles.CortLocalizer.ax1ch=[handles.ax1;get(handles.ax1,'children')];
set(handles.CortLocalizer.ax1ch,'ButtonDownFcn',@SelectElecLine)

guidata(hObject,handles)

% --- Executes on button press in rb_Depths.
function rb_Depths_Callback(hObject, eventdata, handles)
% hObject    handle to rb_Depths (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_Depths
vischoice={'off','on'};
vis=get(hObject,'value');
for i = 1:length(handles.El.He)
    if isgraphics(handles.El.He(i))
        set(handles.El.He(i),'Visible',vischoice{vis+1})
    end
end

% --- Executes on button press in rb_CortL.
function rb_CortL_Callback(hObject, eventdata, handles)
% hObject    handle to rb_CortL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLABmanual
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_CortL
vischoice={'off','on'};
vis=get(hObject,'value');
if ~isempty(handles.CortLocalizer.LH1), set(handles.CortLocalizer.LH1,'Visible',vischoice{vis+1}); end
if ~isempty(handles.CortLocalizer.LH2), set(handles.CortLocalizer.LH2,'Visible',vischoice{vis+1}); end


% --- Executes on button press in rb_CortC.
function rb_CortC_Callback(hObject, eventdata, handles)
% hObject    handle to rb_CortC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_CortC
vischoice={'off','on'};
vis=get(hObject,'value');
if ~isempty(handles.CortLocalizer.CEH), set(handles.CortLocalizer.CEH,'Visible',vischoice{vis+1}); end


% --------------------------------------------------------------------
function mn_Coregistration_Callback(hObject, eventdata, handles)
% hObject    handle to mn_Coregistration (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mn_Camera_Callback(hObject, eventdata, handles)
% hObject    handle to mn_Camera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mn_SaveCamera_Callback(hObject, eventdata, handles)
% hObject    handle to mn_SaveCamera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname, pname, ~] = uiputfile('*.mat', 'Save Camera Position');
CameraPosition.cp=get(handles.ax1,'CameraPosition');
CameraPosition.ct=get(handles.ax1,'CameraTarget');
CameraPosition.cva=get(handles.ax1,'CameraViewAngle');
CameraPosition.uv=get(handles.ax1,'CameraUpVector');
save(fullfile(pname,fname),'CameraPosition')

% ax1=handles.ax1;
% Hp=handles.Cortex.Hp
% CEH=handles.CortLocalizer.CEH;
% save('TEMP','ax1','Hp','CEH')


% --------------------------------------------------------------------
function handles = mn_LoadCamera_Callback(hObject, eventdata, handles)
% hObject    handle to mn_LoadCamera (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fsroot = dbs_getrecentfsfolder; 
ptroot = dbs_getrecentptfolder;
try
    fnames = strtrim(string(ls(fullfile(fsroot, 'Electrode_Locations', 'CamPos-*.mat'))));
    fname = fnames(end);
    load(fullfile(fsroot, 'Electrode_Locations', fname));
catch
    [fname, pname, ~] = uigetfile('*.mat', 'Choose Camera Position File');
    load(fullfile(pname,fname))
end

handles.Camera=CameraPosition;
set(handles.ax1,'CameraPosition',handles.Camera.cp,...
'CameraTarget',handles.Camera.ct,...
'CameraViewAngle',handles.Camera.cva,...
'CameraUpVector',handles.Camera.uv)
% set(handles.tx_XLoc,'string',num2str(CameraPosition.cp(1)));

guidata(hObject,handles)

% --------------------------------------------------------------------
function mn_ClearCort_Callback(hObject, eventdata, handles)
% hObject    handle to mn_ClearCort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.CortLocalizer.CEH)
    delete(handles.CortLocalizer.CEH);
    handles.CortLocalizer.CEH=[];
    handles.CortLocalizer.ElecLoc={};
    handles.CortLocalizer.ElecLocCort={};
    handles.CortLocalizer.numelec=0;
end
guidata(hObject,handles)


% --------------------------------------------------------------------
function mn_Export_Callback(hObject, eventdata, handles)
% hObject    handle to mn_Export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function mn_ExpBrain_Callback(hObject, eventdata, handles)
% hObject    handle to mn_ExpBrain (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=1;
Hp=handles.Cortex.Hp;
CEH1=handles.CortLocalizer.CEH;
CEH=gobjects(1);
for i=1:length(CEH1)
    CEH(i)=CEH1(i);
end
[fname, pname] = uiputfile('*.mat', 'Export Brain to:');
save(fullfile(pname,fname),'CEH','Hp')


% --------------------------------------------------------------------
function mn_ExpElCoord_Callback(hObject, eventdata, handles)
% hObject    handle to mn_ExpElCoord (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=1;
CortElecLoc=handles.CortLocalizer.ElecLocCort;

% if handles.SurfFlip == 1
%     cfcn = @(A) [A(1)*-1, A(2), A(3)];
%     CortElecLoc = cellfun(cfcn, CortElecLoc,'uniformoutput',false);
% end
% [fname, pname] = uiputfile('*.mat', 'Export Electrode Locs to:');
% save(fullfile(pname,fname),'CortElecLoc')

mn_SaveCamera_Callback(hObject, eventdata, handles);

% ----- write landmarks
lm_tbl_fname = handles.FluoroLocalizer.OptimizationSolutionElectrodesFileName; 
lm = handles.FluoroLocalizer.Landmarks; 
% lm = lm(lm.coordframe=="recon", :);
writetable(lm, lm_tbl_fname, 'Delimiter', '\t', 'FileType', 'text'); 
f = msgbox("Saved successfully!");



% --- Executes on button press in bt_MarkFid1.
function bt_MarkFid1_Callback(hObject, eventdata, handles)
% hObject    handle to bt_MarkFid1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FidLocalizer.gcach=[gca;get(gca,'children')];
if handles.FidLocalizer.Status1==0
    set(handles.FidLocalizer.gcach,'ButtonDownFcn',@DrawFidLine)
    handles.FidLocalizer.Status1=1;
    set(handles.tx_LocMode,'string','Cort Localizer Mode ON','ForegroundColor',[1,0,0])
elseif handles.FidLocalizer.Status1==1
    set(handles.FidLocalizer.gcach,'ButtonDownFcn',[]);
    handles.FidLocalizer.Status1=0;
    set(handles.tx_LocMode,'string','Cort Localizer Mode OFF','ForegroundColor',[0,0,0])
end

guidata(hObject,handles)

% --- Executes on button press in bt_MarkFid2.
function bt_MarkFid2_Callback(hObject, eventdata, handles)
% hObject    handle to bt_MarkFid2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.FidLocalizer.gcach=[gca;get(gca,'children')];
if handles.FidLocalizer.Status2==0
    set(handles.FidLocalizer.gcach,'ButtonDownFcn',@SelectFidLine)
    handles.FidLocalizer.Status2=1;
    set(handles.tx_LocMode,'string','Cort Localizer Mode ON','ForegroundColor',[1,0,0])
elseif handles.FidLocalizer.Status2==1
    set(handles.FidLocalizer.gcach,'ButtonDownFcn',[]);
    handles.FidLocalizer.Status2=0;
    set(handles.tx_LocMode,'string','Cort Localizer Mode OFF','ForegroundColor',[0,0,0])
end
guidata(hObject,handles)

function DrawFidLine(hObject, eventdata, handles)
handles=guidata(hObject); hold on
pos=get(gca,'CurrentPoint');
handles.FidLocalizer.ElecLoc{end+1}=pos;
handles.FidLocalizer.LH1(end+1)=plot3(pos(:,1),pos(:,2),pos(:,3),'LineWidth',6,'Marker','o');
guidata(hObject,handles)


function SelectFidLine(hObject, eventdata, handles)
handles=guidata(hObject);
handles.FidLocalizer.co=gco;
handles.FidLocalizer.LH1_i=find(ismember(handles.FidLocalizer.LH1,handles.FidLocalizer.co));
if size(handles.FidLocalizer.LH1_i,2)==1
    set(handles.FidLocalizer.gcach,'ButtonDownFcn',@IntersectFidLine)
    guidata(hObject,handles)
else
    wh=warndlg(sprintf('Try again\nClick Electrode Line'));
    uiwait(wh);
end


function IntersectFidLine(hObject, eventdata, handles)
handles=guidata(hObject);
co=handles.FidLocalizer.co;
col=get(co,'Color');
LH1_i=handles.FidLocalizer.LH1_i;
pos=get(gca,'CurrentPoint');
try
    if ~isempty(handles.FidLocalizer.ElecLoc2{LH1_i})
        delete(handles.FidLocalizer.LH2(LH1_i))
    end
catch
end
handles.FidLocalizer.ElecLoc2{LH1_i}=pos;
handles.FidLocalizer.LH2(LH1_i)=plot3(pos(:,1),pos(:,2),pos(:,3),'LineWidth',6,'Marker','o','MarkerSize',10,'Color',col);
guidata(hObject,handles)

A=handles.FidLocalizer.ElecLoc{LH1_i};
B=handles.FidLocalizer.ElecLoc2{LH1_i};
A1=A(1,:);A2=A(2,:);
B1=B(1,:);B2=B(2,:);

nA = dot(cross(B2-B1,A1-B1),cross(A2-A1,B2-B1));
nB = dot(cross(A2-A1,A1-B1),cross(A2-A1,B2-B1));
d = dot(cross(A2-A1,B2-B1),cross(A2-A1,B2-B1));
A0 = A1 + (nA/d)*(A2-A1);
B0 = B1 + (nB/d)*(B2-B1);

handles.FidLocalizer.SkullFid{LH1_i}=A0';
handles.FidLocalizer.SEH(LH1_i)=plot3(A0(1),A0(2),A0(3),'.','MarkerSize',50,'Color',col);

handles.FidLocalizer.ax1ch=[handles.ax1;get(handles.ax1,'children')];
set(handles.FidLocalizer.ax1ch,'ButtonDownFcn',@SelectElecLine)

set(co,'visible','off')
set(handles.FidLocalizer.LH2(LH1_i),'visible','off')

guidata(hObject,handles)

% --------------------------------------------------------------------
function mn_LoadCortRecon_Callback(hObject, eventdata, handles)
% hObject    handle to mn_LoadCortRecon (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[CR_fname,CR_pname,~] = uigetfile('*.mat');
load(fullfile(CR_pname,CR_fname))
%handles.Cortex.cortex=cortex;  NEEDED FOR MANUAL COREGISTRATION, REMOVED
%FOR NOW
axes(handles.ax1)
delete(handles.Cortex.Hp)
handles.Cortex.Hp = patch('vertices',cortex.vert,'faces',cortex.tri(:,[1 3 2]),...
    'facecolor',[.65 .65 .65],'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .25);
camlight('headlight','infinite');
axis equal
set(gca,'DataAspectRatioMode','manual','PlotBoxAspectRatioMode','manual');
set(gca,'camerapositionmode','manual','cameratargetmode','manual','cameraupvectormode','manual','cameraviewanglemode','manual')
set(gca,'projection','perspective')

% els = cell2mat(CortElecLoc');
% hold on; plot3(els(:,1),els(:,2),els(:,3),'.','color','r','markersize',10)
% handles.Cortex.Hp.Facealpha = 1;

%handles.Cortex.OSkull=patch('vertices',skull.vert,'faces',skull.tri(:,[1 3 2]),'facecolor',[.65 .65 .65],'edgecolor','none','FaceAlpha',0.2);
%handles.Cortex.ISkull=patch('vertices',skull_inner.vert,'faces',skull_inner.tri(:,[1 3 2]),'facecolor',[.65 .65 .65],'edgecolor','none','FaceAlpha',0.2);

guidata(hObject,handles)


% --------------------------------------------------------------------
function mn_LoadCortHull_Callback(hObject, eventdata, handles)
% hObject    handle to mn_LoadCortHull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname,~] = uigetfile('*.mat');
load(fullfile(pname,fname))
handles.Cortex.hull = mask_indices;
handles.Cortex.Hh = plot3(handles.Cortex.hull(:,1),handles.Cortex.hull(:,2),...
    handles.Cortex.hull(:,3),'.','color','m');

pause(0.5)

set(handles.Cortex.Hh,'visible','off')

guidata(hObject,handles)





% --------------------------------------------------------------------
function mn_LoadSkull_Callback(hObject, eventdata, handles)
% hObject    handle to mn_LoadSkull (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname,pname,~] = uigetfile('*.mat');
load(fullfile(pname,fname))

axes(handles.ax1);

% if exist('skull')
%     if isfield(skull,'faces'), skull.tri=skull.faces; end
% elseif exist('skin')
%     skull=skin; clear skin
%     if isfield(skull,'faces'), skull.tri=skull.faces; end
% end
% if isfield(handles,'Skull')
%     if isfield(handles.Skull,'Hp'), delete(handles.Skull.Hp), end
% end
% 
% s=size(skull.vert);
% if s(1)<s(2), skull.vert=skull.vert'; end
% s=size(skull.tri);
% if s(1)<s(2), skull.tri=skull.tri'; end
% 
% a=[-1,0,0;0,-1,0;0,0,1];
% skull.vert = [skull.vert*a];

handles.Skull.Hp = patch('vertices',skull.vert,'faces',skull.tri(:,[1 3 2]),'facecolor',[.65 .65 .65],'edgecolor','none',...
'facelighting', 'gouraud', 'specularstrength', .25);
set(handles.Skull.Hp,'parent',handles.ax1)

axis equal
camlight('headlight','infinite')
set(gca,'DataAspectRatioMode','manual','PlotBoxAspectRatioMode','manual');
set(gca,'camerapositionmode','manual','cameratargetmode','manual','cameraupvectormode','manual','cameraviewanglemode','manual')

guidata(hObject,handles)



function ed_Cmd_Callback(hObject, eventdata, handles)
% hObject    handle to ed_Cmd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ed_Cmd as text
%        str2double(get(hObject,'String')) returns contents of ed_Cmd as a double


% --- Executes during object creation, after setting all properties.
function ed_Cmd_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ed_Cmd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%-----Executes During Key Release on Command Edit Line
function CmdKeyRelease(hObject,eventdata,handles) 
if strcmp(eventdata.Key,'return')
    cmd=get(hObject,'string'); %cmd=cmd{1};
    if ~strcmp(cmd,'Enter Command')
        handles=guidata(hObject);
        eval(cmd);
        set(hObject,'string','Enter Command')
        guidata(hObject,handles)
    end
end


% --------------------------------------------------------------------
function mn_ExpSkullFid_Callback(hObject, eventdata, handles)
% hObject    handle to mn_ExpSkullFid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SkullFiducials=handles.FidLocalizer.SkullFid;
[fname, pname, ~] = uiputfile('*.mat', 'Save Skull Fiducials','Skull_Fiducials.mat');
save(fullfile(pname,fname),'SkullFiducials');


% --------------------------------------------------------------------
function mn_Load_Callback(hObject, eventdata, handles)
% hObject    handle to mn_Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mn_Load_MRIFid_Callback(hObject, eventdata, handles)
% hObject    handle to mn_Load_MRIFid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname,~] = uigetfile('*.mat');
load(fullfile(pname,fname))
%handles.FidLocalizer.MRIFid=MRIFiducials;
[~,I]=sort(fiducial_locations(:,1));%Sorts in ascending order
fiducial_locations=fiducial_locations(I,:);
for i=1:3
    handles.FidLocalizer.MRIFid{i}= fiducial_locations(i,:)';
    %handles.FidLocalizer.MRIFid{i}= fiducial_locations(i,:)';

end
axes(handles.ax1); hold on
for i=1:length(handles.FidLocalizer.MRIFid)
    A0=handles.FidLocalizer.MRIFid{i};
    handles.FidLocalizer.MEH(i)=plot3(A0(1),A0(2),A0(3),'.','MarkerSize',50,'Color','b');
end

guidata(hObject,handles)



% --------------------------------------------------------------------
function mn_Load_SkullFid_Callback(hObject, eventdata, handles)
% hObject    handle to mn_Load_SkullFid (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname,~] = uigetfile('*.mat');
load(fullfile(pname,fname))

if exist('SkullFiducials','var')
    handles.FidLocalizer.SkullFid=SkullFiducials;
    axes(handles.ax3); hold on
    for i=1:length(handles.FidLocalizer.SkullFid)
        A0=handles.FidLocalizer.SkullFid{i};
        handles.FidLocalizer.SEH(i)=plot3(A0(1),A0(2),A0(3),'.','MarkerSize',50,'Color','r');
    end
end
if exist('fiducial_locations','var')
        if isfield(handles.FidLocalizer,'SEH2');
            delete(handles.FidLocalizer.SEH2)
        end
        axes(handles.ax1);
        handles.ax1.Visible='off';
        nf = size(fiducial_locations,2);
        colormap jet
        cmap = colormap;
        c=cmap(1:floor(length(cmap)/nf):end,:);
        handles.FidLocalizer.SEH2 = scatter3(fiducial_locations(:,1),fiducial_locations(:,2),fiducial_locations(:,3),50,c,'+','linewidth',3);
end


guidata(hObject,handles)


% --- Executes on button press in bt_CoReg_MR_CT.
function bt_CoReg_MR_CT_Callback(hObject, eventdata, handles)
% hObject    handle to bt_CoReg_MR_CT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a=1;
set(handles.Skull.Hp,'parent',handles.ax1)
set(handles.FidLocalizer.SEH,'parent',handles.ax1);

if ~isempty(handles.FidLocalizer.SEH) && ~isempty(handles.FidLocalizer.MEH)
    SF=handles.FidLocalizer.SkullFid;
    MF=handles.FidLocalizer.MRIFid;
    SF0=[SF{1},SF{2},SF{3}];
    MF0=[MF{1},MF{2},MF{3}];
    
    
    % %Step 1 - Scale using averaged distances
    % tmp=[1,2;2,3;1,3];
    % for i=1:3
    %     a=tmp(i,1); b=tmp(i,2);
    %     dS(i)= sqrt( (SF0(1,a)-SF0(1,b))^2 + (SF0(2,a)-SF0(2,b))^2 + (SF0(3,a)-SF0(3,b))^2 );
    %     dM(i)= sqrt( (MF0(1,a)-MF0(1,b))^2 + (MF0(2,a)-MF0(2,b))^2 + (MF0(3,a)-MF0(3,b))^2 );
    % end
    % t1=mean(dM./dS);
    % M=makehgtform('scale',t1);
    % MM{1}=M;
    % MF0=bsxfun(@plus,M(1:3,1:3)*MF0,M(1:3,4));
    %Short Circuit
    MM{1}=[1,0,0,0;0,1,0,0;0,0,1,0;0,0,0,1];
    
    
    %Step 2 - Translate Both to Move Point 1 to Origin
    %First MF
    t2=-MF0(:,1);
    M = makehgtform('translate',t2);
    MM{2}=M;
    MF1=bsxfun(@plus,M(1:3,1:3)*MF0,M(1:3,4));
    
    %Next SF
    t2=-SF0(:,1);
    M = makehgtform('translate',t2);
    SM{1}=M;
    SF1=bsxfun(@plus,M(1:3,1:3)*SF0,M(1:3,4));
    
    
    %Step 3 - Rotate (2) onto z axis
    
    ssc = @(v) [0 -v(3) v(2); v(3) 0 -v(1); -v(2) v(1) 0];
    RU = @(A,B) eye(3) + ssc(cross(A,B)) + ssc(cross(A,B))^2*(1-dot(A,B))/(norm(cross(A,B))^2);
    %http://math.stackexchange.com/questions/180418/calculate-rotation-matrix-to-align-vector-a-to-vector-b-in-3d/476311#476311
    
    %First MF2
    u=MF1(:,2)/norm(MF1(:,2)); v=[0,0,1]';
    rm=RU(u,v);
    MF2=rm*MF1;
    MM{3}=[[rm;0,0,0],[0;0;0;1]];
    
    %Next SF2
    u=SF1(:,2)/norm(SF1(:,2)); v=[0,0,1]';
    rm=RU(u,v);
    SF2=rm*SF1;
    SM{2}=[[rm;0,0,0],[0;0;0;1]];
    
    
    %Step 4 - Project 3 onto (x,y) and rotate both 3's onto positive x axis
    %First MF3
    x=MF2(1,3); y=MF2(2,3);
    th=abs(atan(y/x));
    if x>0 && y>0, th=-th; end
    if x>0 && y<0, th=th; end
    if x<0 && y<0, th=-th-pi; end
    if x<0 && y>0, th=th+pi; end
    M = makehgtform('zrotate',th);
    MM{4}=M;
    MF3=M(1:3,1:3)*MF2;
    
    %Next SF3
    x=SF2(1,3); y=SF2(2,3);
    th=abs(atan(y/x));
    if x>0 && y>0, th=-th; end
    if x>0 && y<0, th=th; end
    if x<0 && y<0, th=-th-pi; end
    if x<0 && y>0, th=th+pi; end
    M = makehgtform('zrotate',th);
    SM{3}=M;
    SF3=M(1:3,1:3)*SF2;
    
    
    % Step 5 - Apply transforms to surfaces
    Svert=[get(handles.Skull.Hp,'vertices')]';
    Cvert=[get(handles.Cortex.Hp,'vertices')]';
    
    % for i=1:length(MM)
    %     Cvert=bsxfun(@plus,MM{i}(1:3,1:3)*Cvert,MM{i}(1:3,4));
    % end
    for i=1:length(SM)
        Svert=bsxfun(@plus,SM{i}(1:3,1:3)*Svert,SM{i}(1:3,4));
    end
    SFend=SF3;
    for i=length(MM):-1:1
        M=MM{i}^-1;
        Svert=bsxfun(@plus,M(1:3,1:3)*Svert,M(1:3,4));
        SFend=bsxfun(@plus,M(1:3,1:3)*SFend,M(1:3,4));
    end
    set(handles.Skull.Hp,'Vertices',Svert')
    set(handles.Cortex.Hp,'Vertices',Cvert')
    handles.Cortex.cortex.vert=Cvert';
    
    axes(handles.ax1)
    
    % elecmatrix=handles.El.elecmatrix';
    % for i=1:length(MM)
    %     elecmatrix=bsxfun(@plus,MM{i}(1:3,1:3)*elecmatrix,MM{i}(1:3,4));
    % end
    % elecmatrix=elecmatrix';
    % handles.El.elecmatrix=elecmatrix;
    % delete(handles.El.He);
    % handles.El.He=plot3(elecmatrix(:,1),elecmatrix(:,2),elecmatrix(:,3),'.','MarkerSize',35);
    
    delete(handles.FidLocalizer.SEH)
    delete(handles.FidLocalizer.MEH)
    
    MFend=MF0; SFend=SFend;
    
    handles.FidLocalizer.MEH=plot3(MFend(1,:),MFend(2,:),MFend(3,:),'.','MarkerSize',30,'color','b');
    handles.FidLocalizer.SEH=plot3(SFend(1,:),SFend(2,:),SFend(3,:),'.','MarkerSize',30,'color','r');
    
end
set(handles.Cortex.Hp,'FaceColor',[1,.65,.65]);

guidata(hObject,handles)


% --- Executes on button press in bt_TrimCT.
function bt_TrimCT_Callback(hObject, eventdata, handles)
% hObject    handle to bt_TrimCT (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ax1)
axes(handles.ax3), axis on, set(gca,'color','none')
xlabel('x'), zlabel('z'), ylabel('y')
set(handles.tx_TrimMode,'ForegroundColor','r','string','Click Left Ear')
handles.Trim.step=1;
handles.Trim.h = [handles.ax3; get(handles.ax3,'children')];
set(handles.Trim.h,'ButtonDownFcn',@TrimCT_ButtonPress)
view([0,90])
set(handles.Camera.cth,'visible','off')
handles.Skull.vert=get(handles.Skull.Hp,'Vertices');
handles.Skull.tri=get(handles.Skull.Hp,'Faces');
guidata(hObject,handles)
a=1;

function TrimCT_ButtonPress (hObject,eventdata,handles)
handles=guidata(hObject);
cp=get(handles.ax3,'CurrentPoint');

if handles.Trim.step==1
    voi=cp(1,1);
    handles.Trim.ind{handles.Trim.step}=(handles.Skull.vert(:,1)<voi);
    set(handles.tx_TrimMode,'string','Click Right Ear','ForegroundColor','r')
elseif handles.Trim.step==2
    voi=cp(1,1);
    handles.Trim.ind{handles.Trim.step}=(handles.Skull.vert(:,1)>voi);
    set(handles.tx_TrimMode,'string','Click Tip of Nose','ForegroundColor','r')

elseif handles.Trim.step==3
    voi=cp(1,2);%y for nasion
    handles.Trim.ind{handles.Trim.step}=(handles.Skull.vert(:,2)>voi);
    set(handles.tx_TrimMode,'string','CT Trim Finished','ForegroundColor','k')
    
    ind2=cell(3,1);
    for i=1:length(handles.Trim.ind)
        ind1=handles.Trim.ind{i};
        x=1:length(ind1);
        face=x(ind1);
        dum=reshape(handles.Skull.tri,[],1);
        ind1=ismember(dum,face);
        ind1=reshape(ind1,[],3);
        ind1=sum(ind1,2);
        
        ind2{i}=ind1>2;
    end
    
    ind=logical(ind2{1}.*ind2{2}.*ind2{3});
    delete(handles.Skull.Hp)
    handles.Skull.Hp = patch('Vertices',handles.Skull.vert,'Faces',handles.Skull.tri(ind,[1,3,2]),'facecolor',[.65 .65 .65],'edgecolor','none',...
        'facelighting', 'gouraud', 'specularstrength', .25);
    set(handles.Skull.Hp,'parent',handles.ax3)
    
    axis equal
    handles.Camera.cth=cameratoolbar;

    camlight('headlight','infinite')
    set(gca,'DataAspectRatioMode','manual','PlotBoxAspectRatioMode','manual');
    set(gca,'camerapositionmode','manual','cameratargetmode','manual','cameraupvectormode','manual','cameraviewanglemode','manual')

    
    a=1;
end

handles.Trim.step=handles.Trim.step+1;

guidata(handles.figure1,handles)


% --- Executes on button press in pb_NudgeCortex.
function pb_NudgeCortex_Callback(hObject, eventdata, handles)
% hObject    handle to pb_NudgeCortex (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.Nudge.mode==1
    set(handles.tx_NudgeCortex,'string','Nudge Mode Off','ForegroundColor',[0,0,0])
    handles.Nudge.mode=0;
    ch=get(handles.ax1,'children');
    set([handles.ax1; ch],'ButtonDownFcn',[] )
    
    
else
    
    handles.Nudge.mode=1;
    set(handles.tx_NudgeCortex,'string','Nudge Mode On','ForegroundColor',[1,0,0]);
    ch=get(handles.ax1,'children');
    set([handles.ax1; ch],'ButtonDownFcn',@NudgeCortex )
    axes(handles.ax1);
end

guidata(hObject,handles)




function NudgeCortex(hObject, eventdata, handles)
handles=guidata(hObject);
set(handles.tx_NudgeCortex,'string','Now Click Destination','ForegroundColor',[1,0,0]);
cp1=get(handles.ax1,'CurrentPoint');
handles.Nudge.cp1=cp1(1,:);
handles.Nudge.p1h=plot3(cp1(1,1),cp1(1,2),cp1(1,3),'.','MarkerSize',20,'Color',[0,1,1]);
ch=get(handles.ax1,'children');
set([handles.ax1; ch],'ButtonDownFcn',@NudgeCortex2 )


a=1;
guidata(hObject,handles);


function NudgeCortex2(hObject, eventdata, handles)
handles=guidata(hObject);
cp2=get(handles.ax1,'CurrentPoint');
cp2=cp2(1,:);
cp1=handles.Nudge.cp1;
cpd=cp2-cp1;

if ~isfield(handles.Nudge,'vert')
    handles.Nudge.vert=get(handles.Cortex.Hp,'vertices');
end

vert=get(handles.Skull.Hp,'vertices');
vert=bsxfun(@minus,vert,cpd);
set(handles.Skull.Hp,'vertices',vert)

delete(handles.Nudge.p1h);
set(handles.tx_NudgeCortex,'string','Nudge Mode Off','ForegroundColor',[0,0,0])
handles.Nudge.mode=0;
ch=get(handles.ax1,'children');
set([handles.ax1; ch],'ButtonDownFcn',[] )

guidata(hObject,handles)


% --- Executes on selection change in lb_NudgeDir.
function lb_NudgeDir_Callback(hObject, eventdata, handles)
% hObject    handle to lb_NudgeDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns lb_NudgeDir contents as cell array
%        contents{get(hObject,'Value')} returns selected item from lb_NudgeDir


% --- Executes during object creation, after setting all properties.
function lb_NudgeDir_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lb_NudgeDir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function sl_NudgeScale_Callback(hObject, eventdata, handles)
% hObject    handle to sl_NudgeScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.tx_NudgeScale,'String',num2str(get(hObject,'value')))
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function sl_NudgeScale_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_NudgeScale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pb_NudgeL.
function pb_NudgeL_Callback(hObject, eventdata, handles)
% hObject    handle to pb_NudgeL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dir=get(handles.lb_NudgeDir,'value');
sc=get(handles.sl_NudgeScale,'value')*handles.Nudge.ScaleMax;
vert=get(handles.Skull.Hp,'vertices');
vert(:,dir)=vert(:,dir)-sc;
set(handles.Skull.Hp,'vertices',vert);
if dir==1, a='xdata',
elseif dir==2, a='ydata',
elseif dir==3, a='zdata', end
set(handles.FidLocalizer.SEH2,a, get(handles.FidLocalizer.SEH2,a)-sc  )

guidata(hObject,handles);


% --- Executes on button press in pb_NudgeR.
function pb_NudgeR_Callback(hObject, eventdata, handles)
% hObject    handle to pb_NudgeR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dir=get(handles.lb_NudgeDir,'value');
sc=get(handles.sl_NudgeScale,'value')*handles.Nudge.ScaleMax;
vert=get(handles.Skull.Hp,'vertices');
vert(:,dir)=vert(:,dir)+sc;
set(handles.Skull.Hp,'vertices',vert);
if dir==1, a='xdata',
elseif dir==2, a='ydata',
elseif dir==3, a='zdata', end
set(handles.FidLocalizer.SEH2,a, get(handles.FidLocalizer.SEH2,a)+sc  )

guidata(hObject,handles);


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
val = get(hObject,'value');
if val ==1
    set(handles.Cortex.Hh,'visible','on','color',[.65,.65,.65])
else
    set(handles.Cortex.Hh,'visible','off')
end


% --------------------------------------------------------------------
function mn_Image_Callback(hObject, eventdata, handles)
% hObject    handle to mn_Image (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mn_FlipFlH_Callback(hObject, eventdata, handles)
% hObject    handle to mn_FlipFlH (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles.Fl,'HI')
   handles.Fl.I = flip(handles.Fl.I,2) ;
   set(handles.Fl.HI, 'cdata', handles.Fl.I);
end


% --------------------------------------------------------------------
function mn_FlipSurfX_Callback(hObject, eventdata, handles)
% hObject    handle to mn_FlipSurfX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
I = eye(3);
I(1)=-1;
if handles.SurfFlip == 0, handles.SurfFlip=1; else handles.SurfFlip = 0; end

set(handles.Cortex.Hp,'Vertices', get(handles.Cortex.Hp,'vertices') * I)
set(handles.Skull.Hp,'Vertices', get(handles.Skull.Hp,'Vertices') * I)
set(handles.Cortex.Hh,'xdata', get(handles.Cortex.Hh,'xdata') * -1)
for i = 1:length(handles.El.He)
   set(handles.El.He(i),'XData', get(handles.El.He(i),'Xdata') * -1)
end
set(handles.FidLocalizer.SEH2,'XData', get(handles.FidLocalizer.SEH2,'XData') * -1)

guidata(hObject,handles)


% --------------------------------------------------------------------
function mn_SetXCamPos_Callback(hObject, eventdata, handles)
% hObject    handle to mn_SetXCamPos (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = get(handles.ax1,'cameraposition');
set(handles.ax1,'cameraposition',[750,a(2),a(3)]);


% --- Executes on slider movement.
function sl_XCamLoc_Callback(hObject, eventdata, handles)
% hObject    handle to sl_XCamLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% a = get(handles.ax1,'cameraposition');
% if get(handles.rb_elecSideL,'value'), a(1) = 1*a(1); end
% set(handles.ax1,'cameraposition', [get(hObject,'value'),a(2),a(3)]);
% set(handles.tx_XLoc,'string',num2str(get(hObject,'value')));

% Latane Bullock 2024 01 04
% this function was previosly used to change the depth of the camera--we'll
% now use it to change the view angle
set(handles.ax1,'CameraViewAngle', get(hObject,'value'));

guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function sl_XCamLoc_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sl_XCamLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function mn_LoadAll_Callback(hObject, eventdata, handles)
% hObject    handle to mn_LoadAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.ax1), cla, hold on

% load other images
fsroot = dbs_getrecentfsfolder;
ptroot = dbs_getrecentptfolder;
[folders, files] = dbs_subdir(fsroot);
els_dir = folders(~cellfun(@isempty,strfind(folders,'lectrode')));
[~,els_files] = dbs_subdir(char(els_dir));
%load([fsroot, filesep, cell2mat(files(~cellfun(@isempty,strfind(files,'cortex'))))]) % cortex_indiv.mat
load([fsroot, filesep, 'cortex_indiv.mat']) % cortex_indiv.mat
load([fsroot, filesep, cell2mat(files(~cellfun(@isempty,strfind(files,'hull'))))])
try
load(fullfile(fsroot,'CT_reg','skull.mat'))
catch
load(fullfile(fsroot,'skull.mat'))
end

try
% load(char(fullfile(els_dir,els_files(~cellfun(@isempty,strfind(els_files,'electrodes'))))))
load(fullfile(fsroot, 'Electrode_Locations', 'depthelectrodes.mat'))
catch
disp('Choose Depth Electrodes')
[depth_electrodes,elsPath,ext] = dbs_uigetfile(fsroot,'Choose Depth Electrodes');
load(fullfile(elsPath,[depth_electrodes,ext]))
end

try
% load(char(fullfile(els_dir,els_files(~cellfun(@isempty,strfind(els_files,'Pin'))))))
load(fullfile(fsroot, 'Electrode_Locations', 'MRI_Fiducials.mat'))
catch
disp('Choose PinTips');
[PinTips,pinPath,ext] = dbs_uigetfile(fsroot,'Choose PinTips');
load(fullfile(pinPath,[PinTips,ext]))
end

handles.Cortex.hull = mask_indices;
% a=[-1,0,0;0,-1,0;0,0,1];
% skull.vert=skull.vert'; skull.tri=skull.tri';
% skull.vert = [skull.vert*a];

handles.Cortex.Hp = patch('vertices',cortex.vert,'faces',cortex.tri(:,[1 3 2]),'facecolor',[1 .65 .65],'edgecolor','none');
handles.Skull.Hp = patch('vertices',skull.vert,'faces',skull.tri(:,[1 3 2]),'facecolor',[.65 .65 .65],'edgecolor','none');
% handles.Cortex.Hp = patch('vertices',cortex.vert,'faces',cortex.tri(:,[1 3 2]),'facecolor',[1 .65 .65],'edgecolor','none',...
%     'facelighting', 'gouraud', 'specularstrength', .25);
% handles.Skull.Hp = patch('vertices',skull.vert,'faces',skull.tri(:,[1 3 2]),'facecolor',[.65 .65 .65],'edgecolor','none',...
%     'facelighting', 'gouraud', 'specularstrength', .25);
handles.Cortex.Hh = plot3(handles.Cortex.hull(:,1),handles.Cortex.hull(:,2),handles.Cortex.hull(:,3),'.','color','m');
set(handles.Cortex.Hh,'visible','off')
camlight('headlight'); lighting gouraud;

elecmatrix1 = elecmatrix(elecmatrix(:,1)<0,:);
elecmatrix2 = elecmatrix(elecmatrix(:,1)>0,:);
if ~isempty(elecmatrix1)
    [~,si] = sort(elecmatrix1(:,3));
    elecmatrix1 = elecmatrix1(si,:);
%     for i=1:3
%         elecmatrix1(:,i) = smooth(elecmatrix1(:,i),3) ;
%     end
    handles.El.He(1) = plot3(elecmatrix1(:,1),elecmatrix1(:,2),elecmatrix1(:,3),'linewidth',3,'color','r');
end
if ~isempty(elecmatrix2)
     [~,si] = sort(elecmatrix2(:,3));
    elecmatrix2 = elecmatrix2(si,:);
%     for i=1:3
%         elecmatrix2(:,i) = smooth(elecmatrix2(:,i),3) ;
%     end
    handles.El.He(2) = plot3(elecmatrix2(:,1),elecmatrix2(:,2),elecmatrix2(:,3),'linewidth',3,'color','b');
end

if isfield(handles.FidLocalizer,'SEH2');
    delete(handles.FidLocalizer.SEH2)
end
axes(handles.ax1);
nf = size(fiducial_locations,2);

colormap jet
cmap = colormap;
c=cmap(1:floor(length(cmap)/nf):end,:);
[~,si] = sortrows(fiducial_locations, [1 2]);
fiducial_locations=fiducial_locations(si,:);
handles.FidLocalizer.SEH2 = scatter3(fiducial_locations(:,1),fiducial_locations(:,2),fiducial_locations(:,3),50,c,'+','linewidth',3);
set(handles.FidLocalizer.SEH2,'cdata',[1,0,0;1,0,0;0,0,1;0,0,1]);


% populate the landmarks table with fiducial locations
% handles.Landmarks = table(); 
lm = handles.FluoroLocalizer.Landmarks;
nlm = []; % new landmark 
nlm.coordframe = "recon"; % 'recon' or 'fluoro'
for i = 1:height(fiducial_locations)
    nlm.pos = fiducial_locations(i, :);
    if nlm.pos(1)>0; nlm.hemi = "right";
    else nlm.hemi = "left"; end
    if nlm.pos(2)>0; nlm.name = "pintip_front";
    else nlm.name = "pintip_occ"; end
    if isempty(lm)
        lm = struct2table(nlm);
    else
        lm = [lm; struct2table(nlm)];
    end
end
% add dbslead_end and dbslead_mid
nlm.hemi = "left"; 
nlm.name = "dbslead_bottom"; 
elecpos = elecmatrix1;
[~, is] = min(elecpos(:,3)); 
nlm.pos = elecpos(is, :); 
lm = [lm; struct2table(nlm)];
nlm.name = "dbslead_top"; 
[~, is] = max(elecpos(:,3)); 
nlm.pos = elecpos(is, :); 
lm = [lm; struct2table(nlm)];
nlm.name = "dbslead_mid"; 
[~, is] = min(abs(elecpos(:,3) - mean(elecpos(:,3)))); 
nlm.pos = elecpos(is, :); 
lm = [lm; struct2table(nlm)];

nlm.hemi = "right"; 
nlm.name = "dbslead_bottom"; 
elecpos = elecmatrix2;
[~, is] = min(elecpos(:,3)); 
nlm.pos = elecpos(is, :); 
lm = [lm; struct2table(nlm)];
nlm.name = "dbslead_top"; 
[~, is] = max(elecpos(:,3)); 
nlm.pos = elecpos(is, :); 
lm = [lm; struct2table(nlm)];
nlm.name = "dbslead_mid"; 
[~, is] = min(abs(elecpos(:,3) - mean(elecpos(:,3)))); 
nlm.pos = elecpos(is, :); 
lm = [lm; struct2table(nlm)];
lm = unique(lm); 
handles.FluoroLocalizer.Landmarks = lm;
guidata(hObject,handles);




axis equal
set(gca,'DataAspectRatioMode','manual','PlotBoxAspectRatioMode','manual');
set(gca,'camerapositionmode','manual','cameratargetmode','manual','cameraupvectormode','manual','cameraviewanglemode','manual')
set(gca,'projection','perspective')

xloc = get(handles.sl_XCamLoc,'value');

% set(handles.ax1, 'ylim',[-250,250],'xlim',[-100,100],'zlim',[-100,100])

camlight('headlight','infinite');

%% FLUORO IMAGE
handles = mn_LoadFluoro_Callback(hObject, eventdata, handles); 
% handles = plot_image_3d_space(hObject, eventdata, handles, []); 


%% try to load Camera position
% handles = mn_LoadCamera_Callback(hObject, eventdata, handles); 

% datacursormode on
axes(handles.ax1)
scatter3(0,0,0,100, 'square', 'filled'); % plot 0,0,0 just as a reference point
% plot coordinate frame
ijk = [1 0 0; 0 1 0; 0 0 1]*25; 
orig = zeros(3, 3); 
quiver3(orig(:, 1), orig(:, 2), orig(:, 3), ...
        ijk(:, 1),  ijk(:, 2),  ijk(:, 3)); 
pts = orig + ijk; 
scatter3(pts(:, 1), pts(:, 2), pts(:, 3), 100, ijk, 'square', 'filled'); % plot 0,0,0 just as a reference point



% Set stack
set(handles.bt_FluoroStack,'value',0)
set(handles.bt_SkullStack,'value',0)
axes(handles.ax1);


% sl_Cort_Callback(hObject, eventdata, handles); 
% sl_Fluoro_Callback(hObject, eventdata, handles); 
% sl_Skull_Callback(hObject, eventdata, handles); 

set(handles.Cortex.Hp,'FaceAlpha',0)
set(handles.Skull.Hp,'FaceAlpha',0)
set(handles.Fl.HI,'AlphaData',0)


guidata(hObject,handles)


function handles = plot_image_3d_space(hObject, eventdata, handles, transmat) 
% optionally, pass transmat that specifies a transformation to a new
% coordinate frame (4x4)
% plot fluoro in 3d along with skull

axes(handles.ax1)
hold on; 

if isfield(handles.FluoroLocalizer, 'hS_fluoro'); delete(handles.FluoroLocalizer.hS_fluoro); end
if isfield(handles.FluoroLocalizer, 'hS_transparent'); delete(handles.FluoroLocalizer.hS_transparent); end

% flip = double(handles.FluoroLocalizer.emitterSide=="left")*2-1; 
xoff = (handles.FluoroLocalizer.emitter2screen - handles.FluoroLocalizer.emitter2nose); 
width = handles.FluoroLocalizer.width; 

corners = [
           xoff  width/2   width/2; % T L
           xoff  width/2  -width/2; % B L
           xoff -width/2   width/2; % T R
           xoff -width/2  -width/2; % B R
          ]; 


if ~isempty(transmat)
    corners = (transmat*padarray(corners', [1 0], 1, 'post')); 
    corners = corners(1:3, :)';
end
 

seg1 = corners(1:2, :); 
seg2 = corners(3:4, :); 
xImage = [seg1(:, 1) seg2(:, 1)]; 
yImage = [seg1(:, 2) seg2(:, 2)]; 
zImage = [seg1(:, 3) seg2(:, 3)]; 
handles.FluoroLocalizer.hS_fluoro = surf(xImage,yImage,zImage,...    % hS: handle to the surface
     'CData',handles.Fl.I,...
     'FaceColor','texturemap', 'FaceAlpha', 0.9);

% plot a transparent mesh--this is required to be able to clck on the locations in the fluoro
[z, y] = meshgrid(-(width/2):1:(width/2)); % Generate x and y data
x = zeros(size(y, 1)) + xoff-1; % Generate z data
% handles.FluoroLocalizer.hS_stransparent = surf(x, y, z, 'FaceAlpha', 0, 'EdgeColor', 'none'); % hS: handle to the surface

delete(findall(handles.ax1, 'type', 'light'))
camlight(handles.ax1, 'headlight', 'infinite')

guidata(hObject, handles)


% --------------------------------------------------------------------
function mn_CD_Callback(hObject, eventdata, handles)
% hObject    handle to mn_CD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
dir = uigetdir;
cd(dir);


% --- Executes on button press in LoadElectrodeLocs.
function LoadElectrodeLocs_Callback(hObject, eventdata, handles)
% hObject    handle to LoadElectrodeLocs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[CR_fname,CR_pname,~] = uigetfile('*.mat');
load(fullfile(CR_pname,CR_fname))
%handles.Cortex.cortex=cortex;  NEEDED FOR MANUAL COREGISTRATION, REMOVED
%FOR NOW
axes(handles.ax1)

MarkerSize = 50;
if length(CortElecLoc)>8
   MarkerSize=25;
end

els = cell2mat(CortElecLoc');
hold on; plot3(els(:,1),els(:,2),els(:,3),'.','color',[0 175/255 0],'markersize',MarkerSize)
hold on; plot3(els(:,1),els(:,2),els(:,3),'.','color',[0 175/255 0],'markersize',MarkerSize)
%handles.Cortex.Hp.Facealpha = 1;

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function mn_IndividualItems_Callback(hObject, eventdata, handles)
% hObject    handle to mn_IndividualItems (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function mn_CortElecLocs_Callback(hObject, eventdata, handles)
% hObject    handle to mn_CortElecLocs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[fname,pname,~] = uigetfile('*.mat');
load(fullfile(pname,fname),'CortElecLoc');
handles.CortLocalizer.ElecLocCort = CortElecLoc;
if ~isempty (handles.CortLocalizer.CEH), delete(findobj('type','line')), handles.CortLocalizer.CEH = [];  end
for i = 1:length(handles.CortLocalizer.ElecLocCort)
    A0 = handles.CortLocalizer.ElecLocCort{i};
    handles.CortLocalizer.CEH(i)=plot3(A0(1),A0(2),A0(3),'.','MarkerSize',12,'Color','b');
end

guidata(hObject,handles)


% --- Executes on button press in InterpolateElectrodes.
function InterpolateElectrodes_Callback(hObject, eventdata, handles)
% hObject    handle to InterpolateElectrodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fsroot = dbs_getrecentfsfolder;
load([fsroot filesep 'hull.mat'], 'mask_indices')

N_recent_contacts = str2num(handles.N_recent_contacts.String);

cfg = [];
cfg.CortElecLoc = handles.CortLocalizer.ElecLocCort((end-N_recent_contacts+1):end);
cfg.mask_indices = mask_indices;
cfg.Nrow = str2num(handles.Nrow.String);
cfg.NperR = str2num(handles.NperR.String);

[Elec] = Elec_equalization_Fcn(cfg);

delete(handles.CortLocalizer.CEH((end-N_recent_contacts+1):end));
handles.CortLocalizer.CEH = ...
    handles.CortLocalizer.CEH(1:(end-N_recent_contacts));
handles.CortLocalizer.numelec = ...
    handles.CortLocalizer.numelec - N_recent_contacts + length(Elec);
handles.CortLocalizer.ElecLocCort((end-N_recent_contacts+1):end) = [];
handles.CortLocalizer.ElecLocCort = cat( 2, handles.CortLocalizer.ElecLocCort, Elec );
handles.CortLocalizer.ElecLoc((end-N_recent_contacts+1):end) = [];
handles.CortLocalizer.ElecLoc = cat( 2, handles.CortLocalizer.ElecLoc, Elec );

axes(handles.ax1)
MarkerSize = 50;
if length(Elec)>8
   MarkerSize=25;
end

for i = 1:length(Elec)
    els = Elec{i};
    hold on;
    handles.CortLocalizer.CEH(end+1) = ...
        plot3(els(1),els(2),els(3),'.','color',[0 175/255 0],'markersize',MarkerSize);
end

guidata(hObject,handles)








function Nrow_Callback(hObject, eventdata, handles)
% hObject    handle to Nrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nrow as text
%        str2double(get(hObject,'String')) returns contents of Nrow as a double


% --- Executes during object creation, after setting all properties.
function Nrow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nrow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NperR_Callback(hObject, eventdata, handles)
% hObject    handle to NperR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NperR as text
%        str2double(get(hObject,'String')) returns contents of NperR as a double


% --- Executes during object creation, after setting all properties.
function NperR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NperR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on key press with focus on bt_CoRegCort1 and none of its controls.
function bt_CoRegCort1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to bt_CoRegCort1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)


% --- If Enable == 'on', executes on mouse press in 5 pixel border.
% --- Otherwise, executes on mouse press in 5 pixel border or over bt_CoRegCort1.
function bt_CoRegCort1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to bt_CoRegCort1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function N_recent_contacts_Callback(hObject, eventdata, handles)
% hObject    handle to N_recent_contacts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N_recent_contacts as text
%        str2double(get(hObject,'String')) returns contents of N_recent_contacts as a double


% --- Executes during object creation, after setting all properties.
function N_recent_contacts_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N_recent_contacts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MarkFluoro.
function MarkFluoro_Callback(hObject, eventdata, handles)
% hObject    handle to MarkFluoro (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set(handles.figure1,'pointer','crosshair');

set(gcf,'CurrentCharacter',char(1));

h=datacursormode;
h.Enable = 'on';

set(h,'DisplayStyle','datatip');

waitfor(gcf,'CurrentCharacter',char(32));

pos = getCursorInfo(h);
pos = pos.Position; 

set(handles.figure1,'pointer','arrow');

% % handles.FluoroLocalizer.hDataCursor = datacursormode; 
% point = drawpoint(handles.ax1); 
% 
% % % dcs = obj.DataCursors;
% pos = point.Position;   %Position of 1st cursor


% output_txt{1} = ['X: ', num2str(pos(1))];
% output_txt{2} = ['Y: ', num2str(pos(2))]; %this is the text next to the cursor
% 
% 
% prompt = {'Which landmark is this? {Fluoro, Recon}_{F, O, DBS}{L,R}'};
% prompt = {'space={recon, fluoro}', ...
%           'hemi={left, right}', ...
%         'name={pintip_front, pintip_occ}'};
% dlgtitle = 'Landmark id';
% % fieldsize = [1 40];
% % definput = {'30'};
% % opts.Interpreter = 'tex';
% answer = listdlg(prompt,dlgtitle);


space = {"coordframe", "Choose space:", ["fluoro", "recon"]}; 
hemi = {"hemi", "Choose hemisphere", ["uncertain", "left", "right"]};
name = {"name", "Choose location ID: ", ["pintip_front", "pintip_occ", "dbslead_top", "dbslead_mid", "dbslead_bottom"]}; 
dlgs = {space, hemi, name}; 
answer = {}; 
for i=1:3
    [indx,tf] = listdlg('PromptString', dlgs{i}{2}, 'ListString', dlgs{i}{3}, 'SelectionMode','single');
    answer{i} = dlgs{i}{3}(indx);
end
nlm = []; 
nlm.coordframe = answer{1};
nlm.hemi = answer{2};
nlm.name = answer{3};
nlm.pos = pos;
handles.FluoroLocalizer.Landmarks = [handles.FluoroLocalizer.Landmarks; struct2table(nlm)]; 


% list = {'Fluoro_FR','Fluoro_FL','Fluoro_OR', 'Fluoro_OL',...                   
%         'Fluoro_DBS1L','Fluoro_DBS2L'};
% [indx,tf] = listdlg('ListString',list);
% landmark_id = list(find(indx, 1)); 
% 
% handles.Landmarks.(landmark_id{1}) = pos; 

scatter3(handles.ax1, pos(1), pos(2), pos(3),100, 'square', 'filled'); % plot 0,0,0 just as a reference point

guidata(hObject, handles)

h.Enable = 'off';

fprintf('Created new landmark "%s"\n', answer{1}); 
handles.FluoroLocalizer.Landmarks




% 
% % handles.CortLocalizer.ax1ch=[handles.ax1;get(handles.ax1,'children')];
% if handles.FluoroLocalizer.Status1==0
%     handles.FluoroLocalizer.Status1=1;
% 
%     axes(handles.ax1)
% %     set(handles.CortLocalizer.ax1ch,'ButtonDownFcn',@DrawElecLine)
%     set(handles.tx_LocMode,'string','Mark fluoro Mode ON','ForegroundColor',[1,0,0])
% %     set(handles.figure1,'pointer','crosshair')
% %     handles.FluoroLocalizer.hP = drawpoint; 
% 
% %     set(handles.FluoroLocalizer.hDataCursor, 'UpdateFcn', {@labeldtips, hObject, handles});
%     
% elseif handles.FluoroLocalizer.Status1==1
%     handles.FluoroLocalizer.Status1=0;
% 
% %     set(handles.CortLocalizer.ax1ch,'ButtonDownFcn',[]);
%     set(handles.tx_LocMode,'string','Mark fluoro Mode OFF','ForegroundColor',[0,0,0])
%     set(handles.figure1,'pointer','arrow')
% end

guidata(hObject, handles)


function output_txt  = labeldtips(src, eventdata, hObject, handles)
% Display an observation's Y-data and label for a data tip
% obj          Currently not used (empty)
% event_obj    Handle to event object

% handles = guidata(hObject); 

% % dcs = obj.DataCursors;
pos = src.Position;   %Position of 1st cursor


% output_txt{1} = ['X: ', num2str(pos(1))];
% output_txt{2} = ['Y: ', num2str(pos(2))]; %this is the text next to the cursor
% 
% 
% prompt = {'Which landmark is this? {Fluoro, Recon}_{F, O, DBS}{L,R}'};
prompt = {'space={recon, fluoro}', ...
          'hemi={left, right}', ...
        'name={pintip_front, pintip_occ}'};
dlgtitle = 'Landmark id';
% fieldsize = [1 40];
% definput = {'30'};
% opts.Interpreter = 'tex';
answer = inputdlg(prompt,dlgtitle);

nlm = []; 
nlm.coordframe = answer{1};
nlm.hemi = answer{2};
nlm.name = answer{3};
nlm.pos = pos;
handles.FluoroLocalizer.Landmarks = [handles.FluoroLocalizer.Landmarks; struct2table(nlm)]; 
% 
% list = {'Fluoro_FR','Fluoro_FL','Fluoro_OR', 'Fluoro_OL',...                   
%         'Fluoro_DBS1L','Fluoro_DBS2L'};
% [indx,tf] = listdlg('ListString',list);
% landmark_id = list(find(indx, 1)); 
% 
% handles.Landmarks.(landmark_id{1}) = pos; 

scatter3(handles.ax1, pos(1), pos(2), pos(3),100, 'square', 'filled'); % plot 0,0,0 just as a reference point

guidata(hObject, handles)

fprintf('Created new landmark "%s"\n', answer{1}); 
handles.FluoroLocalizer.Landmarks


% function project_points(hObject, eventdata, handles)
% 
% cp = get(handles.ax1, 'CameraPosition'); 
% u = handles.Landmarks{2, 'pos'} - cp;
% N = cp; 
% n = [300, 0, 0];
% M = [300, 0, 0];
% pos = line_plane_intersection(u, N, n, M);
% scatter3(handles.ax1, pos(1), pos(2), pos(3),100, 'square', 'filled'); % plot 0,0,0 just as a reference point




% --- Executes on button press in pushbutton20.
function Button_Optimize_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% h for 'home'
% i for initial 
% p for 'prime' 


% handles = bt_SaveCamera_Callback(hObject, eventdata, handles);
bt_SaveCamera_Callback(hObject, eventdata, handles);

% % matched points if we know ID

lm = unique(handles.FluoroLocalizer.Landmarks);
lm = lm(ismember(lm.name, handles.FluoroLocalizer.LandmarksNames), :); % filter out non-anchor landmarks (ie, ecog strips)
nlm = []; nlm.coordframe = 'fluoro'; nlm.hemi = "uncertain"; nlm.pos = [nan nan nan]; 
if sum(lm.coordframe=="fluoro" & lm.name=="pintip_front")<2
    nlm.name = "pintip_front";
    lm = [lm; struct2table(nlm)]; 
end
if sum(lm.coordframe=="fluoro" & lm.name=="pintip_occ")<2
    nlm.name = "pintip_occ";
    lm = [lm; struct2table(nlm)]; 
end
assert(sum(lm.hemi=="uncertain")==4); 
assert(sum(lm.name=="pintip_front")==4); 
assert(sum(lm.name=="pintip_occ")==4); 
% this should ensure that dbslead_bottom is at the top of the list--this is KEY to an assumption in the projection_cost_fcn that the first row is the DBS lead
lm = sortrows(lm, 'name');
assert(lm.name(1)=="dbslead_bottom")


hemis = ["left", "right"]; 
hemis_perm = [hemis([1 2]) hemis([1 2]); 
              hemis([2 1]) hemis([1 2]); 
              hemis([1 2]) hemis([2 1]); 
              hemis([2 1]) hemis([2 1])]; 
% perms = []; 
for perm_idx = 4:-1:1

lm_recon = lm(lm.coordframe=="recon", :); 
lm_fluoro = lm(lm.coordframe=="fluoro", :); 
idxs = lm_fluoro.hemi=="uncertain" ;
lm_fluoro.hemi(idxs) = hemis_perm(perm_idx, :); 
lm_matched = join(lm_fluoro, lm_recon, 'Keys', {'hemi', 'name'}); 
h.landmarks_fluoro = lm_matched.pos_lm_fluoro'; % in 3d space
h.landmarks_recon = lm_matched.pos_lm_recon'; % in 3d space
% OR--unmatched points
% h.landmarks_fluoro = lm_fluoro.pos'; % in 3d space
% h.landmarks_recon = lm_recon.pos'; % in 3d space

 % pad for extra dimension of homogenous coordinate matrices
h.landmarks_fluoro = padarray(h.landmarks_fluoro, [1 0], 1, 'post'); 
h.landmarks_recon = padarray(h.landmarks_recon, [1 0], 1, 'post'); 

h.cam_pos = [-handles.FluoroLocalizer.emitter2nose; 0; 0; 1]; 
h.cam_target = [0; 0; 0; 1]; 
h.cam_upvec = [0; 0; 1; 1]; 

fluoro_depth = handles.FluoroLocalizer.emitter2screen - ...
               handles.FluoroLocalizer.emitter2nose; % distance from camtarget to fluoro
h.fluoro_norm = (h.cam_target - h.cam_pos); 
h.fluoro_norm = h.fluoro_norm/norm(h.fluoro_norm);
h.fluoro_norm(4) = 0; % don't translate normal vector to plane
% flip = double(handles.FluoroLocalizer.emitterSide=="left")*2-1; 
h.fluoro_origin = (h.fluoro_norm*fluoro_depth - h.cam_target); 
h.fluoro_origin(4) = 1; 

% used in calculating distance for landmarks NOT in fluoro field of view 
h.radius_of_view = handles.FluoroLocalizer.fov_radius;

bt_ResetCamera_Callback(hObject, eventdata, handles); 
i.cam_pos = [get(handles.ax1, 'CameraPosition')'; 1]; 
i.cam_target = [get(handles.ax1, 'CameraTarget')'; 1];
i.cam_upvec = [get(handles.ax1, 'CameraUpVector')'; 1];
% i.cam_pos = h.cam_pos; 
% i.cam_target = h.cam_target; 
% i.cam_upvec = h.cam_upvec; 

i.rotmat = cam2rigidtransmat(i.cam_pos(1:3), i.cam_target(1:3), i.cam_upvec(1:3));
i.pyr = rotmat2angles(i.rotmat(1:3, 1:3)); 

% Set parameters ranges for the optimization
range_rot = 2; 
range_trans = 100; 
pitch = optimvar("pitch","LowerBound",i.pyr(1)-range_rot,"UpperBound",i.pyr(1)+range_rot);
yaw = optimvar("yaw","LowerBound",i.pyr(2)-range_rot,"UpperBound", i.pyr(2)+range_rot);
roll = optimvar("roll","LowerBound",i.pyr(3)-range_rot,"UpperBound",i.pyr(3)+range_rot);
tx = optimvar("tx","LowerBound",-range_trans,"UpperBound",range_trans);
ty = optimvar("ty","LowerBound",-range_trans,"UpperBound",range_trans);
tz = optimvar("tz","LowerBound",-range_trans,"UpperBound",range_trans);
alpha = optimvar("alpha","LowerBound", 0.8,"UpperBound",1.2);
% alpha = 1; 


% % % -------- test initial value--------
% i.alpha = 1;
% [cost, p, rigidtransmat] = projection_cost_fcn(i.pyr(1), i.pyr(2), i.pyr(3), i.cam_target(1), i.cam_target(2), i.cam_target(3), i.alpha, h);
% handles.Camera.cp = p.cam_pos(1:3);
% handles.Camera.ct = p.cam_target(1:3);
% handles.Camera.uv = [0 0 1];
% bt_ResetCamera_Callback(hObject, eventdata, handles)
% 
% handles = plot_image_3d_space(hObject, eventdata, handles, rigidtransmat); 
% h_Sc = gobjects(); 
% for ip = 1:width(p.landmarks_fluoro)
% h_Sc(end+1) = scatter3(p.landmarks_fluoro(1,ip)-2, ...
%          p.landmarks_fluoro(2,ip), ...
%          p.landmarks_fluoro(3,ip), 100, 'green', 'filled'); 
% drawnow; pause(0.5)
% h_Sc(end+1) = scatter3(p.landmarks_recon_proj(1,ip)-2, ...
%          p.landmarks_recon_proj(2,ip), ...
%          p.landmarks_recon_proj(3,ip), 100, 'r', 'filled'); 
% drawnow; pause(1)
% end
% 
% cost
% fprintf('Initial point'); 
% delete(h_Sc); 
% delete([handles.FluoroLocalizer.hS_fluoro]); 
% % % -------- test initial value--------



% Set initial starting point for the solver
initialPoint.pitch = i.pyr(1);
initialPoint.yaw = i.pyr(2);
initialPoint.roll = i.pyr(3);
initialPoint.tx = i.cam_target(1);
initialPoint.ty = i.cam_target(2);
initialPoint.tz = i.cam_target(3);
initialPoint.alpha = 1;

% Create problem
problem = optimproblem;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp');

% Define problem objective
problem.Objective = fcn2optimexpr(@projection_cost_fcn,pitch,yaw,roll,tx,ty,...
    tz,alpha,h);

% Display problem information
show(problem);

[sol,fval,exitflag,output,lambda] = solve(problem, initialPoint, "Options", options); 

% Clear variables
% clearvars pitch yaw roll tx ty tz alpha

% plot solution: plot fluoro, then move camera
% T = angles2rotmat([sol.pitch sol.yaw sol.roll]); 
% T = [T        [sol.tx; sol.ty; sol.tz]; 
%      0 0 0     1];

[cost, p, affinetransmat, dist] = projection_cost_fcn(sol.pitch, sol.yaw, sol.roll, sol.tx, sol.ty, sol.tz, sol.alpha, h);
% [cost, p, rigidtransmat] = projection_cost_fcn(pitch, yaw, roll, tx, ty, tz, alpha, h);
set(handles.ax1, 'CameraPosition', p.cam_pos(1:3))
set(handles.ax1, 'CameraTarget', p.cam_target(1:3))
set(handles.ax1, 'CameraUpVector', affinetransmat(1:3,1:3)*[0; 0; 1]); 

[handles, answer] = validate_camera(hObject, eventdata, handles);

sol
cost
dist
perm_idx
drawnow; shg; 

switch answer
    case "Save"
        save(handles.FluoroLocalizer.OptimizationSolutionMATFileName, 'sol', 'cost', "dist", 'h', 'i', 'p', 'lm_matched');
%         lm_matched.sol_landmarks_recon_proj = p.landmarks_recon_proj';
%         lm_matched.sol_landmarks_fluoro = p.landmarks_fluoro';
%         lm_matched.sol_dist = dist;
%         writetable(lm_matched, ...
%             handles.FluoroLocalizer.OptimizationSolutionFileName, ...
%             'Delimiter', '\t', 'FileType', 'text');
        
        
        fprintf('Optimization finished\n');
        
        break
    case "No"
        % nothing, continue
    case "Cancel search"
        break
    otherwise
        error('Answer not recognized') 
end


end


function [handles, answer] = validate_camera(hObject, eventdata, handles)

i.cam_pos = [get(handles.ax1, 'CameraPosition')'; 1];
i.cam_target = [get(handles.ax1, 'CameraTarget')'; 1];
i.cam_upvec = [get(handles.ax1, 'CameraUpVector')'; 1];

i.rotmat = cam2rigidtransmat(i.cam_pos(1:3), i.cam_target(1:3), i.cam_upvec(1:3));
i.pyr = rotmat2angles(i.rotmat(1:3, 1:3));

% [cost, p, affinetransmat, dist] = projection_cost_fcn(i.pyr(1), i.pyr(2), i.pyr(3), i.rotmat(1, 4), i.rotmat(1, 4), i.rotmat(1, 4), alpha, h);
% [cost, p, rigidtransmat] = projection_cost_fcn(pitch, yaw, roll, tx, ty, tz, alpha, h);

lm = handles.FluoroLocalizer.Landmarks; 
lm_fluoro = lm(lm.coordframe=="fluoro", :);
pts = lm_fluoro.pos'; 
lm_fluoro.sol_pos = (i.rotmat*padarray(pts, [1 0], 1, 'post'))'; 
lm_fluoro.sol_pos = lm_fluoro.sol_pos(:, 1:3); 

% find intersection of current camera angle, projected through recon, onto
% fluoro
lm_recon = lm(lm.coordframe=="recon", :);
lm_recon.sol_pos = nan([height(lm_recon), 3]); 
for ip = 1:height(lm_recon)
lm_recon.sol_pos(ip, :) = line_plane_intersection(lm_recon.pos(ip, :)' - i.cam_pos(1:3), ...
                                            i.cam_pos(1:3), ...                          
                                            i.cam_target(1:3)-i.cam_pos(1:3), ... % norm
                                            lm_fluoro.sol_pos(1, 1:3)'); % any point on the fluoro
end

lm = [lm_recon; lm_fluoro]; 

% p.affinetransmat = affinetransmat;
% handles = plot_image_3d_space(hObject, eventdata, handles, affinetransmat);
% set(handles.ax1, 'CameraPosition', p.cam_pos(1:3))
% set(handles.ax1, 'CameraTarget', p.cam_target(1:3))
% set(handles.ax1, 'CameraUpVector', affinetransmat(1:3,1:3)*[0; 0; 1]);

handles = plot_image_3d_space(hObject, eventdata, handles, i.rotmat); 

h_validate = gobjects();
h_validate(1) = scatter3(lm.sol_pos(:,1), ...
        lm.sol_pos(:, 2), ...
        lm.sol_pos(:, 3), 100, 'green', 'filled');

% h_validate(end+1) = scatter3(p.landmarks_recon_proj(1,:), ...
%     p.landmarks_recon_proj(2,:), ...
%     p.landmarks_recon_proj(3,:), 100, 'filled', 'MarkerFaceColor', "#7E2F8E");
% 
% npts = size(p.landmarks_recon_proj, 2);
% cam = repmat(p.cam_pos, [1 npts]);
% h_validate(end+1:end+npts) = plot3([cam(1, :); p.landmarks_recon_proj(1, :)], ...
%     [cam(2, :); p.landmarks_recon_proj(2, :)], ...
%     [cam(3, :); p.landmarks_recon_proj(3, :)], 'color', "#7E2F8E", 'linewidth', 2);

handles.FluoroLocalizer.h_validate = h_validate; guidata(hObject, handles);
%     
drawnow; shg;

answer = questdlg('Save landmarks table with current camera view?', '', ...
    'Save', "No", "Cancel search", 'No');
switch answer
    case "Save"
%         save(handles.FluoroLocalizer.OptimizationSolutionMATFileName, 'sol', 'cost', "dist", 'h', 'i', 'p', 'lm_matched');
        handles.FluoroLocalizer.Landmarks = lm;
        guidata(hObject, handles); 

%         mn_SaveCamera_Callback(hObject, eventdata, handles);

    case "No"
        delete(h_validate);
    case "Cancel search"
        delete(h_validate);
end

guidata(hObject, handles)




% --- Executes on button press in pushbutton21.
function pushbutton21_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sfprintf('testing/debug here')

flip = double(handles.FluoroLocalizer.emitterSide=="left")*2-1; 
handles.Camera.cp = [-1*flip*handles.FluoroLocalizer.emitter2nose, 0, 0];
handles.Camera.ct = [0, 0, 0];
handles.Camera.cva = 20;
handles.Camera.uv = [0, 0, 1]; 
bt_ResetCamera_Callback(hObject, eventdata, handles)

bt_SaveCamera_Callback(hObject, eventdata, handles)



handles = validate_camera(hObject, eventdata, handles);



R = cam2rigidtransmat(get(handles.ax1, 'CameraPosition'), get(handles.ax1, 'CameraTarget'), get(handles.ax1, 'CameraUpVector')); 
plot_image_3d_space(hObject, eventdata, handles, R); 
% % plot circle to ensure estimate of fluoro field of view radius is correct
% xo = 0; yo = 0; 
% th = 0:pi/50:2*pi;
% r = handles.FluoroLocalizer.fov_radius;
% y = r * cos(th) + xo;
% z = r * sin(th) + yo;
% x = repmat(handles.FluoroLocalizer.emitter2screen - handles.FluoroLocalizer.emitter2nose, [1 numel(y)]);
% h = plot3(x, y, z);


if isfield(handles.FluoroLocalizer, 'hS_fluoro'); delete(handles.FluoroLocalizer.hS_fluoro); end
if isfield(handles.FluoroLocalizer, 'hS_transparent'); delete(handles.FluoroLocalizer.hS_transparent); end




bt_SaveCamera_Callback(hObject, eventdata, handles);

cammat = [handles.Camera.cp' handles.Camera.ct' handles.Camera.uv'];

rotmat = angles2rotmat([-0.26 0 0]); 
transmat = rotmat2rigidmat(rotmat, []); 
cammat_new = transmat*padarray(cammat, [1, 0], 1, 'post');
cammat_new = num2cell(cammat_new(1:3, :), 1); 
[handles.Camera.cp handles.Camera.ct handles.Camera.uv] = deal(cammat_new{:});
bt_ResetCamera_Callback(hObject, eventdata, handles);




% --- Executes on button press in rb_fluoro_visible.
function rb_fluoro_visible_Callback(hObject, eventdata, handles)
% hObject    handle to rb_fluoro_visible (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rb_fluoro_visible

% plot image on another axis
axes(handles.ax2)
isfluoro = get(handles.rb_fluoro_visible,'Value'); 
set(handles.Fl.HI, 'visible', isfluoro); 

set(handles.bt_FluoroStack,'value',0)
set(handles.bt_SkullStack,'value',0)
axes(handles.ax1);

guidata(hObject, handles)

% set(handles.bt_FluoroStack(2), 'value', ~isfluoro); 



% --------------------------------------------------------------------
function mn_ExportFluoroLandmarks_Callback(hObject, eventdata, handles)
% hObject    handle to mn_ExportFluoroLandmarks (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% 
% % ----- write landmarks
% lm_tbl_fname = handles.FluoroLocalizer.LandmarksFileName; 
% lm = handles.FluoroLocalizer.Landmarks; 
% writetable(lm, lm_tbl_fname, 'Delimiter', '\t', 'FileType', 'text'); 


