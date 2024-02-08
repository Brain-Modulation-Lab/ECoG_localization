function dbs_preplocalization(options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare_for_Localization                                                %
% Michael Randazzo 06/29/2015                                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ptdirectory = char(options.uipatdirs);
fsdirectory = char(options.uifsdir);
cd(fsdirectory)
% Check for Electrode_Locations Folder
if ~exist([char(options.uifsdir) '/Electrode_Locations'])
    mkdir([char(options.uifsdir),'/Electrode_Locations']);
% if ~exist([directory '\Electrode_locations'])  % Windows
%     mkdir([directory,'\Electrode_locations']); % Windows
end

%%%%%%%%%%%% Cortex %%%%%%%%%%%%%%
[cortex.vert_lh,cortex.tri_lh]= read_surf(fullfile(fsdirectory,'surf/lh.pial.T1')); % Reading left side pial surface
[cortex.vert_rh,cortex.tri_rh]= read_surf(fullfile(fsdirectory,'surf/rh.pial.T1')); % Reading right side pial surface

% Generating entire cortex
cortex.vert = [cortex.vert_lh; cortex.vert_rh]; % Combining both hemispheres
cortex.tri = [cortex.tri_lh; (cortex.tri_rh + length(cortex.vert_lh))]; % Combining faces (Have to add to number of faces)

cortex.tri=cortex.tri+1; % freesurfer starts at 0 for indexing

% Reading in MRI parameters
% WJL edit 06/12/2017 
% f=MRIread(fullfile(fsdirectory,'mri/T1.mgz'));
f=MRIread(fullfile(fsdirectory,'mri/T1.nii'));

% Translating into the appropriate space
for k=1:size(cortex.vert,1)
    a=f.vox2ras/f.tkrvox2ras*[cortex.vert(k,:) 1]';
    cortex.vert(k,:)=a(1:3)';
end

save(fullfile(fsdirectory,'cortex_indiv.mat'),'cortex');

% Skull
% Get Mesh
% AB edit 2019/02/21
% PLB edit 2024 02 02 we are creating the skull reconstruction with
% isosurface() and saving to freesurfer/skull.mat in a previous step, 
% so this is unnecessary

%loading skull
[fname,pathname,extension] = dbs_uigetfile(fsdirectory,'Choose Skull Mesh .obj');
[skull.vert,skull.tri] = dbs_getskullobj([pathname filesep fname extension]);

%conforming to cortex structure
if exist('skull')
    if isfield(skull,'faces'), skull.tri=skull.faces; end
elseif exist('skin')
    skull=skin; clear skin
    if isfield(skull,'faces'), skull.tri=skull.faces; end
end
s=size(skull.vert);
if s(1)<s(2), skull.vert=skull.vert'; end
s=size(skull.tri);
if s(1)<s(2), skull.tri=skull.tri'; end

a=[-1,0,0;0,-1,0;0,0,1];
skull.vert = [skull.vert*a];

% figure;
% Hp = patch('Vertices',cortex.vert,'Faces',cortex.tri,...
%     'facecolor',[1 1 1],'edgecolor','none',...
%     'facelighting', 'gouraud', 'specularstrength', .50);
% camlight('headlight','infinite');
% axis off; axis equal
% alpha 1
% hold on
% Hp = patch('Vertices',skull.vert,'Faces',skull.tri,...
%     'facecolor',[1 1 1],'edgecolor','none',...
%     'facelighting', 'gouraud', 'specularstrength', .50);
% camlight('headlight','infinite');
% axis off; axis equal
% alpha 0.2
% hold on

% % Selecting co-registered preop ct from which obj was constructed to
% % extract vox2ras matrix
% [fname,pathname,extension] = dbs_uigetfile(ptdirectory,'co-registered preop ct');
% rpreopct = MRIread([pathname filesep fname extension]);

% Translating skull into the appropriate space
% skull.vert = skull.vert';
% skull.tri = skull.tri';
% for k=1:size(skull.vert,1)
%     a=f.vox2ras/f.tkrvox2ras*[skull.vert(k,:) 1]';
%     skull.vert(k,:)=a(1:3)';
% end
% skull.vert = skull.vert';
% skull.tri = skull.tri';

if ~exist(fullfile(fsdirectory,'ct_reg'))
    mkdir(fullfile(fsdirectory,'ct_reg'))
else
end
disp(['saving    ' fullfile(fsdirectory,'ct_reg','skull.mat') ' .......'])
save(fullfile(fsdirectory,'ct_reg','skull.mat'),'skull')

% Option to prevent overwriting:
% if ~exist(fullfile(directory,fname))
%     try
%         [skull.vert,skull.tri] = dbs_displayobj(fname);
%         save(fullfile(directory,'skull.mat'),'skull')
%     catch
%         cd(directory)
%         [skull.vert,skull.tri] = dbs_displayobj(fname);
%         save(fullfile(directory,'skull.mat'),'skull')
%     end
% else
% end





% Create cortical hull
grayfilename = [char(options.uifsdir) '/mri/t1_class.nii'];
outputdir='hull';
disp('running dbs_gethull.......')
[~,~] = dbs_gethull(grayfilename,outputdir,3,21,.3);
disp('** Process Done')
end

