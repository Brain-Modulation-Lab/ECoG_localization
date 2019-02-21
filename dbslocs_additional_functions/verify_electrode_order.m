% define path variables
SUBJECT='DBS3020';
PATH_DATA='\\136.142.16.9\Nexus\DBS';
PATH_FREESURFER = [PATH_DATA filesep SUBJECT '\Anatomy\FreeSurfer\preop'];
PATH_SYNC = [PATH_DATA filesep SUBJECT '\Preprocessed Data\Sync'];

% load cortical surface and electrode locations
load([PATH_FREESURFER filesep 'cortex_indiv'], 'cortex');
load([PATH_FREESURFER filesep 'Electrode_Locations\CortElecLocL_eq.mat'], 'CortElecLoc');

% load annot tables
electrode    = bml_annot_read([PATH_SYNC filesep 'annot' filesep SUBJECT '_electrode.txt']);

%Displaying cortex
%load('E:\MATLAB\lead_v2.1.5\templates\space\MNI_ICBM_2009b_NLIN_ASYM\cortex\CortexHiRes.mat', 'Vertices', 'Faces');
figure;
Hp = patch('Vertices',cortex.vert,'Faces',cortex.tri,...
    'facecolor',[1 1 1],'edgecolor','none',...
    'facelighting', 'gouraud', 'specularstrength', .50);
camlight('headlight','infinite');
axis off; axis equal
alpha 0.5

% plot electrodes
elec = reshape(cell2mat(CortElecLoc),3,length(CortElecLoc))';
for e=1:length(elec)
   hold on; plot3(elec(e,1), elec(e,2), elec(e,3), 'o', 'color', 'g', 'MarkerSize', 15)
   hold on; text(elec(e,1), elec(e,2), elec(e,3), num2str(e))
end

% flip electrodes within rows, if necessary
CortElecLoc = CortElecLoc([21:-1:1, 42:-1:22, 63:-1:43]);
save([PATH_FREESURFER filesep 'Electrode_Locations\CortElecLocL_eq.mat'], 'CortElecLoc', '-append');

elec = reshape(cell2mat(CortElecLoc),3,length(CortElecLoc))';

% % add columns to electrode table
% if ~any(ismember({'nat_x','nat_y','nat_z'},electrode.Properties.VariableNames))
%     electrode.nat_x(:) = nan;
%     electrode.nat_y(:) = nan;
%     electrode.nat_z(:) = nan;
% end


% add first 63 strip electrode coordinates to electrode table
strip_no = 1;
contact_range = 1:63;
ecog_s1 = contains(electrode.electrode, ['ecog_' num2str(strip_no)]);
electrode.nat_x(ecog_s1) = elec(contact_range,1);
electrode.nat_y(ecog_s1) = elec(contact_range,2);
electrode.nat_z(ecog_s1) = elec(contact_range,3);

% add second 63 strip electrode coordinates to electrode table
strip_no = 2;
contact_range = 64:126;
ecog_s2 = contains(electrode.electrode, ['ecog_' num2str(strip_no)]);
electrode.nat_x(ecog_s2) = elec(contact_range,1);
electrode.nat_y(ecog_s2) = elec(contact_range,2);
electrode.nat_z(ecog_s2) = elec(contact_range,3);

% write new electrode file
writetable(electrode, [PATH_SYNC filesep 'annot' filesep SUBJECT '_electrode.txt'], 'Delimiter', '\t');      
