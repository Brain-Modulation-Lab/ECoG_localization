function dbs_autocoord(options)

% options=ea_resolve_elspec(options);
patientname=options.patientname;
options.prefs=dbs_prefs(patientname);

if ~strcmp(options.patientname,'No Patient Selected') || ~exist(options.patientname) || isempty(options.patientname) % only 3D-rendering viewer can be opened if no patient is selected.
    %     error('No Patient Selected')
end

if options.coregct.do
    [options] = eval(['dbs_coregct', options.coregct.method,'(options)']); % triggers the coregct function and passes the options struct to it.
end

if options.coregctcheck
    % export "control" niftis with wireframe of normal anatomy..
    dbs_show_ctcoregistration(options);
end

if options.preplocal
    dbs_preplocalization(options) %Create .mat Files for Skull and Cortex and Calculate Cortical Hull
end

% Run Fiducials Localizer
if options.fidlocalize.do
    if isempty(options.uipatdirs)
        error('Please Choose a Patient Directory')
    else
        cd(char(options.uipatdirs))
        save('localizer_options.mat','options')
        MRI_Fid_Localizer
        % dbs_fidlocalizer.m
    end
end

% Run Strip Localizer
if options.striplocalize.do
    if isempty(options.uipatdirs)
        error('Please Choose a Patient Directory')
    else
        cd(char(options.uipatdirs))
        save('localizer_options.mat','options')
%         cd(char(options.uifsdir)) % commented out Latane Bullock 2024 05 28
        hGUI_dbs_elec_localizer = DBS_Elec_Localizer;
        assignin('base', 'hGUI_dbs_elec_localizer', hGUI_dbs_elec_localizer);
    end
end







%     if options.normalize.do
%         try load([options.root,options.patientname,filesep,'ea_normmethod_applied']); end
%         if exist('norm_method_applied','var')
%             try
%                 norm_method_applied{end+1}=options.normalize.method;
%             catch
%                 clear norm_method_applied
%                 norm_method_applied{1}=options.normalize.method;
%             end
%         else
%             norm_method_applied{1}=options.normalize.method;
%         end
%         save([options.root,options.patientname,filesep,'ea_normmethod_applied'],'norm_method_applied');
%
%         eval([options.normalize.method,'(options)']); % triggers the normalization function and passes the options struct to it.
%     end


%     if options.dolc % perform lead connectome subroutine..
%         ea_perform_lc(options);
%     end



%     if options.atl.genpt % generate patient specific atlas set
%         ea_ptspecific_atl(options);
%
%     end
%
%     if options.atl.normalize % normalize patient-specific atlas-set.
%         ea_normalize_ptspecific_atl(options)
%     end



%     if options.normalize.check
%
%
%         % export "control" niftis with wireframe of normal anatomy..
%
%         ea_show_normalization(options);
%
%     end


%
%     if options.doreconstruction
%         ea_checkfiles(options);
%
%
%         for side=1:length(options.sides)
%             %try
%             % call main routine reconstructing trajectory for one side.
%             [coords,trajvector{side},trajectory{side},tramat]=ea_reconstruct(patientname,options,options.sides(side));
%
%
%
%             %% refit electrodes starting from first electrode (this is redundant at
%             %% this point).
%
%             coords_mm{side} = ea_map_coords(coords', [options.root,options.prefs.patientdir,filesep,options.prefs.tranii])';
%
%
%             [~,distmm]=ea_calc_distance(options.elspec.eldist,trajvector{side},tramat(1:3,1:3),[options.root,options.prefs.patientdir,filesep,options.prefs.tranii]);
%
%             comp = ea_map_coords([0,0,0;trajvector{side}]', [options.root,options.prefs.patientdir,filesep,options.prefs.tranii])'; % (XYZ_mm unaltered)
%
%             trajvector{side}=diff(comp);
%
%
%
%             normtrajvector{side}=trajvector{side}./norm(trajvector{side});
%
%             for electrode=2:4
%
%                 coords_mm{side}(electrode,:)=coords_mm{side}(1,:)-normtrajvector{side}.*((electrode-1)*distmm);
%
%             end
%             markers(side).head=coords_mm{side}(1,:);
%             markers(side).tail=coords_mm{side}(4,:);
%
%
%             orth=null(normtrajvector{side})*(options.elspec.lead_diameter/2);
%
%             markers(side).x=coords_mm{side}(1,:)+orth(:,1)';
%             markers(side).y=coords_mm{side}(1,:)+orth(:,2)'; % corresponding points in reality
%
%
%
%             coords_mm=ea_resolvecoords(markers,options);
%
%         end
%
%

%     if length(coords_mm)==4 % only one side was processed.
%         if options.sides==1
%             coords_mm=[nan(4,3);coords_mm];
%         elseif options.sides==2
%             coords_mm=[coords_mm;nan(4,3)];
%         end
%     end






% transform trajectory to mm space:
%         for side=1:length(options.sides)
%
%             try
%                 if ~isempty(trajectory{side})
%                     trajectory{side}=ea_map_coords(trajectory{side}', [options.root,patientname,filesep,options.prefs.tranii])';
%                 end
%
%             end
%         end
%
%
%         % save reconstruction results
%         elmodel=options.elmodel;
%         ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,0,options);
%
%     end


%     if options.manualheightcorrection
%         % load reconstruction results
%         %     try
%         %         [coords_mm,trajectory,markers,elmodel,manually_corrected]=ea_load_reconstruction(options);
%         %     catch
%         %         ea_error([patientname,': No reconstruction information found. Please run reconstruction first.']);
%         %     end
%         %     ea_save_reconstruction(coords_mm,trajectory,markers,elmodel,0,options);
%
%
%         mcfig=figure('name',[patientname,': Manual Height Correction'],'numbertitle','off');
%         warning('off');
%         try
%             ea_maximize(mcfig);
%         end
%
%         ea_manualreconstruction(mcfig,patientname,options);
%
%     else
%         ea_write(options)
%     end
%
%     dbs_write(options)
end






