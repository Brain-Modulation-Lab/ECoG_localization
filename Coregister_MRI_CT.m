Subject_Path = '\\136.142.16.9/Nexus/DBS/DBS4079/Anatomy/Strip_DBS4079';

% Coregister and reslice the CT image (source) to the MRI image (reference)
% SPM Commands
matlabbatch{1}.spm.spatial.coreg.estwrite.ref = {[Subject_Path,'/preopmri.nii,1']};
matlabbatch{1}.spm.spatial.coreg.estwrite.source = {[Subject_Path,'/postop_ct.nii,1']};

matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 1;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';

% Executes the SPM commands
spm('defaults', 'FMRI');
spm_jobman('serial',matlabbatch);
