%-TC: 04/22/2013: turn on the AR(1) option

function stats_fmri_fconnect_noscaling_spm12 (datadir, data_type, NUMTRUNC, imagefilter)


spm('Defaults', 'fmri');
statsdir = pwd;

% fMRI stats batchtemplate

load('/batch_stats_spm12.mat');
load '/contrasts.mat';

% Update TR
matlabbatch{1}.spm.stats.fmri_spec.dir = {statsdir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units='secs';
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = 2; 
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;


% Select files (NIFTI/Analyze)
switch data_type
  case 'img'
    files = spm_select('ExtFPList', datadir, ['^',imagefilter,'I.*\.img']);
    nscans       = size(files,1);            
  case 'nii'
    nifti_file = spm_select('ExtFPList', datadir, ['^',imagefilter,'I.*\.nii']);
    V       = spm_vol(deblank(nifti_file));
    if length(V) == 4
      nframes = V.private.dat.dim(4);
      files = spm_select('ExtFPList', datadir, ['^',imagefilter,'I.*\.nii'],1:nframes);
    else
      files = nifti_file;    
    end
    nscans = size(files, 1);
    clear nifti_file V nframes;
end

matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans = {};

for filecnt = (NUMTRUNC(1)+1):(nscans-NUMTRUNC(2))
  nthfile = filecnt - NUMTRUNC(1);
  matlabbatch{1}.spm.stats.fmri_spec.sess(1).scans{nthfile,1} = ...
    deblank(files(filecnt,:)); 
end

matlabbatch{1}.spm.stats.fmri_spec.sess(1).multi{1} ='';

% Update regressors
load(fullfile(statsdir, 'task_design.mat'));
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg=reg_file;
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf=128;
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];

matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mthresh = 0.4;%% old valure is 0.8

% Update SPM.mat directory
matlabbatch{2}.spm.stats.fmri_est.spmmat{1} = fullfile(statsdir,'SPM.mat'); 

% Update contrasts
matlabbatch{3}.spm.stats.con.spmmat{1}= fullfile(statsdir,'SPM.mat');
%build_contrasts(matlabbatch{1}.spm.stats.fmri_spec.sess);

for i = 1:length(contrastNames)
  if (i <= numTContrasts)
    matlabbatch{3}.spm.stats.con.consess{i}.tcon.name   = contrastNames{i};
    matlabbatch{3}.spm.stats.con.consess{i}.tcon.convec = contrastVecs{i};
  elseif (i > numTContrasts)
    matlabbatch{3}.spm.stats.con.consess{i}.fcon.name = contrastNames{i};
    for j = 1:length(contrastVecs{i}(:,1))
      matlabbatch{3}.spm.stats.con.consess{i}.fcon.convec{j} = ...
	  contrastVecs{i}(j,:);
    end
  end
end

save batch_stats matlabbatch;

% Initialize the batch system
spm_jobman('initcfg');
%delete(get(0,'Children'));
% Run analysis
spm_jobman('run', './batch_stats.mat');
delete(get(0,'Children'));
close all;

end