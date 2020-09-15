% -------------- Please specify the individualstats server: resting data, adults and kids, 5 HTT ----------------
%paralist.stats_path     = '/fs/musk2';
paralist.raw_server     = '/home/qinlab/data';
paralist.data_type      = 'nii';

paralist.imagefilter    = 'swcar'; 
paralist.NFRAMES = 232;
paralist.NUMTRUNC = [0,0]; % this parameter is dropped number of volume at beginning and in the end, for task, 4 volume should be dropped refered to liangxiao cc and pnas paper
%NUMTRUNC                = [8,40];%for subject 14-02-17.1_3T2 only

paralist.TR_val = 2;
paralist.bandpass_on = 1;     
paralist.fl = 0.008;
% Upper frequency bound for filtering (in Hz)
paralist.fh = 0.1;

%-white matter and CSF roi files
paralist.wm_csf_roi_file = cell(2,1);
%-white matter and csf rois
paralist.wm_csf_roi_file{1} = '/home/qinlab/SPM/spm8_scripts/correlation4task/white_mask_p08_d1_e1_roi.mat';
paralist.wm_csf_roi_file{2} = '/home/qinlab/SPM/spm8_scripts/correlation4task/csf_mask_p08_d1_e1_roi.mat'; 

% Please specify file name holding subjects to be analyzed
% For one subject list files. For eg.,
paralist.subjlist_file  = {'sublist_TA_36.txt'}; 
 
paralist.ROI_dir        = '/home/qinlab/data/genghaiyang/data/nback/SPM12_results/sub_TA_h18_l18/ppi_networks/ROIs_mat';
paralist.output_folder  = '/home/qinlab/data/genghaiyang/data/nback/SPM12_results/sub_TA_h18_l18/ppi_networks/Results_pearson';
% get roiname:
niidir = dir(paralist.ROI_dir);
niidir = niidir(3:end);
niidir = struct2cell(niidir);
niidir = niidir(1,:);
niidir = niidir';
paralist.roi_list  = niidir;
%fullfile(paralist.roi_nii_folder,'.',niidir);
% ----- Please specify the session name: adults, kids, 5 HTT ----------------
paralist.non_year_dir = [''];
paralist.sess_folder  = 'nback';
paralist.preprocessed_folder = 'smoothed_spm8';