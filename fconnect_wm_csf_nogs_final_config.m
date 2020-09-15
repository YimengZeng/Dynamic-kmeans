%fconnect_wm_csf_nogs_final('fconnect_wm_csf_nogs_final_config.m')

% Parameter configuration for functional connectivity
% -------------------------------------------------------------------------

% Please specify the data type
paralist.data_type = 'nii';

% Please specify the raw data server
paralist.raw_server =  '';

% Please specify the data structure
% If it is a year structure, non_year_folder = [''];
% If it is a non-year structure, specify the folder containing subjects
paralist.non_year_dir = [''];

% Please specify the path for participant folder
paralist.participant_folder ='';
% Please specify the path for subject list file OR a cell array
%paralist.subjectlist = 'subjectslist_new_51-100.txt';

 load('\CBDA_database.mat'); 
 for i = 1:length(CBDA_database)
  paralist.subjectlist{i,1}=CBDA_database{i,1};
 end

% Please specify the session name
paralist.sess_folder = 'resting';

% Please specify the preprocessed output folder
paralist.preprocessed_folder = 'smoothed_spm8';

% Please specify the prefix to the preprocessed images (pipeline)
paralist.imagefilter = 'swcar';

% Please specify the TR of your data (in seconds)
paralist.TR = 2;

% Please specify the option of bandpass filtering.
% Set to 1 to bandpass filter, 0 to skip.
paralist.bandpass_on = 1;     

% Please specify bandpass filter parameters 
% If not bandpassing (i.e. bandpass_on = 0), then these values are ignored.
% Lower frequency bound for filtering (in Hz)
paralist.fl = 0.008;
% Upper frequency bound for filtering (in Hz)
paralist.fh = 0.1;

% Please specify the ROI folders
paralist.ROI_dir = '';

% Please specify the ROI list (full file name with extensions)
paralist.roi_list = 'list.txt';

% Please specify the number of truncated images from the beginning and end
% (unit in SCANS not seconds, a two element vector, 1st slot for the beginning, 
% and 2nd slot for the end, 0 means no truncation)
paralist.NUMTRUNC   = [0 0];

% =========================================================================
