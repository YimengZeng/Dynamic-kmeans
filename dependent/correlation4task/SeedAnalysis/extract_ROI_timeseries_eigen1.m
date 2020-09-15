%Extract time series of a given ROI by taking the first eigen of all voxels
function [timeseries, roi_name] = extract_ROI_timeseries_eigen1(ROIs, datadir, ...
                                  filter, NUMTRUNC, imagefilter, data_type)

% Add path of marsbar functions 
%addpath /home/fmri/fmrihome/SPM/spm8/toolbox/marsbar;

% Unzip images if necessary
unix(sprintf('gunzip -fq %s', fullfile(datadir, [imagefilter, '*.gz'])));

% Select files
switch data_type
  case 'img'
    files = spm_select('ExtFPList', datadir, ['^',imagefilter,'I.*\.img']);
    nscans       = size(files,1);            
  case 'nii'
    nifti_file = spm_select('ExtFPList', datadir, ['^',imagefilter,'I.*\.nii']);
    V       = spm_vol(deblank(nifti_file(1,:)));
    if length(V.private.dat.dim) == 4
      nframes = V.private.dat.dim(4);
      files = spm_select('ExtFPList', datadir, ['^',imagefilter,'I.*\.nii'],1:nframes);
    else
      files = nifti_file;    
    end
    nscans = size(files, 1);
    clear nifti_file V nframes;
end 
if nscans == 0
  error('No image data found');
end

% Truncate selected files
truncpoint = NUMTRUNC+1;
select_data = files(truncpoint:end,:);

% Get ROIs
ROI_list = ROIs;%get_roilist(ROIs);
nrois = length(ROI_list);
if nrois==0
  error('No ROIs specified')
end

% Get Timeseries for each ROI
timeseries = [];
roi_name = cell(nrois,1);
for j = 1:nrois
  roi_obj = maroi(ROI_list{j});
  roi_name{j} = label(roi_obj);
  roi_data_obj = get_marsy(roi_obj, select_data, 'eigen1');
  roi_ts = summary_data(roi_data_obj);
  % Zero-mean and linear detrend if specified
  if filter 
    roi_ts = roi_ts - mean(roi_ts);
    roi_ts = detrend(roi_ts);
  end
  timeseries = [timeseries; roi_ts'];
end

disp('ROI timeseries extraction - Done');
%unix(sprintf('gzip -fq %s', fullfile(datadir, '*.nii')));

end
