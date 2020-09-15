function task_design_FC(subjfile, GLOBAL)
% GLOBAL : 0 when only the roi timeseries is used as a regressor
% GLOBAL : 0 when only the roi timeseries, global timeseries and 
%            movement covariates are all used as a regressor

sess_name = 'ROI_FC';
if(GLOBAL == 0)
%   reg_file  = [roi_ts(subjcnt,:)]; %filename containing regressors
  reg_names = {'ROI_ts'}; % regressor names, ordered according regressor file
                          % structure
	reg_vec   = [1];
else

  reg_names = {'ROI_ts' 'SPMxGX_ts' 'movement_x' 'movement_y' ...
	       'movement_z' 'pitch' 'roll' 'yaw'};
  reg_vec = [1 0 0 0 0 0 0 0];
  
%   reg_names = {'ROI_ts' 'movement_x' 'movement_y' ...
%     'movement_z' 'pitch' 'roll' 'yaw'};
%   reg_vec = [1 0 0 0 0 0 0];
%   

end

reg_file{1} = subjfile;
save task_design.mat sess_name reg_file reg_names reg_vec

