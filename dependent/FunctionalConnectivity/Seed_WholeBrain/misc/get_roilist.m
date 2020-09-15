function ROI_list = get_roilist(ROIs)

if iscell(ROIs) % already a cell array
  ROI_list = ROIs;
elseif ischar(ROIs) % is string, check if valid directory
  if isdir(ROIs)
    [files,subdir] = spm_select('List',ROIs,'.*\.mat$');
    if (isempty(files)) error('Folder contains no ROIs'); end
    for i=1:size(files,1)
      ROI_list{i} = strcat(ROIs, '/', files(i,:));
    end
  end
else
    error('Please enter a valid ROI cell array or folder');
end

end
