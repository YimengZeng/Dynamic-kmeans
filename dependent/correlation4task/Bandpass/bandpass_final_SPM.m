function bandpass_final_SPM(TR,fl,fh,dataloc,input_prefix,image_type)

current_dir = pwd;
cd(dataloc);
files = dir([input_prefix '*.' image_type]);

if ~isempty(files) 
  flength = length(files);
  if flength == 1
    nifti4Dto3D (pwd, input_prefix);
    files = dir([input_prefix '*.' image_type]);
    flength = length(files);
  end
  P = cell(flength,1);
  for i = 1:flength
    % Create a array of image file names
    P{i} = files(i).name; 
  end
else
  sprintf('Could not find image files in %s \n', dataloc);
  exit;
end
P = char(P);
V = spm_vol(P);
data = spm_read_vols(V);

% Sampling Frequency
Fs = 1/TR;              

% Set bandpass filter parameters
% Center Frequency
Fc = 0.5*(fh + fl);
% Filter Order
No = floor(Fs * 2/Fc); 

disp('Filtering ........................................................');
[M,N,S,T] = size(data);
data_tmp = zeros(size(data));
% FIR filter Coefficients
B = getFiltCoeffs(zeros(1,T),Fs,fl,fh,0,No);
A = 1;

for s = 1:S
  % Extract time series of nth slice
  X = data(:,:,s,:);              
  Y = reshape(X,M*N,T)';
  % Filter the data
  Ya = filtfilt(B,A,Y);             
  Yb = reshape(Ya',M,N,T);
  data_tmp(:,:,s,:) = Yb;
end

disp('Writing filtered images ..........................................');
for k = 1:flength
  V(k).fname=['filtered', V(k).fname];
  % Write data volume
  spm_write_vol(V(k),squeeze(data_tmp(:,:,:,k)));
end

disp('Done');
cd(current_dir);

end
