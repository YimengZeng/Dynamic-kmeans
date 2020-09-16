% Written by Yimeng Zeng 5/7/2020
% E-mail:enoch29067@126.com
% qinlab.BNU


% list={'MGZ1001','MGZ1004','MGZ1006','MGZ1009','MGZ1011','MGZ1012','MGZ1013','MGZ1015','MGZ1016','MGZ1017','MGZ1018',...
%    'MGZ1019','MGZ1022','MGZ1023','MGZ1024','MGZ1025','MGZ1026','MGZ1028','MGZ1029','MGZ1030','MGZ1031'};



scan_length=8*60*1000;
for i=1:length(all_data)
    end_p=max(find(all_data{1,i}(:,2)~=0));
    start_p=min(find(all_data{1,i}(:,2)~=0));
    signal{1,i}=all_data{1,i}(end_p-scan_length+1:end_p,1);
end

% Down sample 
  for i=1:length(signal)
       for j=1:240
        d_signal{1,i}(1,j)=mean(signal{1,i}(1+2000*(j-1):2000*j,1));
       end
  end

% Detrend  
  for i=1:length(signal)
       Dd_signal{1,i}=detrend(d_signal{1,i});
  end
  
% Bandpass filter with the same range of 0.008-0.1 Hz
  for i=1:length(Dd_signal)
      FDd_signal{1,i}=bandpass(Dd_signal{1,i},[0.008,0.1],0.5);
  end
    
 final_result = FDd_signal;
    



