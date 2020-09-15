% Written by Yimeng Zeng 5/7/2020
% E-mail:enoch29067@126.com
% qinlab.BNU

%% normalize scr data and exclude subject by their state bias

    for i=1:21
        final_result_n{1,i}=final_result{1,i}/std(final_result{1,i});
    end
  
      exclude_list=[4,7,11,15,16]; % final version 7/9/2020
      
%% average scr data accoding to window length
for n = 1:21
    for i=1:window_number  
        initial_point = 1+(i-1)*movement_step;
        end_point = (i-1)*movement_step+window_length;
      % window_scr{1,n}(1,i)=median(final_result{1,n}(1,initial_point:end_point));
        window_scr_n{1,n}(1,i)=mean(final_result_n{1,n}(1,initial_point:end_point));
    end
end   
       window_scr_n(exclude_list)=[];
% separate state label into individual level

    for i=1:42
        subj_state{1,i}=clt(1,1+(i-1)*window_number:i*window_number);
    end
     
      subj_state(exclude_list)=[];
 
% occupancy calculation  
  for i=1:length(subj_state)
     s1_occ(i,1)=length(find(subj_state{1,i}(:)==1))/196;
     s2_occ(i,1)=length(find(subj_state{1,i}(:)==2))/196;
  end

% subject level mean and pair-t test
 clear s1_mean s2_mean
    for i=1:length(subj_state)
        s1_mean(i,1) = mean(window_scr_n{1,i}(1,find(subj_state{1,i}(:)==1)));
        s2_mean(i,1) = mean(window_scr_n{1,i}(1,find(subj_state{1,i}(:)==2)));
    end
   
    
    [h,p,ci,stats] = ttest(s1_mean,s2_mean)
    
% group level mean and independent t test
for i=1:length(window_scr_n) 
      data_y(1+(i-1)*window_number:i*window_number,1) = window_scr_n{1,i}';
end

for i=1:length(subj_state)
     data_c(1+(i-1)*window_number:i*window_number,1)=subj_state{1,i};
end
    s1_mean_group=data_y(find(data_c(:)==1),1);
    s2_mean_group=data_y(find(data_c(:)==2),1);
    
    [h,p,ci,stats] = ttest2(s1_mean_group,s2_mean_group)
    
% bar plot of scr level
    
    MA(1)=nanmean(s1_mean_group);
    MA(2)=nanmean(s2_mean_group);
    EA(1)=std(s1_mean_group)/sqrt(1120);
    EA(2)=std(s2_mean_group)/sqrt(2096);
    
    MA(1)=nanmean(s1_mean);
    MA(2)=nanmean(s2_mean);
    EA(1)=std(s1_mean)/sqrt(15);
    EA(2)=std(s2_mean)/sqrt(15);
    ss=[1,1.5];
    for s=1:2
      bar(ss(s),MA(s))
      hold on
      errorbar(ss(s),MA(s),EA(s))
      hold on
    end
   
    title('el occupancy')
    
%% Elastic net regression 
    % cell to matrix
for i=1:length(window_scr_n) 
      data_y(1+(i-1)*window_number:i*window_number,1) = window_scr_n{1,i}';
end

 for i=1:length(subj_state)
      data_c(1+(i-1)*window_number:i*window_number,1)=subj_state{1,i};
 end
 
 for i=1:length(subj_state)
    data_x{1,i}=data(1+(i-1)*window_number:i*window_number,:);
end
    data_x(exclude_list)=[];
    
for i=1:length(data_x)
    data_xx(1+(i-1)*window_number:i*window_number,:)=data_x{1,i};
end

x_s1=data_xx(find(data_c(:)==1),:);
y_s1=data_y(find(data_c(:)==1),:);
x_s2=data_xx(find(data_c(:)==2),:);
y_s2=data_y(find(data_c(:)==2),:);


[B1,FitInfo1]= lasso(x_s1,y_s1,'Alpha',0.75,'CV',10);
[B2,FitInfo2]= lasso(x_s2,y_s2,'Alpha',0.75,'CV',10);
[B,FitInfo]= lasso(data_xx,data_y,'Alpha',0.75,'CV',10);
lassoPlot(B1,FitInfo1,'PlotType','CV');
lassoPlot(B2,FitInfo2,'PlotType','CV');
[c1,p1] = corr(x_s1*B1(:,FitInfo1.Index1SE)+FitInfo1.Intercept(1,FitInfo1.Index1SE),y_s1,'Type','Pearson')
[c2,p2] = corr(x_s2*B2(:,FitInfo2.Index1SE)+FitInfo2.Intercept(1,FitInfo2.Index1SE),y_s2,'Type','Pearson')
corr(data_xx*B(:,FitInfo.Index1SE)+FitInfo.Intercept(1,FitInfo.Index1SE),data_y,'Type','Pearson')


  

  %% permutation test 
  for i=1:5000
      rng(i)
      order=randperm(length(data_y));
      order1=order(1:1120);
      order2=order(1121:end);
      data_x_1=data_xx(order1,:);
      data_x_2=data_xx(order2,:);
      [B1,FitInfo1]= lasso(data_x_1,data_y(order1),'Alpha',0.75,'CV',10);
      [B2,FitInfo2]= lasso(data_x_2,data_y(order2),'Alpha',0.75,'CV',10);
      permu_test(i,1)=corr(data_x_1*B1(:,FitInfo1.Index1SE)+FitInfo1.Intercept(1,FitInfo1.Index1SE),data_y(order1),'Type','Pearson');
      permu_test(i,2)=corr(data_x_2*B2(:,FitInfo2.Index1SE)+FitInfo2.Intercept(1,FitInfo2.Index1SE),data_y(order2),'Type','Pearson');
      clear order order1 order2 data_y_ran B1 FitInfo1 B2 FitInfo2 
  end

% bootstrap analysis 
  for i=1:1000
      order3=randsample(length(x_s1),500,0);
      order4=randsample(length(x_s2),500,0);
      [B3,FitInfo3]= lasso(x_s1(order3,:),y_s1(order3),'Alpha',0.75,'CV',10);
      [B4,FitInfo4]= lasso(x_s2(order4,:),y_s2(order4),'Alpha',0.75,'CV',10);
      bootstrap_test(i,1)=corr(x_s1(order3,:)*B3(:,FitInfo3.Index1SE)+FitInfo3.Intercept(1,FitInfo3.Index1SE),y_s1(order3),'Type','Pearson');
      bootstrap_test(i,2)=corr(x_s2(order4,:)*B4(:,FitInfo4.Index1SE)+FitInfo4.Intercept(1,FitInfo4.Index1SE),y_s2(order4),'Type','Pearson');
      clear order3 order4 B3 B4 FitInfo3 FitInfo4
  end
  
  
% for correlation t-test , be caution to z-score before
 
 a3=bootstrap_test(:,1);
 a4=bootstrap_test(:,2);
 a3_z = 1/2*log((1+a3)./(1-a3));
 a4_z = 1/2*log((1+a4)./(1-a4));
 [h,p,ci,stats] = ttest2(a3_z,a4_z(6:end));
 
%% plot linear regression
y=x_s1*B(:,FitInfo.Index1SE)+FitInfo.Intercept(1,FitInfo.Index1SE);
y=x_s1*B1(:,FitInfo1.Index1SE)+FitInfo1.Intercept(1,FitInfo1.Index1SE);
y=x_s2*B2(:,FitInfo2.Index1SE)+FitInfo2.Intercept(1,FitInfo2.Index1SE);
x=data_y_ran;
x=y_s1;
x=y_s2;
xfit = min(x):0.01:max(x); 

for ii = 1:1%size(y,2)
    num = strcat('11',num2str(ii));
     subplot(num); a = plot(x,y(:,ii),'o'); set(a,'MarkerEdgeColor','r','MarkerFaceColor','r');
    [p,s] = polyfit(x,y(:,ii),1); 
    [yfit,dy] = polyconf(p,xfit,s,'predopt','curve','simopt','on'); 
    line(xfit,yfit,'color','k','LineWidth',1)
    hold on
    plot(xfit,yfit-dy,'k:')
    plot(xfit,yfit+dy,'k:')
    p = []; s = [];
end