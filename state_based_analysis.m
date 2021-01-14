% Written by Yimeng Zeng 5/7/2020
% E-mail:enoch29067@126.com
% qinlab.BNU

%% normalize scr data and exclude subject by their state bias

    for i=1:37
        final_result_n{1,i}=final_result{1,i}/std(final_result{1,i});
    end
  
      exclude_list=[7,11,15,16,22,23,24,25,33];  % details in supplementary figurg.S5
      
%% average scr data accoding to window length
for n = 1:length(final_result_n)
    for i=1:window_number  
        initial_point = 1+(i-1)*movement_step;
        end_point = (i-1)*movement_step+window_length;
      % window_scr{1,n}(1,i)=median(final_result{1,n}(1,initial_point:end_point));
        window_scr_n{1,n}(1,i)=mean(final_result_n{1,n}(1,initial_point:end_point));
    end
end   
       window_scr_n(exclude_list)=[];
% separate state label into individual level

    for i=1:length(final_result_n)
        subj_state{1,i}=clt(1,1+(i-1)*window_number:i*window_number);
    end
     
      subj_state(exclude_list)=[];



%% Preparing FCs as input and corresponding SCL as response variable
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

%% optimal parameter set (alpha\lambda) for Elastic-net regression model 
con_vec = linspace(0,1,21);
rng(1)
for i=1:20
    [B1,FitInfo1]= lasso(x_s1,y_s1,'Alpha',con_vec(i+1),'CV',10);
    [B2,FitInfo2]= lasso(x_s2,y_s2,'Alpha',con_vec(i+1),'CV',10);
    inte_elastic{i,1}=B1;
    inte_elastic{i,2}=FitInfo1;
    segre_elastic{i,1}=B2;
    segre_elastic{i,2}=FitInfo2;
    clear B1 B2 FitInfo1 FitInfo2
end

% Choosing the parameter set with minimum standard error
for i=1:20
   inte_comp(i) = inte_elastic{i, 2}.MSE(inte_elastic{i, 2}.IndexMinMSE);
   segre_comp(i) = segre_elastic{i, 2}.MSE(segre_elastic{i, 2}.IndexMinMSE);
end

    inte_alpha = find(inte_comp(:) == min(inte_comp));
    segre_alpha = find(segre_comp(:) == min(segre_comp));
    
    B1 = inte_elastic{inte_alpha,1};
    FitInfo1 = inte_elastic{inte_alpha,2}; 
    B2 = segre_elastic{segre_alpha,1};
    FitInfo2 = segre_elastic{segre_alpha,2}; 
    
    plot(inte_comp)
    hold on
    plot(segre_comp)
    
  % fig.s7  
    lassoPlot(B1,FitInfo1,'PlotType','CV');
    lassoPlot(B2,FitInfo2,'PlotType','CV');
    
  % calcuating the correlation between predicted and real scl values
    [c1,p1] = corr(x_s1*B1(:,FitInfo1.Index1SE)+FitInfo1.Intercept(1,FitInfo1.Index1SE),y_s1,'Type','Pearson')
    [c2,p2] = corr(x_s2*B2(:,FitInfo2.Index1SE)+FitInfo2.Intercept(1,FitInfo2.Index1SE),y_s2,'Type','Pearson')



%% 2nd permutation test for group effect (integration\segregation group)
  
  for i=1:5000
      order1=randperm(length(data_y));
      order2 = order1(1:length(y_s1));
      order3 = order1((length(y_s1)+1):end);
      [B1,FitInfo1]= lasso(data_xx(order2,:),data_y(order2),'Alpha',0.65,'CV',10);
      [B2,FitInfo2]= lasso(data_xx(order3,:),data_y(order3),'Alpha',0.35,'CV',10);

      permu_test_2(i,1)=corr(data_xx(order2,:)*B1(:,FitInfo1.Index1SE)+FitInfo1.Intercept(1,FitInfo1.Index1SE),data_y(order2),'Type','Pearson');
      permu_test_2(i,2)=corr(data_xx(order3,:)*B2(:,FitInfo2.Index1SE)+FitInfo2.Intercept(1,FitInfo2.Index1SE),data_y(order3),'Type','Pearson');
      clear order1 order2 order3 B1 FitInfo1 B2 FitInfo2 
      i
  end
      
      permu_test_2(:,3)=sort(permu_test_2(:,1),'descend');
      permu_test_2(:,4)=sort(permu_test_2(:,2),'descend');
  

  
%% bootstrapping for testing difference between integration and segregation

  for i=1:5000
      order3=randsample(length(x_s1),1000,1);
      order4=randsample(length(x_s2),1000,1);
      [B3,FitInfo3]= lasso(x_s1(order3,:),y_s1(order3),'Alpha',0.65,'CV',10);
      [B4,FitInfo4]= lasso(x_s2(order4,:),y_s2(order4),'Alpha',0.35,'CV',10);
      bootstrap_test(i,3)=corr(x_s1(order3,:)*B3(:,FitInfo3.Index1SE)+FitInfo3.Intercept(1,FitInfo3.Index1SE),y_s1(order3),'Type','Pearson');
      bootstrap_test(i,4)=corr(x_s2(order4,:)*B4(:,FitInfo4.Index1SE)+FitInfo4.Intercept(1,FitInfo4.Index1SE),y_s2(order4),'Type','Pearson');
      clear order3 order4 B3 B4 FitInfo3 FitInfo4
      i
  end

     
    [h,p,ci,stats] =  ttest2(bootstrap_test(:,1),bootstrap_test(:,2));
    bootstrap_test(:,1) = sort(bootstrap_test(:,3),'descend');
    bootstrap_test(:,2) = sort(bootstrap_test(:,4),'descend');
  
% over