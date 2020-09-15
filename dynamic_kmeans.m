% Written by Yimeng Zeng 5/7/2020
% E-mail:enoch29067@126.com
% qinlab.BNU

load('MGZ_seed_6roi_ts.mat');
seed_ts=data;
seed_name=roi_names;
clear data roi_names
load('MGZ_target_20roi_ts.mat');
tar_ts=data;
clear data
%---------------------------------------------------------------------------
% set order parameter
total_time_series = 235; 
bla_order = 5;
cma_order = 6;

% set sliding window  parameter
window_length = 40;
movement_step= 1;

window_number= (total_time_series-window_length)/movement_step+1;

for n = 1:42
    for i=1:window_number    
        initial_point = 1+(i-1)*movement_step;
        end_point = (i-1)*movement_step+window_length;
            for j = 1:20
                corr_bla{1,n}(i,j)=corr(seed_ts{1,n}(bla_order,initial_point:end_point)',tar_ts{1,n}(j,initial_point:end_point)','type','pearson');
                corr_cma{1,n}(i,j)=corr(seed_ts{1,n}(cma_order,initial_point:end_point)',tar_ts{1,n}(j,initial_point:end_point)','type','pearson');
            end
    end
end
% computing for static fc as comparison
for i=1:42
    bla_avg(i,:)=corr_bla{1,i};
    cma_avg(i,:)=corr_cma{1,i};
end

bla_avg_g=mean(bla_avg);
cma_avg_g=mean(cma_avg);
% -----------------------------------------------------------------------
% compute k-means 
% combine cma bla in subject level
% combine all subject data into one group
counter=1;

for i=1:21
    for j=1:window_number
       data(counter,:)=[corr_bla{1,i}(j,:),corr_cma{1,i}(j,:)];
       counter=counter+1;
    end
end

idx_v1=kmeans(data,10,'Display','final','MaxIter',500,'Replicates',1000);


% -------------------------------------------------------------------------------------
%test for best parameter set 

cluster_number_set=2;
clust = zeros(size(data,1),cluster_number_set);

for i=1:cluster_number_set
clust(:,i) = kmeans(data,i,'Display','final','MaxIter',500,'Replicates',1000);
end
% test for cluster quality
eva1=evalclusters(data,clust,'CalinskiHarabasz');
eva2=evalclusters(data,clust,'silhouette');
% plot cluster quality
yyaxis right
plot(eva1.CriterionValues)
yyaxis left
plot(eva2.CriterionValues)

plotyy([1:15],eva1.CriterionValues,[1:15],eva2.CriterionValues)

% get fc state according to their labels
cluster_number=2;
clt=clust(:,cluster_number)';
for i=1:cluster_number
    l{i}=find(clt(:)==i);
end

for i=1:cluster_number
final{i}= mean(data(l{i},:))';
end

 inte_fc = data(l{1},:);
 segre_fc = data(l{2},:);
% inte_fc = data(l{2},:);
% segre_fc = data(l{1},:);


% average fc for bla and cma
avermean(inte_fc(:,1:20))
z_inte = 1/2*log((1+inte_fc)./(1-inte_fc));
z_segre = 1/2*log((1+segre_fc)./(1-segre_fc));

avg_inte(:,1) = mean(z_inte(:,1:20),2);
avg_inte(:,2) = mean(z_inte(:,21:40),2);
avg_segre(:,1) = mean(z_segre(:,1:20),2);
avg_segre(:,2) = mean(z_segre(:,21:40),2);

[h,p,ci,stats] = ttest2(avg_inte(:,1),avg_segre(:,1));
[h,p,ci,stats] = ttest2(avg_inte(:,2),avg_segre(:,2));
[h,p,ci,stats] = ttest(avg_inte(:,1),avg_inte(:,2));
[h,p,ci,stats] = ttest(avg_segre(:,1),avg_segre(:,2));


 MA(1)=nanmean(avg_inte(:,1));
 MA(2)=nanmean(avg_segre(:,1));
 MA(3)=nanmean(avg_inte(:,2));
 MA(4)=nanmean(avg_segre(:,2));
 EA(1)=std(avg_inte(:,1));
 EA(2)=std(avg_segre(:,1));
 EA(3)=std(avg_inte(:,2));
 EA(4)=std(avg_segre(:,2));

   s_location=[1,1.5,2.5,3];
   for s=1:4
      bar(s_location(s),MA(s))
      hold on
      errorbar(s_location(s),MA(s),EA(s))
      hold on
    end



[A,B,r,U,V] = canoncorr(inte_fc(:,1:20),inte_fc(:,21:40)); % not significant

for i=1:20
   % [h,p,ci,stats] = ttest2(inte_fc(:,i),inte_fc(:,i+20));
    [h,p,ci,stats] = ttest2(segre_fc(:,i),segre_fc(:,i+20));
    t_stat(i)= stats.tstat;
    pvalues(i,1) = p;
end
pvalues_corr=mafdr(pvalues,'BHFDR',true); 

for i=1:20
    [h,p,ci,stats] = ttest2((inte_fc(:,i)-inte_fc(:,i+20)),(segre_fc(:,i)-segre_fc(:,i+20)));
        t_stat(i)= stats.tstat;
    pvalues(i) = p;
end
pvalues_corr=mafdr(pvalues,'BHFDR',true); 

% similarity for integration and segregation

for i=1:20
    [r_value1,p_value1]=corr(z_inte(:,i),z_inte(:,i+20),'Type','Pearson');
    [r_value2,p_value2]=corr(z_segre(:,i),z_segre(:,i+20),'Type','Pearson');
    inte_corr(1:2,i)=[r_value1,p_value1];
    segre_corr(1:2,i)=[r_value2,p_value2];
end
  
zz_inte = 1/2*log((1+inte_corr(1,:))./(1-inte_corr(1,:)));
zz_segre = 1/2*log((1+segre_corr(1,:))./(1-segre_corr(1,:)));

[h,p,ci,stats] = ttest(zz_inte,zz_segre)


 MA(1)=nanmean(zz_inte);
 MA(2)=nanmean(zz_segre);
 EA(1)=std(zz_inte)/sqrt(20);
 EA(2)=std(zz_segre)/sqrt(20);

    s_location=[1,1.5];
    for s=1:2
      bar(s_location(s),MA(s))
      hold on
      errorbar(s_location(s),MA(s),EA(s))
      hold on
    end


%% running time lagged cross correlation (CTLL)

a1= mean(z_inte(:,1:20),2);
a2= mean(z_inte(:,21:40),2);
b1= mean(z_segre(:,1:20),2);
b2= mean(z_segre(:,21:40),2);



 lag= 20;
 % CBDA 98 107; MGZ 41 64
 window_num = 41;
 
for j = 1:window_num
    a11 = a1(1+(40)*(j-1):(40*j));
    a22 = a2(1+(40)*(j-1):(40*j));
    for i = 1:lag
        CTLL_inte_before(i)=corr(a22((1+lag-i+1):end),a11(1:(end-lag+i-1)),'Type','Pearson');
        CTll_inte_zero=corr(a22,a11,'Type','Pearson');
        CTLL_inte_after(i)=corr(a11((1+i):end),a22(1:(end-i)),'Type','Pearson');
    end
        CTLL_result = [CTLL_inte_before,CTll_inte_zero,CTLL_inte_after];
        CTLL_window_result_a(j,:) = CTLL_result;
        clear CTLL_inte_before  CTll_inte_zero CTLL_inte_after CTLL_result a11 a22
end

window_num = 64;
for j = 1:window_num
    b11 = b1(1+(40)*(j-1):(40*j));
    b22 = b2(1+(40)*(j-1):(40*j));
    for i = 1:lag
        CTLL_inte_before(i)=corr(b22((1+lag-i+1):end),b11(1:(end-lag+i-1)),'Type','Pearson');
        CTll_inte_zero=corr(b22,b11,'Type','Pearson');
        CTLL_inte_after(i)=corr(b11((1+i):end),b22(1:(end-i)),'Type','Pearson');
    end
        CTLL_result = [CTLL_inte_before,CTll_inte_zero,CTLL_inte_after];
        CTLL_window_result_b(j,:) = CTLL_result;
        clear CTLL_inte_before  CTll_inte_zero CTLL_inte_after CTLL_result b11 b22
end


plot(CTLL_result)
set(gca,'XTickLabel',[-150:150]);
imagesc(CTLL_window_result_a)
imagesc(CTLL_window_result_b)
mean_CTLL(1,:)=mean(CTLL_window_result_a);
mean_CTLL(2,:)=mean(CTLL_window_result_b);

imagesc(mean_CTLL)
imagesc(mean_CTLL(2,:))

plot(mean_CTLL(1,:))
hold on
plot(mean_CTLL(2,:))

%% average WTLCC and t-test
x=[-20:20];
shadedErrorBar(x,mean(CTLL_window_result_b,1),std(CTLL_window_result_b)/sqrt(64),'lineProps','b');
hold on
shadedErrorBar(x,mean(CTLL_window_result_a,1),std(CTLL_window_result_a)/sqrt(41),'lineProps','r');

CTLL_b_z = 1/2*log((1+CTLL_window_result_b)./(1-CTLL_window_result_b));
CTLL_a_z = 1/2*log((1+CTLL_window_result_a)./(1-CTLL_window_result_a));
 
for i=1:41
    [h,p,ci,stats] = ttest(CTLL_b_z(:,i));
    t_result(2,i) = p;
    clear g p ci stats
end
 
pvalues_corr(1,:)=mafdr(t_result(1,:),'BHFDR',true); 


%% normalize scr data and exclude subject by their state bias

    for i=1:21
        final_result_n{1,i}=final_result{1,i}/std(final_result{1,i});
    end
  
  % square root
    for i=1:21
       final_result_n{1,i}=sign(final_result{1,i}).*sqrt(abs(final_result{1,i}));
     end
    
 % plot scr   
    for i=1:20
        subplot(4,5,i)
        plot(final_result_n{1,i});
        ylim([-4,4])
        % plot(final_result_n{1,i});
    end
    
 % functional difference between state and BLA CMA   
 inte=mean(data(find(clt(:)==2),:));
 inte_diff=inte(1:20)-inte(21:40);
 segre=mean(data(find(clt(:)==1),:));
 segre_diff=segre(1:20)-segre(21:40);
 
 for
    
    
    
     exclude_list=[3,5,8]; % for left
     exclude_list=[1,3,5]; % for bilateral 95%
     exclude_list=[1,2,3,5,8]; % for bilateral 90% 
     exclude_list=[1,3,4,6,7,9,11,15,16,17,19,20];
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
      
  clear data_y data_c

% state occupancy for bias default >95% occupancy of subject be removed)
    
 for i=1:length(subj_state)
     s1_occ(i,1)=length(find(subj_state{1,i}(:)==1))/196;
     s2_occ(i,1)=length(find(subj_state{1,i}(:)==2))/196;
 end
        
 for i=1:42
     sub_evolution(i,:)=subj_state{1,i};
 end
 
    imagesc(sub_evolution);
    
 for i=1:196
     sub_mean(1,i)= length(find(sub_evolution(:,i)==1))/42;
     sub_mean(2,i)= length(find(sub_evolution(:,i)==2))/42; 
 end
 plot(sub_mean(1,:));
 hold on 
 plot(sub_mean(2,:));

 % mean occupancy is not good
    plot(s1_occ);
    hold on
    plot(s2_occ);
    
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
    
 % for correlation t-test , be caution to z-score before
 
    a3_z = 1/2*log((1+a3)./(1-a3));
    a4_z = 1/2*log((1+a4)./(1-a4));
    [h,p,ci,stats] = ttest2(a3_z,a4_z(6:end));
    
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

%% group level svm

% cell to matrix
for i=1:length(window_scr_n) 
      data_y(1+(i-1)*window_number:i*window_number,1) = window_scr_n{1,i}';
end

 for i=1:length(subj_state)
      data_c(1+(i-1)*window_number:i*window_number,1)=subj_state{1,i};
 end
 
 % delete roi timeseries vector
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



% svm model
model_s1=fitrsvm(x_s1,y_s1,'Standardize',true,'KFold',10);
model_s1=fitrsvm(x_s1,y_s1,'Standardize',true);
model_s1=fitrsvm(x_s2,y_s2);
model_s1_test=fitrsvm(x_s1,y_s1,'Standardize',true);
model_s1=fitrsvm(data_xx,data_y);

model_s1.ConvergenceInfo.Converged
model_s1.NumIterations


% MSE function
kfoldLoss(model_s1,'Mode','individual')

y_pre=predict(model_s1,data_xx);
% correlation
corr(y_pre, data_y,'Type','Pearson')

% permutation test for svm
  for i=1:10
      rng(i)
      order=randperm(length(y_s1));
      data_y_ran=y_s1(order);
      model_s1_ran=fitrsvm(x_s1,data_y_ran,'Standardize',true);
      predict_y=predict(model_s1_ran,x_s1);
      permu_test(i,1)=corr(predict_y, data_y_ran,'Type','Pearson');
      clear order data_y_ran model_s1_ran predict_y
  end


%% lasso model elstic-net regression
rng(1)
[B1,FitInfo1]= lasso(x_s1,y_s1,'Alpha',0.75,'CV',10);
[B2,FitInfo2]= lasso(x_s2,y_s2,'Alpha',0.75,'CV',10);
[B,FitInfo]= lasso(data_xx,data_y,'Alpha',0.75,'CV',10);
lassoPlot(B1,FitInfo1,'PlotType','CV');
lassoPlot(B2,FitInfo2,'PlotType','CV');
[c1,p1] = corr(x_s1*B1(:,FitInfo1.Index1SE)+FitInfo1.Intercept(1,FitInfo1.Index1SE),y_s1,'Type','Pearson')
[c2,p2] = corr(x_s2*B2(:,FitInfo2.Index1SE)+FitInfo2.Intercept(1,FitInfo2.Index1SE),y_s2,'Type','Pearson')
corr(data_xx*B(:,FitInfo.Index1SE)+FitInfo.Intercept(1,FitInfo.Index1SE),data_y,'Type','Pearson')

B1_sort=sort(abs(B1(:,FitInfo1.Index1SE)),'descend');
thres1=B1_sort((0.2*40),1);
B11=abs(B1(:,FitInfo1.Index1SE));
B11(find(B11(:)>=thres1))=0;

% for bla cma in integration separately

B1_sort=sort(abs(B1(1:20,FitInfo1.Index1SE)),'descend');
thres1=B1_sort((0.4*20),1);
B11=abs(B1(:,FitInfo1.Index1SE));
B11(find(B11>=thres1))=0;


B2_sort=sort(abs(B2(:,FitInfo2.Index1SE)),'descend');
thres2=B2_sort((0.4*20),1);
B22=abs(B2(:,FitInfo2.Index1SE));
B22(find(B22(:)>=thres2))=0;

for i=1:10
 %  order(1,i)=find(B11(:)==B1_sort(i));
   order(2,i)=find(B22(:)==B2_sort(i));
end

find(B22(:)~=0)

inte_segre_diff(1:2,:) = [B11';B22'];
order1=find(inte_segre_diff(1,:)~=0);
order2=find(inte_segre_diff(2,:)~=0);
order3=order1(find(order1(:)==order2(:)));

(B11(order3)+B22(order3))/2

inte_segre_diff(1,order1)

% select strongest edges
system_list=[4,7,3,2,4];
% B1_sort=sort(abs(B1(:,FitInfo1.Index1SE)),'descend');
B2_sort=sort(abs(B2(:,FitInfo2.Index1SE)),'descend');
% MAA=B1_sort(1:10);
MAA=B2_sort(1:10);
ss=[0.5:0.5:5];
  for s=1:10
      bar(ss(s),MAA(s))
      hold on
  end
B222=abs(B2(:,FitInfo2.Index1SE));
find(B222(:)>=MAA(10))

thres1=B1_sort((0.2*40),1);
B11=B1(:,FitInfo1.Index1SE);


B111(1,:)=B11(1:20);
B111(2,:)=B11(21:40);

abs(B111(1,1:4))
    
B111(1,find(B111(1,:)~=0))=1;
B111(2,find(B111(2,:)~=0))=1;
B11(find(B11(:)<thres1))=0;


% permutation test for lasso
  for i=1:100
      rng(i)
      order=randperm(length(y_s1));
      data_y_ran=y_s1(order);
      [B,FitInfo]= lasso(x_s1,data_y_ran,'CV',10);
      permu_test(i,1)=corr(x_s1*B(:,FitInfo.Index1SE)+FitInfo.Intercept(1,FitInfo.Index1SE),data_y_ran,'Type','Pearson');
      clear order data_y_ran B FitInfo 
  end
  a=sort(mean(permu_test,2),'descend');
  
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



% plot linear regression
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