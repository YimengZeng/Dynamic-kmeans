% Written by Yimeng Zeng 5/7/2020
% E-mail:enoch29067@126.com
% qinlab.BNU

% dynamic fcuntional connectivity analysis with sliding window and k-means clustering
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

%% running k-means clustering with different parameter set 

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


% state based analysis of their temporal and spatial properties
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
%% similarity analysis
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
    
%% two sample t test with fdr correction
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


%% similarity for integration and segregation

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
 % Cohort 1 98 107; Cohort 2 76 109
 window_num = 76;
 
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

window_num = 109;
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