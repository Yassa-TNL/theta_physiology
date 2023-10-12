clear all;close all;clc
T = zeros(1000, 4)
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
pre_stim = 0.5;
fs       = 500;
reg      = 'HC';
exp_type = 'tuning_correct'
[subj_list] = get_subj_list(reg);
if strcmp('encoding',exp_type)
    time_idx = [.3*fs 2.5*fs]
else
    time_idx = [.3*fs 1.5*fs]
end

% exemplar
parent_fn='/mnt/yassamri/iEEG/sandra/group_data/IndivSubjData/';
child_fn='LM_reref/normalization_cond_spec_prestim/tuning_correct_onset/';
cd([parent_fn child_fn subj_list{1}])
load(['subj' subj_list{1} '_' reg '_spectrograms_cond_spec_prestim_cue_responseive_yes.mat'])

% load zmap and adjust size
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/LM_reref/cluster_matrices'])
range = 'deltatheta'   %theta delta deltatheta alpha gamma
load([reg '_onset_' exp_type '_groupcluster_' range '.mat'])
zmapthresh_temp = zeros(size(cond2,1),size(cond2,2));
zmapthresh_temp(:,time_idx(1):time_idx(2)) = logical(zmapthresh);
figure;imagesc(zmapthresh)
cntr=0;
for iSubj =1:4

cd([parent_fn child_fn subj_list{iSubj}])
load(['subj' subj_list{iSubj} '_' reg '_spectrograms_cond_spec_prestim_cue_responseive_yes.mat'])

% reset to baseline
temp3 = squeeze(nanmean(cond3,3)); % avg acros trls in chan -> freqXtimeXchan
temp2 = squeeze(nanmean(cond2,3)); 

reset = nanmean(nanmean(temp3(:,0.3*fs:pre_stim*fs,:),3),2);
cond3 = cond3 - reset;

reset = nanmean(nanmean(temp2(:,0.3*fs:pre_stim*fs,:),3),2);
cond2 = cond2 - reset;
   
for iElec = 1:size(cond3,4)
    for iTrl =1:size(cond3,3)
        % cond3
        temp = cond3(:,:,iTrl,iElec);
        if ~isnan(nanmean(temp(logical(zmapthresh_temp))))
        cntr=cntr+1;
        T(cntr,1) = iSubj;%subj
        T(cntr,2) = iElec;%elec
        T(cntr,3) = 3;   %predictor
        T(cntr,4) = nanmean(temp(logical(zmapthresh_temp)));%depend variable           
        end
        clear tmp
    end
end  
    
for iElec = 1:size(cond2,4)
   
    for iTrl =1:size(cond2,3)
        % cond3
        temp = cond2(:,:,iTrl,iElec);
        if ~isnan(nanmean(temp(logical(zmapthresh_temp))))
        cntr=cntr+1;
        T(cntr,1) = iSubj;
        T(cntr,2) = iElec;
        T(cntr,3) = 2;
        T(cntr,4) = nanmean(temp(logical(zmapthresh_temp)));            
        end
        clear tmp
    end
end  

end

T = T(1:cntr,:);% remove zero values
table = array2table(T, 'VariableNames',{ 'Subject', 'Electrode', 'Condition', 'ThetaPower'})
table.Condition = nominal(table.Condition );
writetable(table,'LME_StatsTable_HC_Theta.csv')
%% generate the model

% general linear model 
lme   = fitlme(table,['ThetaPower ~ 1 + Condition + Electrode +  Electrode:Subject)'] , 'FitMethod', 'REML')

% start
lme   = fitlme(table,['ThetaPower ~ 1 + Condition + Electrode + (1 | Subject)'] , 'FitMethod', 'REML')
lme   = fitlme(table,['ThetaPower ~ 1 + Condition + Electrode + (1 | Subject) + (1 | Electrode:Subject)'], 'FitMethod', 'REML')
                            % 1 = estimate intercept =avg acrs obs irrespc of treatment type
 
                            % (1 | Electrode:Subject): randome effect -
                            % diff offsets for each combination of elec and
                            % subject. the effect of theta powerby each
                            % elec-sub comb has a different baseline value.
                            % 
                            
% random slope
lme   = fitlme(table,['ThetaPower ~ 1 + Condition + Electrode + (1 | Subject) + (-1 + Subject | Condition)'], 'FitMethod', 'REML')
lme   = fitlme(table,['ThetaPower ~ 1 + Condition + Electrode + (1 | Subject) + (-1 + Subject | Condition)+(1 + Electrode | Subject)'], 'FitMethod', 'REML')
% bonferoni correc
%%
lme   = fitlme(table,['ThetaPower ~ 1 + Condition + Electrode + (1 | Subject) + (1 + Electrode | Subject)'], 'FitMethod', 'REML')
% random slope: the magnitude of their change from one condition to the
% next is not going to be exactly the same. 
%% 
%lme   = fitlme(table,['ThetaPower ~ 1 + Condition + (Electrode|Subject)'], 'FitMethod', 'REML')
%lme   = fitlme(table,['ThetaPower ~ 1 + Condition + (1 |Subject) + (1 + Condition|Electrode:Subject)'], 'FitMethod', 'REML')

% condition differences % is the second p-value the degree to which cond 2
% differs from cond1?
%pVal = coefTest(lme)
stat = anova(lme)

%% examine slope and intercept variation across subjects
close all
colors = {'r' 'b' 'g' 'c'}
figure;hold on
for iSubj = 1:4
    cond1_x_vals = ones(1,length(T(T(:,1)==iSubj & T(:,3)==3,4)));
    cond2_x_vals = 2*ones(1,length(T(T(:,1)==iSubj & T(:,3)==2,4)));
    
    cond1_y_vals = T(T(:,1)==iSubj & T(:,3)==3,4)';
    cond2_y_vals = T(T(:,1)==iSubj & T(:,3)==2,4)';
    
    plot(cond1_x_vals,cond1_y_vals, 'o', 'MarkerFaceColor', colors{iSubj}, 'MarkerEdgeColor','none')
    plot(cond2_x_vals,cond2_y_vals, 'o', 'MarkerFaceColor', colors{iSubj}, 'MarkerEdgeColor','none')
    xlim([0 3])
    ylim([-2 2])

% fit line and get slope and intercept
p = polyfit([cond1_x_vals cond2_x_vals],[cond1_y_vals cond2_y_vals],1);
m = p(1); % y = mx + b
b = p(2);
plot([cond1_x_vals cond2_x_vals], m*[cond1_x_vals cond2_x_vals]+b,'color', colors{iSubj}, 'LineWidth', 1);
% running stats on significancy of the slope
stats = regstats([cond1_y_vals cond2_y_vals], [cond1_x_vals cond2_x_vals], 'Linear');
b_2 = stats.beta(1);  % intercept value, the same as b.
m_2 = stats.beta(2);  % slope value, the same as m.
b_pval = stats.tstat.pval(1); % p-value of significancy of estimates intercept b 
if b_pval <0.001
    b_pval = '<0.001';
else
    b_pval = ['=' num2str(b_pval)];
end
m_pval = stats.tstat.pval(2); % p-value of significancy of estimates slope m 
if m_pval <0.001
   m_pval = '<0.001'; 
else
    m_pval = ['=' num2str(m_pval)];
end
stats.tstat.t(2); % t-value associated to slope test statistics

xlim([0 3]); xticks([1 2]);xticklabels({'cond3', 'cond2'});
%y =ylim;
y =  [-2 2]
if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
end
%title (['subj' num2str(iSubj)  ' b = ' num2str(b) ', p' b_pval ', m = ' num2str(m) ', p' m_pval])
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
end
print('-clipboard','-dbitmap')



%% examine slope and intercept variation across electrodes w/ in subject
figure;hold on
for iSubj = 4
    for iElec = 1:2
        if isempty(T(T(:,1)==iSubj & T(:,2)==iElec & T(:,3)==3,4))
            continue
        end
    iCond = 3; 
    cond1_x_vals = ones(1,length(T(T(:,1)==iSubj & T(:,2)==iElec & T(:,3)==iCond,4)));
    cond1_y_vals =T(T(:,1)==iSubj & T(:,2)==iElec & T(:,3)==iCond,4)';
    
    iCond = 2; 
    cond2_x_vals = 2*ones(1,length(T(T(:,1)==iSubj & T(:,2)==iElec & T(:,3)==iCond,4)));
    cond2_y_vals =T(T(:,1)==iSubj & T(:,2)==iElec & T(:,3)==iCond,4)';
    
    plot(cond1_x_vals,cond1_y_vals, 'o', 'MarkerFaceColor', colors{iSubj}, 'MarkerEdgeColor','none')
    plot(cond2_x_vals,cond2_y_vals, 'o', 'MarkerFaceColor', colors{iSubj}, 'MarkerEdgeColor','none')
    xlim([0 3])
    ylim([-2 2])

% fit line and get slope and intercept
p = polyfit([cond1_x_vals cond2_x_vals],[cond1_y_vals cond2_y_vals],1);
m = p(1); % y = mx + b
b = p(2);
plot([cond1_x_vals cond2_x_vals], m*[cond1_x_vals cond2_x_vals]+b,'color', colors{iSubj}, 'LineWidth', 1);
% running stats on significancy of the slope
stats = regstats([cond1_y_vals cond2_y_vals], [cond1_x_vals cond2_x_vals], 'Linear');
b_2 = stats.beta(1);  % intercept value, the same as b.
m_2 = stats.beta(2);  % slope value, the same as m.
b_pval = stats.tstat.pval(1); % p-value of significancy of estimates intercept b 
if b_pval <0.001
    b_pval = '<0.001';
else
    b_pval = ['=' num2str(b_pval)];
end
m_pval = stats.tstat.pval(2); % p-value of significancy of estimates slope m 
if m_pval <0.001
   m_pval = '<0.001'; 
else
    m_pval = ['=' num2str(m_pval)];
end
stats.tstat.t(2); % t-value associated to slope test statistics

xlim([0 3]); xticks([1 2]);xticklabels({'cond3', 'cond2'});
y =ylim;
%y =  [-2 2]
if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
end
%title (['subj' num2str(iSubj)  ' b = ' num2str(b) ', p' b_pval ', m = ' num2str(m) ', p' m_pval])
title (['subj' num2str(iSubj) ])

set(gca, 'FontSize', 12, 'FontWeight', 'bold')
    end
end
print('-clipboard','-dbitmap')

%% matlab examples
load flu

flu2 = stack(flu,2:10,'NewDataVarName','FluRate',...
    'IndVarName','Region');
flu2.Date = nominal(flu2.Date);

lme = fitlme(flu2,'FluRate ~ 1 + Region + (1|Date)')
% fixed effect for region, and a random intercept that varies by data
% there are 9 regions, it was estimate the slope for each region