% window 1
group = {'39' '57' '63' '44' '66' }%
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
cd('/mnt/yassamri/iEEG/sandra/group_data')
st = 500
ed = 800

for a = 1:length(group)
    % load(['MTL_theta_' num2str(st) '_' num2str(ed) '_subj_' subj{a} '.mat'])
    %load(['MTL_theta_allchans_subj_' group{a} '.mat'])
    %load(['temporal_lobe_theta_allchans_subj_' group{a} '.mat'])
    %load(['temporal_lobe_theta_' num2str(st) '_' num2str(ed) '_subj_' group{a} '.mat'])
    %load(['frontal_lobe_theta_' num2str(st) '_' num2str(ed) '_subj_' group{a} '.mat'])
    %load(['frontal_lobe_theta_allchans_subj_' group{a} '.mat'])
    
   % load(['MTL_gamma_allchans_subj_' group{a} '.mat'])
   % load(['temporal_lobe_gamma_allchans_subj_' group{a} '.mat'])
    %load(['frontal_lobe_gamma_allchans_subj_' group{a} '.mat'])

    
   % load(['MTL_gamma_' num2str(st) '_' num2str(ed) '_subj_' group{a} '.mat'])
   % load(['temporal_lobe_gamma_' num2str(st) '_' num2str(ed) '_subj_' group{a} '.mat'])
    load(['frontal_lobe_gamma_' num2str(st) '_' num2str(ed) '_subj_' group{a} '.mat'])
    
all_subj_cond_1(a,:) = all_conds(1,:);
all_subj_cond_2(a,:) = all_conds(2,:);
all_subj_cond_3(a,:) = all_conds(3,:);
all_subj_cond_4(a,:) = all_conds(4,:);

end

figure
hold on
stdshade(all_subj_cond_1,.1,'b',[],[] ,freq_name,win)
%h1 = plot(nanmean(all_subj_cond_1,1), 'b', 'LineWidth', 2)
h1 = plot(conv(nanmean(all_subj_cond_1,1),ones(1,win)/win,'same'),'b', 'LineWidth', 2)

stdshade(all_subj_cond_2,.1,'g',[],[] ,freq_name,win)
%h2 = plot(nanmean(all_subj_cond_2,1), 'g', 'LineWidth', 2)
h2 = plot(conv(nanmean(all_subj_cond_2,1),ones(1,win)/win,'same'),'g', 'LineWidth', 2)

stdshade(all_subj_cond_3,.1,'r',[],[] ,freq_name,win)
%h3 = plot(nanmean(all_subj_cond_3,1), 'r', 'LineWidth', 2)
h3 = plot(conv(nanmean(all_subj_cond_3,1),ones(1,win)/win,'same'),'r', 'LineWidth', 2)

stdshade(all_subj_cond_4,.1,'m',[],[] ,freq_name,win)
%h4 = plot(nanmean(all_subj_cond_4,1), 'm', 'LineWidth', 2)
h4 = plot(conv(nanmean(all_subj_cond_4,1),ones(1,win)/win,'same'),'m', 'LineWidth', 2)

xlabel('time from stim onset')
ylabel('z-score power')
legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')
title([ num2str(st)  ' - ' num2str(ed) ' ms channel selectivity'])
%title(['all temporal lobe chans - theta'])
%title(['all temporal lobe chans - theta'])

set(gca, 'FontSize', 12, 'FontWeight', 'bold')
