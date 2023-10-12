% Loop thru chan
stim_onset = 301
strt_time = 200
end_time  = 800
% calc theta power vales in all conds: 4-6 hz .3-.8sec
for chan = 1:chan_counter
    
    cond1_4(chan,:) = [nanmean(nanmean(mn_acrs_trials_cond1(freq<6,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond2( freq<6,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond3( freq<6,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond4( freq<6,stim_onset+strt_time:stim_onset+end_time,chan),1))];
    
end

% chans with cond specificity and percentage
cond_spec_log_vec = zeros(chan_counter,1);
for chan = 1:chan_counter  
  PS_theta_perc_chan(chan) =  cond1_4(chan,3) / sum(cond1_4(chan,:));
  if cond1_4(chan,3) > cond1_4(chan,[1:2 4])
  cond_spec_log_vec (chan) = 1;
  end
end



% highest to lowest theta disc 
[chan_order idx]       = sort(PS_theta_perc_chan, 'descend')
all_channel_disc_order = chan_label(idx)'
chan_theta_and_spec = [idx' cond_spec_log_vec(idx)]
theta_and_spec_chans = chan_theta_and_spec(:,1).*chan_theta_and_spec(:,2)
theta_and_spec_chans_final = theta_and_spec_chans(find(theta_and_spec_chans))
chan_name_discrim_spec = chan_label(theta_and_spec_chans_final)';

%% plot 
% make dir
path_name='/mnt/yassamri/iEEG/sandra'
figures_dir = [path_name '/subj_' num2str(subj) '/figures/theta_selective_chans']
if ~exist(figures_dir, 'dir')
mkdir (figures_dir)
end
cd(figures_dir)

%%
for chan =94%34:length(theta_and_spec_chans_final)

figure
    


suptitle(chan_label{theta_and_spec_chans_final(chan)}(4:end-3))
set(gca, 'FontSize',5)

subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond1(:,:,theta_and_spec_chans_final(chan)))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['repeat: ' num2str(sum(trial_log_1_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize',9)

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond2(:,:,theta_and_spec_chans_final(chan)))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['Pattern Comp - incorr lure: ' num2str(sum(trial_log_2_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize',9)

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond3(:,:,theta_and_spec_chans_final(chan)))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['Pattern Sep - corr lure: ' num2str(sum(trial_log_3_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize', 9)

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond4(:,:,theta_and_spec_chans_final(chan)))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['new: ' num2str(sum(trial_log_4_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize',9)

%saveas(gcf, ['chan' num2str(chan) '.jpg']) 

end