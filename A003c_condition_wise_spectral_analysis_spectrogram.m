% clear all
% close all
% clc

%% indicate variables

subj = '39'
exp_type ='tuning_correct' % {'study_test' 'tuning' 'tuning_correct' 'tuning_incorrect' 'indoor_outdoor'}
freq_analysis = 'spectrogram'

% load info
cd (['/mnt/yassamri/Sandra/subj_' subj ])
if strcmp('spectrogram',freq_analysis)
    load(['baseline_info_spectrogram_notched_artifact_reject_subj_' num2str(subj) '.mat'])
else
    load(['baseline_info_wavelet_32num_20logdb_3hz_300hz_notched_artifact_reject_subj_' subj '.mat'])
end
load(['trial_data_subj' subj '.mat'])
load(['behavior_subj' subj '.mat'])

%% organize data into condition

study_trial_data = trial_data(1:nTrials_study,:,:);
test_trial_data = trial_data(nTrials_study+1:end,:,:);

training_behav_matrix_label
testing_behav_matrix_label

repeat_log_vec       = testing_behav_matrix(:,2)== 0;
new_log_vec          = testing_behav_matrix(:,2)== 0.5;
all_lure_log_vec     = testing_behav_matrix(:,2)>=1;
corr_resp_log_vec    = testing_behav_matrix(:,3)== testing_behav_matrix(:,5);


%diff from repeats --> similar to repeats 5-->1
lure_1_log_vec     =testing_behav_matrix(:,2)== 1;
lure_2_log_vec     =testing_behav_matrix(:,2)== 2;
lure_4_log_vec     =testing_behav_matrix(:,2)== 4;
lure_5_log_vec     =testing_behav_matrix(:,2)== 5;

% indoor/outdoor
indoor_log_vec  = training_behav_matrix(:,3)== 1;
outdoor_log_vec = training_behav_matrix(:,3)== 2;

if strcmp('study_test',exp_type)
    cond1 = study_trial_data;
    cond2 = test_trial_data;
    cond3= [];
    cond4= [];
    cond5= [];
    cond6= [];
elseif strcmp('tuning',exp_type)
    cond1= test_trial_data(repeat_log_vec,:,:);
    cond2= test_trial_data(lure_1_log_vec,:,:);
    cond3= test_trial_data(lure_2_log_vec,:,:);
    cond4= test_trial_data(lure_4_log_vec,:,:);
    cond5= test_trial_data(lure_5_log_vec,:,:);
    cond6= test_trial_data(new_log_vec,:,:);
elseif strcmp('tuning_correct',exp_type)
    cond1= test_trial_data(repeat_log_vec&corr_resp_log_vec ==1,:,:);
    cond2= test_trial_data(lure_1_log_vec&corr_resp_log_vec ==1,:,:);
    cond3= test_trial_data(lure_2_log_vec&corr_resp_log_vec ==1,:,:);
    cond4= test_trial_data(lure_4_log_vec&corr_resp_log_vec ==1,:,:);
    cond5= test_trial_data(lure_5_log_vec&corr_resp_log_vec ==1,:,:);
    cond6= test_trial_data(new_log_vec,:,:);
elseif strcmp('tuning_incorrect',exp_type)
    cond1= test_trial_data(repeat_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond2= test_trial_data(lure_1_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond3= test_trial_data(lure_2_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond4= test_trial_data(lure_4_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond5= test_trial_data(lure_5_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond6= test_trial_data(new_log_vec==1&corr_resp_log_vec ==0,:,:);
elseif strcmp('indoor_outdoor',exp_type)
    cond1= study_trial_data(indoor_log_vec,:,:);
    cond2= study_trial_data(outdoor_log_vec,:,:);
    cond3= [];
    cond4= [];
    cond5= [];
    cond6= [];
end

trial_lengths = [size(cond1,1) size(cond2,1) size(cond3,1) size(cond4,1) size(cond5,1) size(cond6,1)];

% exemplar matrix
[~,f,t,ps] = spectrogram(cond1(1, :,1),win,ovrlp,F,fs);

% initialize variables
norm_freq_acrs_chan_cond_1 = nan(size(ps,1), size(ps,2), size(cond1,1), chan_counter);
norm_freq_acrs_chan_cond_2 = nan(size(ps,1), size(ps,2), size(cond2,1), chan_counter);
norm_freq_acrs_chan_cond_3 = nan(size(ps,1), size(ps,2), size(cond3,1), chan_counter);
norm_freq_acrs_chan_cond_4 = nan(size(ps,1), size(ps,2), size(cond4,1), chan_counter);
norm_freq_acrs_chan_cond_5 = nan(size(ps,1), size(ps,2), size(cond5,1), chan_counter);
norm_freq_acrs_chan_cond_6 = nan(size(ps,1), size(ps,2), size(cond6,1), chan_counter);

trial_log_1_artifact = cell(chan_counter,1);
trial_log_2_artifact = cell(chan_counter,1);
trial_log_3_artifact = cell(chan_counter,1);
trial_log_4_artifact = cell(chan_counter,1);
trial_log_5_artifact = cell(chan_counter,1);
trial_log_6_artifact = cell(chan_counter,1);
%% split trials into conds
tic
for chan =1:chan_counter
    
    % first condition
    parfor trial = 1:trial_lengths(1)
        trial_temp = cond1(trial, :,chan);
        % if trail has artifact = nan
        idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
        idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_1_artifact_temp(trial)= sum(total_idx);
        %  spectrogram
        [~,f,t,ps] = spectrogram(trial_temp,win,ovrlp,F,fs);
        ps = 20*log10(ps);
        
        if sum(total_idx)<1 % if trial is artifact free
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(ps));
            for a = 1:size(ps,1) %loop thru freq
                for b = 1:size(ps,2)%loop thru timepoints
                    norm_freq(a, b) = (ps(a,b)-chan_powr_mn_spectrogram(a,chan))/chan_powr_std_spectrogram(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(ps));
        end
        norm_freq_acrs_chan_cond_1(:, :, trial, chan) =  norm_freq;
    end
    trial_log_1_artifact{chan}=  trial_log_1_artifact_temp';
    
    
    % second condition
    parfor trial = 1:trial_lengths(2)
        trial_temp = cond2(trial, :,chan);
        % if trail has artifact = nan
        idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
        idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_2_artifact_temp(trial)= sum(total_idx);
        %  spectrogram
        [~,f,t,ps] = spectrogram(trial_temp,win,ovrlp,F,fs);
        ps = 20*log10(ps);
        if sum(total_idx)<1 % if trial is artifact free
            
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(ps));
            for a = 1:size(ps,1) %loop thru freq
                for b = 1:size(ps,2)%loop thru timepoints
                    norm_freq(a, b) = (ps(a,b)-chan_powr_mn_spectrogram(a,chan))/chan_powr_std_spectrogram(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(ps));
        end
        norm_freq_acrs_chan_cond_2(:, :, trial, chan) =  norm_freq;
    end
    trial_log_2_artifact{chan}=  trial_log_2_artifact_temp';
    
    % third condition
    parfor trial = 1:trial_lengths(3)
        trial_temp = cond3(trial, :,chan);
        % if trail has artifact = nan
        idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
        idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_3_artifact_temp(trial)= sum(total_idx);
        
        %  spectrogram
        [~,f,t,ps] = spectrogram(trial_temp,win,ovrlp,F,fs);
        ps = 20*log10(ps);
        
        if sum(total_idx)<1 % if trial is artifact free
            
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(ps));
            for a = 1:size(ps,1) %loop thru freq
                for b = 1:size(ps,2)%loop thru timepoints
                    norm_freq(a, b) = (ps(a,b)-chan_powr_mn_spectrogram(a,chan))/chan_powr_std_spectrogram(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(ps));
        end
        norm_freq_acrs_chan_cond_3(:, :, trial, chan) =  norm_freq;
    end
    trial_log_3_artifact{chan}=  trial_log_3_artifact_temp';
    
    % fourth condition
    parfor trial = 1:trial_lengths(4)
        trial_temp = cond4(trial, :,chan);
        % if trail has artifact = nan
        idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
        idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_4_artifact_temp(trial)= sum(total_idx);
        
        %  spectrogram
        [~,f,t,ps] = spectrogram(trial_temp,win,ovrlp,F,fs);
        ps = 20*log10(ps);
        
        if sum(total_idx)<1 % if trial is artifact free
            
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(ps));
            for a = 1:size(ps,1) %loop thru freq
                for b = 1:size(ps,2)%loop thru timepoints
                    norm_freq(a, b) = (ps(a,b)-chan_powr_mn_spectrogram(a,chan))/chan_powr_std_spectrogram(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(ps));
        end
        norm_freq_acrs_chan_cond_4(:, :, trial, chan) =  norm_freq;
    end
    trial_log_4_artifact{chan}=  trial_log_4_artifact_temp';
    
    % fifth condition
    parfor trial = 1:trial_lengths(5)
        trial_temp = cond5(trial, :,chan);
        % if trail has artifact = nan
        idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
        idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_5_artifact_temp(trial)= sum(total_idx);
        
        %  spectrogram
        [~,f,t,ps] = spectrogram(trial_temp,win,ovrlp,F,fs);
        ps = 20*log10(ps);
        
        if sum(total_idx)<1 % if trial is artifact free
            
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(ps));
            for a = 1:size(ps,1) %loop thru freq
                for b = 1:size(ps,2)%loop thru timepoints
                    norm_freq(a, b) = (ps(a,b)-chan_powr_mn_spectrogram(a,chan))/chan_powr_std_spectrogram(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(ps));
        end
        
        norm_freq_acrs_chan_cond_5(:, :, trial, chan) =  norm_freq;
    end
    trial_log_5_artifact{chan}=  trial_log_5_artifact_temp';
    
    % sixth condition
    parfor trial = 1:trial_lengths(6)
        trial_temp = cond6(trial, :,chan);
        % if trail has artifact = nan
        idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
        idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
        
        total_idx = idx_above_thresh+idx_below_thresh;
        trial_log_6_artifact_temp(trial)= sum(total_idx);
        
        %  spectrogram
        [~,f,t,ps] = spectrogram(trial_temp,win,ovrlp,F,fs);
        ps = 20*log10(ps);
        
        if sum(total_idx)<1 % if trial is artifact free
            
            % normalize z-score relative to 1-hour session
            norm_freq = zeros(size(ps));
            for a = 1:size(ps,1) %loop thru freq
                for b = 1:size(ps,2)%loop thru timepoints
                    norm_freq(a, b) = (ps(a,b)-chan_powr_mn_spectrogram(a,chan))/chan_powr_std_spectrogram(a, chan);
                end
            end
            
        else
            norm_freq = nan(size(ps));
        end
        norm_freq_acrs_chan_cond_6(:, :, trial, chan) =  norm_freq;
    end
    trial_log_6_artifact{chan}=  trial_log_6_artifact_temp';
    
    
end
toc
%% num of clean trials
for a = 1:size(trial_log_1_artifact,1)
    chan_trials_log_1 (a)= sum(trial_log_1_artifact{a}==0);
end
total_tiral_log_1 = sum(chan_trials_log_1);

for a = 1:size(trial_log_2_artifact,1)
    chan_trials_log_2 (a)= sum(trial_log_2_artifact{a}==0);
    
end
total_tiral_log_2 = sum(chan_trials_log_2);

for a = 1:size(trial_log_3_artifact,1)
    chan_trials_log_3 (a)= sum(trial_log_3_artifact{a}==0);
end
total_tiral_log_3 = sum(chan_trials_log_3);

for a = 1:size(trial_log_4_artifact,1)
    chan_trials_log_4 (a)= sum(trial_log_4_artifact{a}==0);
    
end
total_tiral_log_4 = sum(chan_trials_log_4);

for a = 1:size(trial_log_5_artifact,1)
    chan_trials_log_5 (a)= sum(trial_log_5_artifact{a}==0);
end
total_tiral_log_5 = sum(chan_trials_log_5);


for a = 1:size(trial_log_6_artifact,1)
    chan_trials_log_6 (a)= sum(trial_log_6_artifact{a}==0);
    
end
total_tiral_log_6 = sum(chan_trials_log_6);

%% average across trials within chans


% initialize matrix for all frequencies
mn_acrs_trials_cond1 = nan(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2), chan_counter);
mn_acrs_trials_cond2 = nan(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2), chan_counter);
mn_acrs_trials_cond3 = nan(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2), chan_counter);
mn_acrs_trials_cond4 = nan(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2), chan_counter);
mn_acrs_trials_cond5 = nan(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2), chan_counter);
mn_acrs_trials_cond6 = nan(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2), chan_counter);

% average across trials within chans
for elec = 1:size(norm_freq_acrs_chan_cond_1,4) %loop thru chans
    mn_acrs_trials_cond1 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_1(:, :, :, elec), 3);
    mn_acrs_trials_cond2 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_2(:, :, :, elec), 3);
    mn_acrs_trials_cond3 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_3(:, :, :, elec), 3);
    mn_acrs_trials_cond4 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_4(:, :, :, elec), 3);
    mn_acrs_trials_cond5 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_5(:, :, :, elec), 3);
    mn_acrs_trials_cond6 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_6(:, :, :, elec), 3);
end

% mn across right and left HC

right_HC_cond1 = nanmean(mn_acrs_trials_cond1(:,:, right_HC), 3);
right_HC_cond2 = nanmean(mn_acrs_trials_cond2(:,:, right_HC), 3);
right_HC_cond3 = nanmean(mn_acrs_trials_cond3(:,:, right_HC), 3);
right_HC_cond4 = nanmean(mn_acrs_trials_cond4(:,:, right_HC), 3);
right_HC_cond5 = nanmean(mn_acrs_trials_cond5(:,:, right_HC), 3);
right_HC_cond6 = nanmean(mn_acrs_trials_cond6(:,:, right_HC), 3);

left_HC_cond1 = nanmean(mn_acrs_trials_cond1(:,:, left_HC), 3);
left_HC_cond2 = nanmean(mn_acrs_trials_cond2(:,:, left_HC), 3);
left_HC_cond3 = nanmean(mn_acrs_trials_cond3(:,:, left_HC), 3);
left_HC_cond4 = nanmean(mn_acrs_trials_cond4(:,:, left_HC), 3);
left_HC_cond5 = nanmean(mn_acrs_trials_cond5(:,:, left_HC), 3);
left_HC_cond6 = nanmean(mn_acrs_trials_cond6(:,:, left_HC), 3);

% mean acrs chans
indiv_freq_chans_cond1 = nanmean(mn_acrs_trials_cond1,3);
indiv_freq_chans_cond2 = nanmean(mn_acrs_trials_cond2,3);
indiv_freq_chans_cond3 = nanmean(mn_acrs_trials_cond3,3);
indiv_freq_chans_cond4 = nanmean(mn_acrs_trials_cond4,3);
indiv_freq_chans_cond5 = nanmean(mn_acrs_trials_cond5,3);
indiv_freq_chans_cond6 = nanmean(mn_acrs_trials_cond6,3);


%% plot across all chans
cd(['/mnt/yassamri/Sandra/subj_' num2str(subj) '/figures/spectrogram/' exp_type])

figure
tickmarks = 1:20:length(freq);
pre_stim  = .5;
post_stim = 2;
mn = -.8
mx = .8
time = linspace(-pre_stim,post_stim,length(t));

% plot first condition
subplot (6,1,1)
imagesc(time,1:length(freq), indiv_freq_chans_cond1)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')

subplot (6,1,2)
imagesc(time, 1:length(freq), indiv_freq_chans_cond2)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim++++')

subplot (6,1,3)
imagesc(time, 1:length(freq), indiv_freq_chans_cond3)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim+++')

subplot (6,1,4)
imagesc(time, 1:length(freq), indiv_freq_chans_cond4)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim++')

subplot (6,1,5)
imagesc(time, 1:length(freq), indiv_freq_chans_cond5)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim+')

subplot (6,1,6)
imagesc(time, 1:length(freq), indiv_freq_chans_cond6)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
saveas(gcf, 'all_chans.jpg') 



%% plot right HC
cd(['/mnt/yassamri/Sandra/subj_' num2str(subj) '/figures/' exp_type])

figure
% plot first condition
subplot (6,1,1)
imagesc(time, 1:length(freq), right_HC_cond1)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')

subplot (6,1,2)
imagesc(time, 1:length(freq), right_HC_cond2)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim++++')

subplot (6,1,3)
imagesc(time, 1:length(freq), right_HC_cond3)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim+++')

subplot (6,1,4)
imagesc(time, 1:length(freq), right_HC_cond4)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim++')

subplot (6,1,5)
imagesc(time, 1:length(freq), right_HC_cond5)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim+')

subplot (6,1,6)
imagesc(time, 1:length(freq), right_HC_cond6)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
saveas(gcf, 'right_HC_chans.jpg') 

%% Left HC
figure
% plot first condition
subplot (6,1,1)
imagesc(time, 1:length(freq), left_HC_cond1)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')

subplot (6,1,2)
imagesc(time, 1:length(freq), left_HC_cond2)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim++++')

subplot (6,1,3)
imagesc(time, 1:length(freq), left_HC_cond3)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim+++')

subplot (6,1,4)
imagesc(time, 1:length(freq), left_HC_cond4)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim++')

subplot (6,1,5)
imagesc(time, 1:length(freq), left_HC_cond5)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('sim+')

subplot (6,1,6)
imagesc(time, 1:length(freq), left_HC_cond6)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
saveas(gcf, 'left_HC_chans.jpg') 

%% plot per chan
mx = 1.5
mn = -1.5
for chan = 1:chan_counter
    
    figure
    fig_name = chan_label(chan)
    suptitle(fig_name{1}(4:end))
    % plot first condition
    subplot (6,1,1)
    imagesc(time, 1:length(freq), mn_acrs_trials_cond1(:,:,chan))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar 
    
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('repeat')
    
    subplot (6,1,2)
    imagesc(time, 1:length(freq),  mn_acrs_trials_cond2(:,:,chan))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('sim++++')
    
    subplot (6,1,3)
    imagesc(time, 1:length(freq), mn_acrs_trials_cond3(:,:,chan))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('sim+++')
    
    subplot (6,1,4)
    imagesc(time, 1:length(freq),  mn_acrs_trials_cond4(:,:,chan))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('sim++')
    
    subplot (6,1,5)
    imagesc(time, 1:length(freq),  mn_acrs_trials_cond5(:,:,chan))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('sim+')
    
    subplot (6,1,6)
    imagesc(time, 1:length(freq),  mn_acrs_trials_cond6(:,:,chan))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('new')
    saveas(gcf, [cell2mat(chan_label(chan)) '.jpg']) 
    close all
end

%% plot trials

for chan = 1:chan_counter
    figure
    fig_name = chan_label(chan)
    title([fig_name{1}(4:8) ' repeats'])
    hold on
    yaxis_shift = 300
    for trial = 1:size(cond1,1)
          trial_temp = cond1(trial, :,chan);
        % if trail has artifact = nan
        idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
        idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
        
        total_idx = idx_above_thresh+idx_below_thresh;
       
       if sum(total_idx)<1 % if trial is artifact free
            % normalize z-score relative to 1-hour sess
            a = a+1
            clean_trials(a,:)=trial_temp;
        plot(-.5:1/fs:2, cond1(trial,:,chan)+(yaxis_shift*a))
        end
    end
    saveas(gcf, [cell2mat(chan_label(chan)) 'time.jpg'])
    close all
end
