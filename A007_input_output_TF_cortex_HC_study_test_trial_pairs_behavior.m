%clear all; close all; clc
subj = '39'

% load data
cd ( ['/media/SSanDra/Pattern_Separation/subj_' subj ])
% load(['behavior_subj' subj '.mat'])
% load(['trial_data_subj_onset_' subj '_select_chan_1.mat'])
% load(['baseline_info_wavelet_32num_20logdb_3hz_350hz_notched_artifact_reject_subj_' subj '_select_chan_1.mat'])
train_images   = cellstr(train_images)
test_images    = cellstr(test_images)
corr_resp_log_vec    = testing_behav_matrix(:,3)== testing_behav_matrix(:,5);


% condition pairs
% study-R

pair_1b        = testing_behav_matrix(testing_behav_matrix(:,2)== 0 & corr_resp_log_vec ==1);
pair_1b_images = test_images(pair_1b)
[IDX1 IDX2]    = ismember(train_images,pair_1b_images);
IDX3    = [find(IDX1) IDX2(IDX1)]
IDX4    = sortrows(IDX3,2)
pair_1a = IDX4(:,1)
pair_1a_images = train_images(pair_1a)
train_images(pair_1a)
test_images(pair_1b)

% study-sim1
pair_2b        = testing_behav_matrix(testing_behav_matrix(:,2)== 1 & corr_resp_log_vec ==0);
pair_2b_images = test_images(pair_2b)
for trial = 1:length(pair_2b_images)
  pair_2a_images{trial,1} =   [pair_2b_images{trial}(1:3) 'a' pair_2b_images{trial}(5:end)]
end
[IDX1 IDX2]    = ismember(train_images,pair_2a_images);
IDX3    = [find(IDX1) IDX2(IDX1)]
IDX4    = sortrows(IDX3,2)
pair_2a = IDX4(:,1)
pair_2a_images = train_images(pair_2a)
train_images(pair_2a)
test_images(pair_2b)

% study-sim2
pair_3b        = testing_behav_matrix(testing_behav_matrix(:,2)== 2 & corr_resp_log_vec ==0);
pair_3b_images = test_images(pair_3b)
for trial = 1:length(pair_3b_images)
  pair_3a_images{trial,1} =   [pair_3b_images{trial}(1:3) 'a' pair_3b_images{trial}(5:end)]
end
[IDX1 IDX2]    = ismember(train_images,pair_3a_images);
IDX3    = [find(IDX1) IDX2(IDX1)]
IDX4    = sortrows(IDX3,2)
pair_3a = IDX4(:,1)
pair_3a_images = train_images(pair_3a)
train_images(pair_3a)
test_images(pair_3b)

% study-sim3
pair_4b        = testing_behav_matrix(testing_behav_matrix(:,2)== 4 & corr_resp_log_vec ==1);
pair_4b_images = test_images(pair_4b)
for trial = 1:length(pair_4b_images)
  pair_4a_images{trial,1} =   [pair_4b_images{trial}(1:3) 'a' pair_4b_images{trial}(5:end)]
end
[IDX1 IDX2]    = ismember(train_images,pair_4a_images);
IDX3    = [find(IDX1) IDX2(IDX1)]
IDX4    = sortrows(IDX3,2)
pair_4a = IDX4(:,1)
pair_4a_images = train_images(pair_4a)
train_images(pair_4a)
test_images(pair_4b)

% study-sim4
pair_5b        = testing_behav_matrix(testing_behav_matrix(:,2)== 5 & corr_resp_log_vec ==1);
pair_5b_images = test_images(pair_5b)
for trial = 1:length(pair_5b_images)
  pair_5a_images{trial,1} =   [pair_5b_images{trial}(1:3) 'a' pair_5b_images{trial}(5:end)]
end
[IDX1 IDX2]    = ismember(train_images,pair_5a_images);
IDX3    = [find(IDX1) IDX2(IDX1)]
IDX4    = sortrows(IDX3,2)
pair_5a = IDX4(:,1)
pair_5a_images = train_images(pair_5a)
train_images(pair_5a)
test_images(pair_5b)


% study-sim6
pair_6b  = testing_behav_matrix(testing_behav_matrix(:,2) == .5  & corr_resp_log_vec ==1);
pair_6a  = randperm(length(train_images),length(pair_6b))

%%

% organize trials: study and test
study_trial_data = trial_data(1:nTrials_study,:,:);
test_trial_data  = trial_data(nTrials_study+1:end,:,:);

cond1 = study_trial_data;
cond2 = test_trial_data;

trial_lengths = [size(cond1,1) size(cond2,1)];
    
% wavelet params
dt = 1/fs;
NumVoices = 32;
a0 = 2^(1/NumVoices);
wavCenterFreq = 6/(2*pi);
minscale = wavCenterFreq/(maxfreq*dt); 
maxscale = wavCenterFreq/(minfreq*dt); 
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale)); 
scales   = a0.^(minscale:maxscale).*dt;
freq     = wavCenterFreq./(fs*scales.*dt);

% exemplar matrix
cwt = cwtft({cond1(1, :,1),dt},...
    'scales',scales,'wavelet','morl');
cwt_power_exemp = 10*log10(abs(cwt.cfs).^2);

% initialize variables
norm_freq_acrs_chan_cond_1 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond1,1), chan_counter);
norm_freq_acrs_chan_cond_2 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond2,1), chan_counter);

trial_log_1_artifact = cell(chan_counter,1);
trial_log_2_artifact = cell(chan_counter,1);
clear trial
%% spectral analysis

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
    
    % wavelet spectrogram
    cwt = cwtft({trial_temp,dt},...
        'scales',scales,'wavelet','morl');
    cwt_power = 10*log10(abs(cwt.cfs).^2);
           
    if sum(total_idx)<1 % if trial is artifact free
        
        % normalize z-score relative to 1-hour session
        norm_freq = zeros(size(cwt_power));
        for a = 1:size(cwt_power,1) %loop thru freq
            for b = 1:size(cwt_power,2)%loop thru timepoints
                norm_freq(a, b) = (cwt_power(a,b)-chan_powr_mn(a,chan))/chan_powr_std(a, chan);
            end
        end
        
    else
        norm_freq = nan(size(cwt_power));
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
    % wavelet spectrogram
    cwt = cwtft({trial_temp,dt},...
        'scales',scales,'wavelet','morl');
    cwt_power = 10*log10(abs(cwt.cfs).^2);
    
    if sum(total_idx)<1 % if trial is artifact free
        
        % normalize z-score relative to 1-hour session
        norm_freq = zeros(size(cwt_power));
        for a = 1:size(cwt_power,1) %loop thru freq
            for b = 1:size(cwt_power,2)%loop thru timepoints
                norm_freq(a, b) = (cwt_power(a,b)-chan_powr_mn(a,chan))/chan_powr_std(a, chan);
            end
        end
        
    else
        norm_freq = nan(size(cwt_power));
    end
   
    norm_freq_acrs_chan_cond_2(:, :, trial, chan) =  norm_freq;
end
trial_log_2_artifact{chan}=  trial_log_2_artifact_temp';
end
toc
%% split data into bands: band_data:  5bands X timepoints X trials X chans
band_num  = 6;
band_data_1 = nan(band_num, size(norm_freq_acrs_chan_cond_1,2), size(norm_freq_acrs_chan_cond_1,3), size(norm_freq_acrs_chan_cond_1,4));
band_data_2 = nan(band_num, size(norm_freq_acrs_chan_cond_2,2), size(norm_freq_acrs_chan_cond_2,3), size(norm_freq_acrs_chan_cond_2,4));

% banded study data
band_data_1 (1,:,:,:) = nanmean(norm_freq_acrs_chan_cond_1(find(freq<4 & freq>1),:,:,:),1);
band_data_1 (2,:,:,:) = nanmean(norm_freq_acrs_chan_cond_1(find(freq<8 & freq>4),:,:,:),1);
band_data_1 (3,:,:,:) = nanmean(norm_freq_acrs_chan_cond_1(find(freq<15 & freq>8),:,:,:),1);
band_data_1 (4,:,:,:) = nanmean(norm_freq_acrs_chan_cond_1(find(freq<32 & freq>16),:,:,:),1);
band_data_1 (5,:,:,:) = nanmean(norm_freq_acrs_chan_cond_1(find(freq<80 & freq>32),:,:,:),1);
band_data_1 (6,:,:,:) = nanmean(norm_freq_acrs_chan_cond_1(find(freq>80),:,:,:),1);

% banded test data
band_data_2 (1,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq<4 & freq>1),:,:,:),1);
band_data_2 (2,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq<8 & freq>4),:,:,:),1);
band_data_2 (3,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq<15 & freq>8),:,:,:),1);
band_data_2 (4,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq<32 & freq>16),:,:,:),1);
band_data_2 (5,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq<80 & freq>32),:,:,:),1);
band_data_2 (6,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq>80),:,:,:),1);

%%
chan_name = 'POLRHH3Ref'
EC = 29 % IR 39
HC = find(ismember(chan_label, chan_name)) % IR 39:42

stim_onset   = 501;
NC_start     = stim_onset+90;
NC_end       = stim_onset+300;
HC_start     = stim_onset+200; % old was 110
HC_end       = stim_onset+900; % old was 330

% pair 1 generate EC and CA3 vectors for study and test in 
good_trials = [];
for trial = 1:length(pair_1a)
    % study
   cond_1a_EC(trial,:)= nanmean(band_data_1(:,NC_start:NC_end, pair_1a(trial), EC),2);
   cond_1a_HC(trial,:)= nanmean(band_data_1(:,HC_start:HC_end, pair_1a(trial), HC),2);
   
   % test
   cond_1b_EC(trial,:)= nanmean(band_data_2(:,NC_start:NC_end, pair_1b(trial), EC),2);
   cond_1b_HC(trial,:)= nanmean(band_data_2(:,HC_start:HC_end, pair_1b(trial), HC),2);

   % find clean trials
  if sum(isnan([cond_1a_EC(trial,:) cond_1a_HC(trial,:) cond_1b_EC(trial,:) cond_1b_HC(trial,:)]))==0
      good_trials = [good_trials trial]
  end
end

% good trials
cond_1a_EC = cond_1a_EC(good_trials,:);
cond_1b_EC = cond_1b_EC(good_trials,:);

cond_1a_HC = cond_1a_HC(good_trials,:);
cond_1b_HC = cond_1b_HC(good_trials,:);

corr_val_all_trials_pair_1a = zeros(size(cond_1a_EC,1),1);
corr_val_all_trials_pair_1b = zeros(size(cond_1a_EC,1),1);

for trial = 1:size(cond_1a_EC,1) 
    % corr of input
    [corr_val sig]= corr([cond_1a_EC(trial,:) ; cond_1b_EC(trial,:)]')
    corr_val_all_trials_pair_1a (trial) = [corr_val(2)]
    % corr of output
    [corr_val sig]= corr([cond_1a_HC(trial,:) ; cond_1b_HC(trial,:)]')
    corr_val_all_trials_pair_1b (trial) = [corr_val(2)]
end

% pair 2
good_trials = [];
for trial = 1:length(pair_2a)
   % study
   cond_2a_EC(trial,:)= nanmean(band_data_1(:,NC_start:NC_end, pair_2a(trial), EC),2);
   cond_2a_HC(trial,:)= nanmean(band_data_1(:,HC_start:HC_end, pair_2a(trial), HC),2);
   
   % test
   cond_2b_EC(trial,:)= nanmean(band_data_2(:,NC_start:NC_end, pair_2b(trial), EC),2);
   cond_2b_HC(trial,:)= nanmean(band_data_2(:,HC_start:HC_end, pair_2b(trial), HC),2);

   % find clean trials
  if sum(isnan([cond_2a_EC(trial,:) cond_2a_HC(trial,:) cond_2b_EC(trial,:) cond_2b_HC(trial,:)]))==0
      good_trials = [good_trials trial]
  end
end

% good trials
cond_2a_EC = cond_2a_EC(good_trials,:);
cond_2b_EC = cond_2b_EC(good_trials,:);

cond_2a_HC = cond_2a_HC(good_trials,:);
cond_2b_HC = cond_2b_HC(good_trials,:);

corr_val_all_trials_pair_2a = zeros(size(cond_2a_EC,1),1);
corr_val_all_trials_pair_2b = zeros(size(cond_2a_EC,1),1);

for trial = 1:size(cond_2a_EC,1) 
    % corr of input
    [corr_val sig]= corr([cond_2a_EC(trial,:) ; cond_2b_EC(trial,:)]')
    corr_val_all_trials_pair_2a (trial) = [corr_val(2)]
    % corr of output
    [corr_val sig]= corr([cond_2a_HC(trial,:) ; cond_2b_HC(trial,:)]')
    corr_val_all_trials_pair_2b (trial) = [corr_val(2)]
end
% pair 3
good_trials = [];
for trial = 1:length(pair_3a)
    % study
   cond_3a_EC(trial,:)= nanmean(band_data_1(:,NC_start:NC_end, pair_3a(trial), EC),2);
   cond_3a_HC(trial,:)= nanmean(band_data_1(:,HC_start:HC_end, pair_3a(trial), HC),2);
   
   % test
   cond_3b_EC(trial,:)= nanmean(band_data_2(:,NC_start:NC_end, pair_3b(trial), EC),2);
   cond_3b_HC(trial,:)= nanmean(band_data_2(:,HC_start:HC_end, pair_3b(trial), HC),2);

   % find clean trials
  if sum(isnan([cond_3a_EC(trial,:) cond_3a_HC(trial,:) cond_3b_EC(trial,:) cond_3b_HC(trial,:)]))==0
      good_trials = [good_trials trial]
  end
end

% good trials
cond_3a_EC = cond_3a_EC(good_trials,:);
cond_3b_EC = cond_3b_EC(good_trials,:);

cond_3a_HC = cond_3a_HC(good_trials,:);
cond_3b_HC = cond_3b_HC(good_trials,:);

corr_val_all_trials_pair_3a = zeros(size(cond_3a_EC,1),1);
corr_val_all_trials_pair_3b = zeros(size(cond_3a_EC,1),1);

for trial = 1:size(cond_3a_EC,1) 
    % corr of input
    [corr_val sig]= corr([cond_3a_EC(trial,:) ; cond_3b_EC(trial,:)]')
    corr_val_all_trials_pair_3a (trial) = [corr_val(2)]
    % corr of output
    [corr_val sig]= corr([cond_3a_HC(trial,:) ; cond_3b_HC(trial,:)]')
    corr_val_all_trials_pair_3b (trial) = [corr_val(2)]
end

% pair 4
good_trials = [];
for trial = 1:length(pair_4a)
    % study
   cond_4a_EC(trial,:)= nanmean(band_data_1(:,NC_start:NC_end, pair_4a(trial), EC),2);
   cond_4a_HC(trial,:)= nanmean(band_data_1(:,HC_start:HC_end, pair_4a(trial), HC),2);
   
   % test
   cond_4b_EC(trial,:)= nanmean(band_data_2(:,NC_start:NC_end, pair_4b(trial), EC),2);
   cond_4b_HC(trial,:)= nanmean(band_data_2(:,HC_start:HC_end, pair_4b(trial), HC),2);

   % find clean trials
  if sum(isnan([cond_4a_EC(trial,:) cond_4a_HC(trial,:) cond_4b_EC(trial,:) cond_4b_HC(trial,:)]))==0
      good_trials = [good_trials trial]
  end
end

% good trials
cond_4a_EC = cond_4a_EC(good_trials,:);
cond_4b_EC = cond_4b_EC(good_trials,:);

cond_4a_HC = cond_4a_HC(good_trials,:);
cond_4b_HC = cond_4b_HC(good_trials,:);

corr_val_all_trials_pair_4a = zeros(size(cond_4a_EC,1),1);
corr_val_all_trials_pair_4b = zeros(size(cond_4a_EC,1),1);

for trial = 1:size(cond_4a_EC,1) 
    % corr of input
    [corr_val sig]= corr([cond_4a_EC(trial,:) ; cond_4b_EC(trial,:)]')
    corr_val_all_trials_pair_4a (trial) = [corr_val(2)]
    % corr of output
    [corr_val sig]= corr([cond_4a_HC(trial,:) ; cond_4b_HC(trial,:)]')
    corr_val_all_trials_pair_4b (trial) = [corr_val(2)]
end

% pair 5
good_trials = [];
for trial = 1:length(pair_5a)
    % study
   cond_5a_EC(trial,:)= nanmean(band_data_1(:,NC_start:NC_end, pair_5a(trial), EC),2);
   cond_5a_HC(trial,:)= nanmean(band_data_1(:,HC_start:HC_end, pair_5a(trial), HC),2);
   
   % test
   cond_5b_EC(trial,:)= nanmean(band_data_2(:,NC_start:NC_end, pair_5b(trial), EC),2);
   cond_5b_HC(trial,:)= nanmean(band_data_2(:,HC_start:HC_end, pair_5b(trial), HC),2);

   % find clean trials
  if sum(isnan([cond_5a_EC(trial,:) cond_5a_HC(trial,:) cond_5b_EC(trial,:) cond_5b_HC(trial,:)]))==0
      good_trials = [good_trials trial]
  end
end

% good trials
cond_5a_EC = cond_5a_EC(good_trials,:);
cond_5b_EC = cond_5b_EC(good_trials,:);

cond_5a_HC = cond_5a_HC(good_trials,:);
cond_5b_HC = cond_5b_HC(good_trials,:);

corr_val_all_trials_pair_5a = zeros(size(cond_5a_EC,1),1);
corr_val_all_trials_pair_5b = zeros(size(cond_5a_EC,1),1);

for trial = 1:size(cond_5a_EC,1) 
    % corr of input
    [corr_val sig]= corr([cond_5a_EC(trial,:) ; cond_5b_EC(trial,:)]')
    corr_val_all_trials_pair_5a (trial) = [corr_val(2)]
    % corr of output
    [corr_val sig]= corr([cond_5a_HC(trial,:) ; cond_5b_HC(trial,:)]')
    corr_val_all_trials_pair_5b (trial) = [corr_val(2)]
end

% pair 6
good_trials = [];
for trial = 1:length(pair_6a)
    % study
   cond_6a_EC(trial,:)= nanmean(band_data_1(:,NC_start:NC_end, pair_6a(trial), EC),2);
   cond_6a_HC(trial,:)= nanmean(band_data_1(:,HC_start:HC_end, pair_6a(trial), HC),2);
   
   % test
   cond_6b_EC(trial,:)= nanmean(band_data_2(:,NC_start:NC_end, pair_6b(trial), EC),2);
   cond_6b_HC(trial,:)= nanmean(band_data_2(:,HC_start:HC_end, pair_6b(trial), HC),2);

   % find clean trials
  if sum(isnan([cond_6a_EC(trial,:) cond_6a_HC(trial,:) cond_6b_EC(trial,:) cond_6b_HC(trial,:)]))==0
      good_trials = [good_trials trial]
  end
end

% good trials
cond_6a_EC = cond_6a_EC(good_trials,:);
cond_6b_EC = cond_6b_EC(good_trials,:);

cond_6a_HC = cond_6a_HC(good_trials,:);
cond_6b_HC = cond_6b_HC(good_trials,:);

corr_val_all_trials_pair_6a = zeros(size(cond_6a_EC,1),1);
corr_val_all_trials_pair_6b = zeros(size(cond_6a_EC,1),1);

for trial = 1:size(cond_6a_EC,1) 
    % corr of input
    [corr_val sig]= corr([cond_6a_EC(trial,:) ; cond_6b_EC(trial,:)]')
    corr_val_all_trials_pair_6a (trial) = [corr_val(2)]
    % corr of output
    [corr_val sig]= corr([cond_6a_HC(trial,:) ; cond_6b_HC(trial,:)]')
    corr_val_all_trials_pair_6b (trial) = [corr_val(2)]
end


% input output corr
input_corr = [mean(corr_val_all_trials_pair_1a) mean(corr_val_all_trials_pair_2a) mean(corr_val_all_trials_pair_3a)...
    mean(corr_val_all_trials_pair_4a) mean(corr_val_all_trials_pair_5a)  mean(corr_val_all_trials_pair_6a) ]

output_corr = [mean(corr_val_all_trials_pair_1b) mean(corr_val_all_trials_pair_2b) mean(corr_val_all_trials_pair_3b)...
    mean(corr_val_all_trials_pair_4b) mean(corr_val_all_trials_pair_5b) mean(corr_val_all_trials_pair_6b) ]
%
figure;
plot([1 2 3 4 5 6], output_corr, 'r*', 'MarkerSize', 14)
title(['Output: ' chan_name])
xlabel('corre of incr diff stim pairs', 'FontSize', 12,'FontWeight', 'bold' )
ylabel('Corresponding HC pairwase correl', 'FontSize', 12,'FontWeight', 'bold' )
ylim([-1 1])
line([0 6], [0 0], 'Color', 'k')
set(gca, 'XTick', [1:6], 'XTickLabel', {'stud-rep', 'stud-sim1', 'stud-sim2', 'stud-sim3', 'stud-sim4', 'rand stud-new'}, 'XTickLabelRotation', 70)
set(gca, 'FontSize', 13, 'FontWeight', 'bold')

% figure;
% plot([1 2 3 4 5 6], input_corr, 'r*', 'MarkerSize', 14)
% title('Input: EC')