clear all; close all; clc
subj = '44'

% load data
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
cd ( ['/mnt/yassamri/iEEG/sandra/subj_' subj ])
load(['behavior_subj' subj '.mat'])
load(['trial_data_subj_onset_' subj '_select_chan_3.mat'])
load(['baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_' subj '_select_chan_3.mat'])
train_images          = cellstr(train_images)
test_images           = cellstr(test_images)
corr_resp_log_vec     = testing_behav_matrix(:,3)== testing_behav_matrix(:,5);


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

if subj == '57'
    chan_counter = 70
end
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
% split data into bands: band_data:  5bands X timepoints X trials X chans
% make band index matrix
all_bands = [freq<4 & freq>1 ;freq<8 & freq>4;freq<15 & freq>8;freq<32 & freq>16;freq>32];
band_name = {'delta' 'theta' 'alpha' 'beta' 'gamma'};
band_num  = 1;
maxlag = 200 %length(200:800)
% stim timing         
stim_onset   = 501;
NC_start     = stim_onset+90;
NC_end       = stim_onset+300;
HC_start     = stim_onset+200; % old was 110
HC_end       = stim_onset+800; % old was 330
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/figures/idea3c_univariate_output_correlations_tuning_curve'])
%%
for chan = 1:chan_counter
        chan_name = chan_label{chan};

    for band = 1:length(band_name)       
        band_idx = all_bands(band,:);
        
        % study and test conds
        band_data_1 = nan(band_num, size(norm_freq_acrs_chan_cond_1,2), size(norm_freq_acrs_chan_cond_1,3), size(norm_freq_acrs_chan_cond_1,4));
        band_data_2 = nan(band_num, size(norm_freq_acrs_chan_cond_2,2), size(norm_freq_acrs_chan_cond_2,3), size(norm_freq_acrs_chan_cond_2,4));
        
        % banded study data
        band_data_1 (1,:,:,:) = nanmean(norm_freq_acrs_chan_cond_1(band_idx,:,:,:),1);
        
        % banded test data
        band_data_2 (1,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(band_idx,:,:,:),1);

       
        %% pair 1
        good_trials = [];
        for trial = 1:length(pair_1a)
            % study
            cond_1a(trial,:)= band_data_1(:,HC_start:HC_end, pair_1a(trial), chan);
            
            % test
            cond_1b(trial,:)= band_data_2(:,HC_start:HC_end, pair_1b(trial), chan);
            
            % find clean trials
            if sum(isnan([cond_1a(trial,:) cond_1b(trial,:)]))==0
                good_trials = [good_trials trial];
            end
        end
        
        % good trials
        cond_1a = cond_1a(good_trials,:);
        cond_1b = cond_1b(good_trials,:);

        % corr in chan chan 
        corr_val_all_trials_pair_1a = zeros(size(cond_1a,1),1);        
        for trial = 1:size(cond_1a,1)
            corr_val_all_trials_pair_1a (trial) = unique(max(xcorr(cond_1a(trial,:), cond_1b(trial,:), maxlag,'coeff')));
        end
        
        %% pair 2
         good_trials = [];
        for trial = 1:length(pair_2a)
            % study
            cond_2a(trial,:)= band_data_1(:,HC_start:HC_end, pair_2a(trial), chan);
            
            % test
            cond_2b(trial,:)= band_data_2(:,HC_start:HC_end, pair_2b(trial), chan);
            
            % find clean trials
            if sum(isnan([cond_2a(trial,:) cond_2b(trial,:)]))==0
                good_trials = [good_trials trial];
            end
        end
        
        % good trials
        cond_2a = cond_2a(good_trials,:);
        cond_2b = cond_2b(good_trials,:);

        % corr in chan chan 
        corr_val_all_trials_pair_2a = zeros(size(cond_2a,1),1);        
        for trial = 1:size(cond_2a,1)
            corr_val_all_trials_pair_2a (trial) = unique(max(xcorr(cond_2a(trial,:), cond_2b(trial,:), maxlag,'coeff')));
        end
        
        
        
        %% pair 3
          good_trials = [];
        for trial = 1:length(pair_3a)
            % study
            cond_3a(trial,:)= band_data_1(:,HC_start:HC_end, pair_3a(trial), chan);
            
            % test
            cond_3b(trial,:)= band_data_2(:,HC_start:HC_end, pair_3b(trial), chan);
            
            % find clean trials
            if sum(isnan([cond_3a(trial,:) cond_3b(trial,:)]))==0
                good_trials = [good_trials trial];
            end
        end
        
        % good trials
        cond_3a = cond_3a(good_trials,:);
        cond_3b = cond_3b(good_trials,:);

        % corr in chan chan 
        corr_val_all_trials_pair_3a = zeros(size(cond_3a,1),1);        
        for trial = 1:size(cond_3a,1)
            corr_val_all_trials_pair_3a (trial) = unique(max(xcorr(cond_3a(trial,:), cond_3b(trial,:), maxlag,'coeff')));
        end

        
        %% pair 4
         good_trials = [];
        for trial = 1:length(pair_4a)
            % study
            cond_4a(trial,:)= band_data_1(:,HC_start:HC_end, pair_4a(trial), chan);
            
            % test
            cond_4b(trial,:)= band_data_2(:,HC_start:HC_end, pair_4b(trial), chan);
            
            % find clean trials
            if sum(isnan([cond_4a(trial,:) cond_4b(trial,:)]))==0
                good_trials = [good_trials trial];
            end
        end
        
        % good trials
        cond_4a = cond_4a(good_trials,:);
        cond_4b = cond_4b(good_trials,:);

        % corr in chan chan 
        corr_val_all_trials_pair_4a = zeros(size(cond_4a,1),1);        
        for trial = 1:size(cond_4a,1)
            corr_val_all_trials_pair_4a (trial) = unique(max(xcorr(cond_4a(trial,:), cond_4b(trial,:), maxlag,'coeff')));
        end

        
        %% pair 5
         good_trials = [];
        for trial = 1:length(pair_5a)
            % study
            cond_5a(trial,:)= band_data_1(:,HC_start:HC_end, pair_5a(trial), chan);
            
            % test
            cond_5b(trial,:)= band_data_2(:,HC_start:HC_end, pair_5b(trial), chan);
            
            % find clean trials
            if sum(isnan([cond_5a(trial,:) cond_5b(trial,:)]))==0
                good_trials = [good_trials trial];
            end
        end
        
        % good trials
        cond_5a = cond_5a(good_trials,:);
        cond_5b = cond_5b(good_trials,:);

        % corr in chan chan 
        corr_val_all_trials_pair_5a = zeros(size(cond_5a,1),1);        
        for trial = 1:size(cond_5a,1)
            corr_val_all_trials_pair_5a (trial) = unique(max(xcorr(cond_5a(trial,:), cond_5b(trial,:), maxlag,'coeff')));
        end

        
        %% pair 6
        good_trials = [];
        for trial = 1:length(pair_6a)
            % study
            cond_6a(trial,:)= band_data_1(:,HC_start:HC_end, pair_6a(trial), chan);
            
            % test
            cond_6b(trial,:)= band_data_2(:,HC_start:HC_end, pair_6b(trial), chan);
            
            % find clean trials
            if sum(isnan([cond_6a(trial,:) cond_6b(trial,:)]))==0
                good_trials = [good_trials trial];
            end
        end
        
        % good trials
        cond_6a = cond_6a(good_trials,:);
        cond_6b = cond_6b(good_trials,:);

        % corr in chan chan 
        corr_val_all_trials_pair_6a = zeros(size(cond_6a,1),1);        
        for trial = 1:size(cond_6a,1)
            corr_val_all_trials_pair_6a (trial) = unique(max(xcorr(cond_6a(trial,:), cond_6b(trial,:), maxlag,'coeff')));
        end

        
        %  output corr
        output_corr = [mean(corr_val_all_trials_pair_1a) mean(corr_val_all_trials_pair_2a) mean(corr_val_all_trials_pair_3a)...
            mean(corr_val_all_trials_pair_4a) mean(corr_val_all_trials_pair_5a)  mean(corr_val_all_trials_pair_6a) ];

        
        % plot
        figure('visible','off');
        plot([1 2 3 4 5 6], output_corr, 'r*', 'MarkerSize', 14)
        title(['Output: ' chan_name ' ' band_name{band}])
        ylabel('HC pairwase correl', 'FontSize', 12,'FontWeight', 'bold' )
        ylim([-1 1])
        line([0 6], [0 0], 'Color', 'k')
        set(gca, 'XTick', [1:6], 'XTickLabel', {'stud-rep', 'stud-sim1', 'stud-sim2', 'stud-sim3', 'stud-sim4', 'rand stud-new'}, 'XTickLabelRotation', 70)
        set(gca, 'FontSize', 13, 'FontWeight', 'bold')
        saveas(gcf, [band_name{band} '_' num2str(maxlag) '_' chan_name '.jpg']) 
    end
end
