clear all; close all; clc
subj = '39'

% load data
addpath('/media/SSanDra/Pattern_Separation/analysis_pipeline_final')
cd ( ['/media/SSanDra/Pattern_Separation/subj_' subj ])
load(['behavior_subj' subj '.mat'])
load(['trial_data_subj_onset_' subj '_select_chan_3.mat'])
load(['baseline_info_wavelet_32num_20logdb_3hz_350hz_notched_artifact_reject_subj_' subj '_select_chan_3.mat'])
train_images   = cellstr(train_images)
test_images    = cellstr(test_images)
corr_resp_log_vec = testing_behav_matrix(:,3)== testing_behav_matrix(:,5);


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


% study-sim6
pair_2b  = testing_behav_matrix(testing_behav_matrix(:,2) == .5  & corr_resp_log_vec ==1);
pair_2a  = randperm(length(train_images),length(pair_2b))

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
% spectral analysis

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
cd /media/SSanDra/Pattern_Separation/subj_39/figures/idea3d_encode_retrieve_optimal_univariate_measure
all_bands = [freq>1 & freq<4; freq>4 & freq<8 ; freq>8 & freq<15  ; freq>1 & freq<15 ; freq>16&freq<32 ;freq>32];
band_name = {'delta 1-4' 'theta 4-8' 'alpha 8-15' '1-15' 'beta 16-32' 'gamma >32'}; 
band_num  = 1;
maxlag    = 200; %length(200:800)
% stim timing         
stim_onset   = 501;
HC_start     = stim_onset+110; % old was 110
HC_end       = stim_onset+900; % old was 330

for chan = 22:chan_counter
        chan_name = chan_label{chan};

        for band = 1:length(band_name)
            band_idx = all_bands(band,:);
            
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
            
            
            %  output corr
            output_corr = [mean(corr_val_all_trials_pair_1a) mean(corr_val_all_trials_pair_2a)];
            
            
            % plot
            figure('visible','off');
            plot([ 1 2 ], output_corr, 'r*', 'MarkerSize', 14)
            title(['Output: ' chan_name ' ' band_name{band}])
            ylabel('HC pairwase correl', 'FontSize', 12,'FontWeight', 'bold' )
            line([0 3], [0 0], 'Color', 'k')
            set(gca, 'XTick', [1:2], 'XTickLabel', {'stud-rep', 'rand stud-new'})
            set(gca, 'FontSize', 13, 'FontWeight', 'bold')
            saveas(gcf, [band_name{band} '_' num2str(maxlag) '_' chan_name '.jpg'])
            
            
            
        end
end
