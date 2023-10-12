% input-outptu TF: cortex-HC
clear all; close all;clc
subj = '39'

% load data
cd ( ['/media/SSanDra/Pattern_Separation/subj_' subj ])
load(['behavior_subj' subj '.mat'])
load(['trial_data_subj_onset_' subj '_select_chan_3.mat'])
load(['baseline_info_wavelet_32num_20logdb_3hz_350hz_notched_artifact_reject_subj_' subj '_select_chan_3.mat'])

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
scales = a0.^(minscale:maxscale).*dt;
freq = wavCenterFreq./(fs*scales.*dt);

% exemplar matrix
cwt = cwtft({cond1(1, :,1),dt},...
    'scales',scales,'wavelet','morl');
cwt_power_exemp = 10*log10(abs(cwt.cfs).^2);

% initialize variables
norm_freq_acrs_chan_cond_1 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond1,1), chan_counter);
norm_freq_acrs_chan_cond_2 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond2,1), chan_counter);

trial_log_1_artifact = cell(chan_counter,1);
trial_log_2_artifact = cell(chan_counter,1);

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

%% split data into bands: band_data:  5bands X timepoints X trials X chans
band_num  = 6;
band_data_2 = nan(band_num, size(norm_freq_acrs_chan_cond_2,2), size(norm_freq_acrs_chan_cond_2,3), size(norm_freq_acrs_chan_cond_2,4));

band_data_2 (1,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq<4 & freq>1),:,:,:),1);
band_data_2 (2,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq<8 & freq>4),:,:,:),1);
band_data_2 (3,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq<15 & freq>8),:,:,:),1);
band_data_2 (4,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq<32 & freq>16),:,:,:),1);
band_data_2 (5,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq>32 & freq<80),:,:,:),1);
band_data_2 (6,:,:,:) = nanmean(norm_freq_acrs_chan_cond_2(find(freq>80),:,:,:),1);

%% cortex vs. HC groupings
for HC_CHAN  = 2
HC_chans     = zeros(1,length(chan_label));
cortex_chans = ones(1,length(chan_label));
hc_chan_idx  = [39:42 51 52 86 87 95 96] % 39
HC_chans([hc_chan_idx])     = 1;
HC_chans     = find(HC_chans)
cortex_chans([hc_chan_idx]) = 0;
cortex_chans = find(cortex_chans)

 cortex_chans = 30 % use this line to only plot 1 cortical elec

stim_onset   = 501;
cortex_start = stim_onset+90;
cortex_end   = stim_onset+300;
HC_start     = stim_onset+110;
HC_end       = stim_onset+330;

% test phase
% loop thru each trial and extract it from all chans
all_cortex_trials = zeros(trial_lengths(2),length(cortex_chans)*size(band_data_2,1));
for trial = 1:trial_lengths(2)
    trial_temp = [];
    for chan = 1:length(cortex_chans)       
      trial_temp= [trial_temp ; nanmean(band_data_2(:,cortex_start:cortex_end,trial,cortex_chans(chan)),2)];
    end
    all_cortex_trials(trial,:) = trial_temp;
end

all_HC_trials = zeros(trial_lengths(2),size(band_data_2,1));
for trial = 1:trial_lengths(2)
    all_HC_trials(trial,:) = nanmean(band_data_2(:,HC_start:HC_end,trial,HC_chans(HC_CHAN)),2);
end



% remove same bad trials from both matrices
NC_bad = [];
for a = 1:size(all_cortex_trials,1)
    
    if sum(isnan(all_cortex_trials(a,:)))>0
        NC_bad = [NC_bad a];
    end
    
end

HC_bad = [];
for a = 1:size(all_HC_trials,1)
    
    if sum(isnan(all_HC_trials(a,:)))>0
        HC_bad = [HC_bad a];
    end
    
end

all_bad     = unique( [HC_bad NC_bad]);
trials_keep = ones (1, size(all_cortex_trials,1));
trials_keep(all_bad) = 0;

% clean matrices - remove nan trials
all_cortex_trials = all_cortex_trials(logical(trials_keep),:);
all_HC_trials     = all_HC_trials(logical(trials_keep),:);

% measure pair wise correl in cortex
rho_NC = corr(all_cortex_trials');

% measure pair wise coreel in HC
rho_HC = corr(all_HC_trials');

all_corr_NC=[];
all_corr_HC=[];
for a = 1:size(rho_NC,1)
    all_corr_NC = [all_corr_NC rho_NC(a,a+1:end)];
    all_corr_HC = [all_corr_HC rho_HC(a,a+1:end)];
end

% sort cortex from hi to low
[increaing_NC_corr I]= sort (all_corr_NC);


%
figure;
plot(increaing_NC_corr, all_corr_HC(I), '*r')
xlabel('EC pairwise correl', 'FontSize', 12,'FontWeight', 'bold' )
ylabel('Corresponding HC pairwase correl', 'FontSize', 12,'FontWeight', 'bold' )
A = [increaing_NC_corr ;all_corr_HC(I)]'
[R,P] = corrcoef(A);
corr_cals(HC_CHAN,:) = [R(2) P(2) ];
end