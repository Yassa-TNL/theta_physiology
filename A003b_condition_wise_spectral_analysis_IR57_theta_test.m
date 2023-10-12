%% indicate variables
clear all ;close all; clc
theta_test = 'yes'
subj = '57'
freq_analysis = 'wavelet'
select_chan   = 3
lock ='onset' %'response' % 
exp_type ='tuning_correct' % {'study_test' 'tuning' 'tuning_correct' 'tuning_incorrect' 'indoor_outdoor'}
minfreq = 3;
maxfreq = 200;
cd (['/mnt/yassamri/iEEG/sandra/subj_' subj])
if strcmp('39', subj)
    load('baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_39_select_chan_3.mat')
    %load('trial_data_subj_onset_39_ref__select_chan_3.mat')
elseif strcmp('44', subj)
    load('normalization_entire_recording_ref__baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_44_select_chan_3.mat')
   % load('trial_data_subj_onset_44_ref__select_chan_3.mat')

elseif strcmp('84', subj)
    load('baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_84.mat')
    % load('trial_data_subj84.mat')
    load('timestamps.mat')
elseif strcmp('83', subj)
    load('normalization_entire_recording_ref__baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_83_select_chan_3.mat')
    load('trial_data_subj_onset_83_ref__select_chan_3.mat')
elseif strcmp('63', subj)
    load('normalization_entire_recording_ref__baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_63_select_chan_3.mat')
   % load('trial_data_subj_onset_63_ref__select_chan_3.mat')
elseif strcmp('66', subj)
    load('normalization_entire_recording_baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_66_select_chan_3.mat')
   % load('trial_data_subj_onset_66_ref__select_chan_3.mat')

elseif strcmp('57', subj)
    load('baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_57_select_chan_3')
    %load('trial_data_subj_onset_57_ref__select_chan_3.mat')

end

if strcmp('response', lock)
    load(['trial_data_subj_response_' subj '_ref__select_chan_3.mat']);
    % remove trils w/ no resp
    bd_trls = zeros(1,length(responses));
    bd_trls(responses<.2) = 1;
    trial_data(logical(bd_trls),:,:) = nan;
elseif strcmp('onset', lock)
    load(['trial_data_subj_onset_' subj '_ref__select_chan_3.mat'])
    if strcmp('yes', theta_test)
        load(['trial_data_subj_onset_' subj '_ref__select_chan_3_theta_onset_test.mat'])
    end
end
load(['behavior_subj' subj '.mat'])



tic


% organize data into condition
study_trial_data = trial_data(1:nTrials_study,:,:);
test_trial_data  = trial_data(nTrials_study+1:end,:,:);

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
lures              =[1 2 4 5];
diff               =[1 2 ];
lure_diff          =logical(ismember(testing_behav_matrix(:,2),diff));
lure_all           =logical(ismember(testing_behav_matrix(:,2),lures));

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
    trial_lengths = [size(cond1,1) size(cond2,1) size(cond3,1) size(cond4,1) size(cond5,1) size(cond6,1)];

elseif strcmp('tuning_correct',exp_type)
    
    cond1= test_trial_data(repeat_log_vec&corr_resp_log_vec ==1,:,:);
    cond2= test_trial_data(lure_all==1 & corr_resp_log_vec ==0,:,:);  %pattern comp
    cond3= test_trial_data(lure_all==1 & corr_resp_log_vec ==1,:,:); %lure 
    cond4= test_trial_data(new_log_vec & corr_resp_log_vec ==1,:,:);
    trial_lengths = [size(cond1,1) size(cond2,1) size(cond3,1) size(cond4,1) ];

%     cond1= test_trial_data(repeat_log_vec&corr_resp_log_vec ==1,:,:);
%     cond2= test_trial_data(lure_1_log_vec&corr_resp_log_vec ==1,:,:);
%     cond3= test_trial_data(lure_2_log_vec&corr_resp_log_vec ==1,:,:);
%     cond4= test_trial_data(lure_4_log_vec&corr_resp_log_vec ==1,:,:);
%     cond5= test_trial_data(lure_5_log_vec&corr_resp_log_vec ==1,:,:);
%     cond6= test_trial_data(new_log_vec,:,:);
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
norm_freq_acrs_chan_cond_3 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond3,1), chan_counter);
norm_freq_acrs_chan_cond_4 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond4,1), chan_counter);

trial_log_1_artifact = cell(chan_counter,1);
trial_log_2_artifact = cell(chan_counter,1);
trial_log_3_artifact = cell(chan_counter,1);
trial_log_4_artifact = cell(chan_counter,1);

if strcmp('tuning',exp_type)
    trial_log_5_artifact = cell(chan_counter,1);
    trial_log_6_artifact = cell(chan_counter,1);
    norm_freq_acrs_chan_cond_5 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond5,1), chan_counter);
    norm_freq_acrs_chan_cond_6 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond6,1), chan_counter);
end

%%
cntr=0;
for chan =1:MTL_chan_idx %chan_counter    
    cntr = cntr+1;
    % first condition
parfor trial = 1:trial_lengths(1)
   
    trial_temp = cond1(trial, :,chan);
    
    % if trail has artifact = nan
    idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
    idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
    total_idx = idx_above_thresh+idx_below_thresh;
    
    % indicate a nan trial as an artifact    
    if sum(isnan(trial_temp))==length(trial_temp)
        total_idx = ones(1,length(trial_temp));
    end
    
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
        
    norm_freq_acrs_chan_cond_1(:, :, trial, cntr) =  norm_freq;
end
trial_log_1_artifact{cntr}=  trial_log_1_artifact_temp';


    % second condition
parfor trial = 1:trial_lengths(2)
   
    trial_temp = cond2(trial, :,chan);
    
    % if trail has artifact = nan
    idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
    idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
    total_idx        = idx_above_thresh+idx_below_thresh;
    
    % indicate a nan trial as an artifact    
    if sum(isnan(trial_temp))==length(trial_temp)
        total_idx = ones(1,length(trial_temp));
    end
    
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
   
    norm_freq_acrs_chan_cond_2(:, :, trial, cntr) =  norm_freq;
end
trial_log_2_artifact{cntr}=  trial_log_2_artifact_temp';

    % third condition
parfor trial = 1:trial_lengths(3)
   trial_temp = cond3(trial, :,chan);
   
    % if trail has artifact = nan
    idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
    idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
    total_idx = idx_above_thresh+idx_below_thresh;
    
    % indicate a nan trial as an artifact    
    if sum(isnan(trial_temp))==length(trial_temp)
        total_idx = ones(1,length(trial_temp));
    end
    
    trial_log_3_artifact_temp(trial)= sum(total_idx); 
    
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
    
    norm_freq_acrs_chan_cond_3(:, :, trial, cntr) =  norm_freq;
end
trial_log_3_artifact{cntr}=  trial_log_3_artifact_temp';



    % fourth condition
parfor trial = 1:trial_lengths(4)
    trial_temp = cond4(trial, :,chan);
    
    % if trail has artifact = nan
    idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
    idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
    total_idx = idx_above_thresh+idx_below_thresh;
    
    % indicate a nan trial as an artifact    
    if sum(isnan(trial_temp))==length(trial_temp)
        total_idx = ones(1,length(trial_temp));
    end
    
    trial_log_4_artifact_temp(trial)= sum(total_idx); 
    
    
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
    
    norm_freq_acrs_chan_cond_4(:, :, trial, cntr) =  norm_freq;
end

trial_log_4_artifact{cntr}=  trial_log_4_artifact_temp';
if strcmp('tuning',exp_type)
    parfor trial = 1:trial_lengths(5)
        trial_temp = cond5(trial, :,chan);
        
        % if trail has artifact = nan
        idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
        idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
        total_idx = idx_above_thresh+idx_below_thresh;
        
        if sum(isnan(trial_temp))==length(trial_temp)
            total_idx = ones(1,length(trial_temp));
        end
        
        trial_log_5_artifact_temp(trial)= sum(total_idx);
        
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
        
        norm_freq_acrs_chan_cond_5(:, :, trial, chan) =  norm_freq;
    end
    trial_log_5_artifact{chan}=  trial_log_5_artifact_temp';
    
    parfor trial = 1:trial_lengths(6)
        trial_temp = cond6(trial, :,chan);
        
        % if trail has artifact = nan
        idx_above_thresh = (trial_temp)>chan_artifact_thresh(chan,2);
        idx_below_thresh = (trial_temp)<chan_artifact_thresh(chan,1);
        total_idx = idx_above_thresh+idx_below_thresh;
        
        if sum(isnan(trial_temp))==length(trial_temp)
            total_idx = ones(1,length(trial_temp));
        end
        
        trial_log_6_artifact_temp(trial)= sum(total_idx);
        
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
        
        norm_freq_acrs_chan_cond_6(:, :, trial, chan) =  norm_freq;
    end
    trial_log_6_artifact{chan}=  trial_log_6_artifact_temp';
end
end



% num of clean trials
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


% average across trials within chans
edge_points = 200;

% initialize matrix for all frequencies
mn_acrs_trials_cond1 = nan(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2)-2*edge_points, chan_counter);
mn_acrs_trials_cond2 = nan(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2)-2*edge_points, chan_counter);
mn_acrs_trials_cond3 = nan(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2)-2*edge_points, chan_counter);
mn_acrs_trials_cond4 = nan(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2)-2*edge_points, chan_counter);

% average across trials within chans
for elec = 1:size(norm_freq_acrs_chan_cond_1,4) %loop thru chans
    mn_acrs_trials_cond1 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_1(:, edge_points+1:end-edge_points, :, elec), 3);
    mn_acrs_trials_cond2 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_2(:, edge_points+1:end-edge_points, :, elec), 3);
    mn_acrs_trials_cond3 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_3(:, edge_points+1:end-edge_points, :, elec), 3);
    mn_acrs_trials_cond4 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_4(:, edge_points+1:end-edge_points, :, elec), 3);
end

% mean acrs chans
indiv_freq_chans_cond1 = nanmean(mn_acrs_trials_cond1,3);
indiv_freq_chans_cond2 = nanmean(mn_acrs_trials_cond2,3);
indiv_freq_chans_cond3 = nanmean(mn_acrs_trials_cond3,3);
indiv_freq_chans_cond4 = nanmean(mn_acrs_trials_cond4,3);


% remove edge effects
indiv_freq_chans_cond1 = indiv_freq_chans_cond1(:, edge_points+1:end-edge_points);
indiv_freq_chans_cond2 = indiv_freq_chans_cond2(:, edge_points+1:end-edge_points);
indiv_freq_chans_cond3 = indiv_freq_chans_cond3(:, edge_points+1:end-edge_points);
indiv_freq_chans_cond4 = indiv_freq_chans_cond4(:, edge_points+1:end-edge_points);

% mn across right and left HC
if select_chan== 1 || 2
    
    right_HC_cond1 = nanmean(mn_acrs_trials_cond1(:,:, logical(right_chans)), 3);
    right_HC_cond2 = nanmean(mn_acrs_trials_cond2(:,:, logical(right_chans)), 3);
    right_HC_cond3 = nanmean(mn_acrs_trials_cond3(:,:, logical(right_chans)), 3);
    right_HC_cond4 = nanmean(mn_acrs_trials_cond4(:,:, logical(right_chans)), 3);

    
    left_HC_cond1 = nanmean(mn_acrs_trials_cond1(:,:, logical(left_chans)), 3);
    left_HC_cond2 = nanmean(mn_acrs_trials_cond2(:,:, logical(left_chans)), 3);
    left_HC_cond3 = nanmean(mn_acrs_trials_cond3(:,:, logical(left_chans)), 3);
    left_HC_cond4 = nanmean(mn_acrs_trials_cond4(:,:, logical(left_chans)), 3);

end


if strcmp('tuning', exp_type)

    for a = 1:size(trial_log_5_artifact,1)
        chan_trials_log_5 (a)= sum(trial_log_5_artifact{a}==0);
        
    end
    total_tiral_log_5 = sum(chan_trials_log_5);
    
    for a = 1:size(trial_log_6_artifact,1)
        chan_trials_log_6 (a)= sum(trial_log_6_artifact{a}==0);
        
    end
    total_tiral_log_6 = sum(chan_trials_log_6);
    
    mn_acrs_trials_cond5 = nan(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2)-2*edge_points, chan_counter);
    mn_acrs_trials_cond6 = nan(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2)-2*edge_points, chan_counter);
    
    for elec = 1:chan_counter
        mn_acrs_trials_cond5 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_5(:, edge_points+1:end-edge_points, :, elec), 3);
        mn_acrs_trials_cond6 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_6(:, edge_points+1:end-edge_points, :, elec), 3);
    end
    
    indiv_freq_chans_cond5 = nanmean(mn_acrs_trials_cond5,3);
    indiv_freq_chans_cond6 = nanmean(mn_acrs_trials_cond6,3);
    
    indiv_freq_chans_cond5 = indiv_freq_chans_cond5(:, edge_points+1:end-edge_points);
    indiv_freq_chans_cond6 = indiv_freq_chans_cond6(:, edge_points+1:end-edge_points);
    
    if select_chan== 1 || 2
        right_HC_cond5 = nanmean(mn_acrs_trials_cond5(:,:, logical(right_chans)), 3);
        right_HC_cond6 = nanmean(mn_acrs_trials_cond6(:,:, logical(right_chans)), 3);
        
        left_HC_cond5 = nanmean(mn_acrs_trials_cond5(:,:, logical(left_chans)), 3);
        left_HC_cond6 = nanmean(mn_acrs_trials_cond6(:,:, logical(left_chans)), 3);
    end
    
end

% figure params

mn = -.4
mx = .4
tickmarks = 1:20:length(freq);

% get start and end trial indices
if strcmp ('onset',lock)
    if strcmp('yes',theta_test)
    pre_stim  = 2;
    else
    pre_stim  = 0.5;
    end
    post_stim = 2; % changed to 1 sec. ISI is 1 sec?
    strt_time = on_idx' - pre_stim*fs;
    end_time  = on_idx' + post_stim*fs;
     
elseif strcmp ('offset',lock)
    pre_stim = 1;
    post_stim = .5;
    strt_time = off_idx' - pre_stim*fs
    end_time  = off_idx' + post_stim*fs;
end

cond_num = length(trial_lengths);


toc
%% plot across all chans
path_name='/mnt/yassamri/iEEG/sandra'
figures_dir = [path_name '/subj_' num2str(subj) '/figures/' norm{normalization} '/' exp_type '_' lock]
if ~exist(figures_dir, 'dir')
mkdir ([path_name '/subj_' num2str(subj) '/figures/' norm{normalization} '/' exp_type '_' lock])
end
cd([path_name '/subj_' num2str(subj) '/figures/' norm{normalization} '/' exp_type '_' lock])

%%
% plot first condition
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), indiv_freq_chans_cond1)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), indiv_freq_chans_cond2)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('++++ sim')
%title('Pattern Comp: incorr lure')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), indiv_freq_chans_cond3)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('+++ sim')
%title('Pattern Sep: corr lure')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), indiv_freq_chans_cond4)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('++ sim')
%title('new')

if strcmp('tuning', exp_type)
subplot (cond_num,1,5)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), indiv_freq_chans_cond5)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('+ sim')
%title('Pattern Sep: corr lure')

subplot (cond_num,1,6)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), indiv_freq_chans_cond6)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')

end
saveas(gcf, 'all_chans.jpg') 



%% plot right chans

figure
% plot first condition
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), right_HC_cond1)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), right_HC_cond2)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('Pattern Comp: incorr lure')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), right_HC_cond3)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('Pattern Sep: corr lure')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), right_HC_cond4)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')

saveas(gcf, 'right_HC_chans.jpg') 

%% Left chans
figure
% plot first condition
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), left_HC_cond1)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), left_HC_cond2)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('pattern Comp: incorr lure')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), left_HC_cond3)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('Pattern Sep: corr lure')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), left_HC_cond4)
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
mx = .8
mn = -1
for chan =1:5

ax = figure 
%('visible', 'off')
suptitle([chan_label(chan) 'mean RT: ' num2str(mean(responses(nTrials_study+1:end)))...
    ' min RT: ' num2str(min(responses(nTrials_study+1:end)))])
%suptitle('Inferior Frontal Gyrus')
%suptitle('CA3 / DG')

subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond1(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['repeat: ' num2str(sum(trial_log_1_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond2(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['Pattern Comp - incorr lure: ' num2str(sum(trial_log_2_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond3(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['Pattern Sep - corr lure: ' num2str(sum(trial_log_3_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond4(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['new: ' num2str(sum(trial_log_4_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

set(gcf, 'PaperUnits', 'inches');
x_width=7 ;y_width=12
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf, [chan_label{chan} '_chan' num2str(chan) '.bmp']) 
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
close all
end




%% larger plots

loop_string = {'LTI9' 'LTI10'}
chan =  find(~cellfun(@isempty,regexp(chan_label, loop_string{elec})))
suptitle(['LIFG ' chan_label{chan}(4:end-3) ])
chan = chan(1);


mx   = 1
mn   = -1
for chan = 1
f=figure('visible', 'off')
suptitle(chan_label{chan}(4:end-3))
subplot(1,4,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond1(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['repeat: ' num2str(sum(trial_log_1_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

subplot(1,4,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond2(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['Pattern Comp - incorr lure: ' num2str(sum(trial_log_2_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

subplot(1,4,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond3(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['Pattern Sep - corr lure: ' num2str(sum(trial_log_3_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

subplot(1,4,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond4(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['new: ' num2str(sum(trial_log_4_artifact{chan}==0)) ' trials'])
colormap jet
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
x_width=22 ;y_width=8

set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

saveas(gcf, ['chan' num2str(chan) '.bmp']) 
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
close all

end


%% subtraction plots


mx = .81
figure
suptitle(chan_label{chan}(4:end-3))

subplot(3,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq),(mn_acrs_trials_cond3(:,:,chan)-mn_acrs_trials_cond1(:,:,chan)))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['lure corr - repeat '])
colormap jet
set(gca, 'FontSize',9)

subplot(3,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq),mn_acrs_trials_cond3(:,:,chan)-mn_acrs_trials_cond2(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['lure corr - lure incorr '])
colormap jet
set(gca, 'FontSize',9)

subplot(3,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq),mn_acrs_trials_cond3(:,:,chan)-mn_acrs_trials_cond4(:,:,chan))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(['lure corr - new ' ])
colormap jet
set(gca, 'FontSize',9)
%% plot trials
cond_temp = []
cond_temp = cond3;
for chan =16:19
    figure
    fig_name = chan_label(chan)
    title(fig_name{1}(4:8))
    hold on
    yaxis_shift = 300
    
     plot(-.5:1/fs:2, nanmean(cond1(:,:,chan),1),'b', 'LineWidth', 2)
     plot(-.5:1/fs:2, nanmean(cond2(:,:,chan),1),'g', 'LineWidth', 2)
     plot(-.5:1/fs:2, nanmean(cond3(:,:,chan),1),'r', 'LineWidth', 2)
     plot(-.5:1/fs:2, nanmean(cond4(:,:,chan),1),'m', 'LineWidth', 2)
legend({'repeat', 'incorr lure', 'corr lure', 'new'})

    
%     for trial = 1:size(cond_temp,1)
%         
%         plot(-.5:1/fs:2, cond_temp(trial,:,chan)+(yaxis_shift*trial))
%        % plot(-.5:1/fs:1, cond_temp(trial,:,chan))
% 
%     end
%     saveas(gcf, [cell2mat(chan_label(chan)) 'time.jpg'])
%     close all
end
