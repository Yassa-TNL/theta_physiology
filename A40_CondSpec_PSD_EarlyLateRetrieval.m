clear all;close all;clc
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
locks       = {'onset'  'response'}
reg         =  'OFC_FRO_TEMP_sig'%'TEMP' 'OFC_FRO_TEMP' OFC_FRO_TEMP_sig 'OFC_FRO_TEMP_CING_INS_EC_sig'
fn_ext      = 'cue_resp' %
phase       = 'retrieval'  %'encoding' 'retrieval'
DS          = 'yes';
ref         = 'LM';
clinical    = 'yes';
fs          = 500;
subj_list   = {'39' '44'  '57'  '63' '66' '84' '85' '87'};
cd('/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/LM_reref/tuning_correct_onset')
if strcmp('OFC_FRO_TEMP_CING_INS_EC_sig',reg)
load('subj_ALL_NC_specificty')
elseif strcmp('OFC_FRO_TEMP_sig',reg)
load('subj_NC_specificty.mat')
end

if strcmp('encoding',phase)
    exp_type    = phase;
else
    exp_type    = 'tuning_correct';
end

chan_psd_allsubj_onset_cond2    = cell(1,length(subj_list));
chan_psd_allsubj_response_cond2 = cell(1,length(subj_list));
chan_psd_allsubj_onset_cond3    = cell(1,length(subj_list));
chan_psd_allsubj_response_cond3 = cell(1,length(subj_list));
chan_psd_allsubj_onset_cond2_dB    = cell(1,length(subj_list));
chan_psd_allsubj_response_cond2_dB = cell(1,length(subj_list));
chan_psd_allsubj_onset_cond3_dB    = cell(1,length(subj_list));
chan_psd_allsubj_response_cond3_dB = cell(1,length(subj_list));

% loop thru onset and response
for iLock = 1:length(locks)
    lock = locks{iLock};
    
    % loop thru subj
    for sub_counter = 1:length(subj_list)
        subj =  subj_list{sub_counter};
        % indicate whether to pull from clin or research recordings
        if strcmp('84', subj) || strcmp('85', subj) || strcmp('87', subj)
            if strcmp('yes', clinical)
                fn_nm = '_clinical';
            else
                fn_nm = '_research';
            end
        else
            fn_nm = '';
        end
        cd(['/mnt/yassamri/iEEG/sandra/subj_' subj_list{sub_counter}])
        
        % load baseline
        if strcmp('39', subj) || strcmp('44', subj) || strcmp('57', subj) || strcmp('63', subj) || strcmp('66', subj) || strcmp('83', subj)
            load(['normalization_entire_recording_ref_' ref '_baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_' subj '_select_chan_3_fs_500.mat'])
        elseif strcmp('84', subj) || strcmp('85', subj) || strcmp('87', subj)
            load(['normalization_entire_recording_ref_' ref '_baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_' subj '_select_chan_3' fn_nm '_fs_500.mat'])
        end
        
        % load trial data
        if strcmp('yes',DS)
            load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3' fn_nm '_fs_500.mat'])
        elseif strcmp('',DS)
            load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3_NotDownsampled.mat'])
        end
        
        % get chans
        [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
            ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj);
        
        if strcmp('OFC',reg)
            chan_idx = OFC_chan_idx;
        elseif strcmp('FRO',reg)
            chan_idx = fro_chan_idx;
        elseif strcmp('TEMP',reg)
            chan_idx = temp_chan_idx;
        elseif strcmp('CING',reg)
            chan_idx = cingulate_chan_idx;
        elseif strcmp('INS',reg)
            chan_idx = insula_chan_idx;
        elseif strcmp('EC',reg)
            chan_idx = EC_chan_idx;
        elseif strcmp('HC',reg)
            chan_idx = HC_chan_idx;
        elseif strcmp('CA1',reg)
            chan_idx = CA1_chan_idx;
        elseif strcmp('CA3',reg)
            chan_idx = CA3_chan_idx;
        elseif strcmp('NC',reg)
            chan_idx = NC_chan_idx;
        elseif strcmp('OFC_FRO_TEMP', reg) || strcmp('OFC_FRO_TEMP_sig', reg)
             chan_idx = [OFC_chan_idx fro_chan_idx temp_chan_idx];
        elseif strcmp('OFC_FRO_TEMP_CING_INS_EC_sig', reg)
             chan_idx = [OFC_chan_idx fro_chan_idx temp_chan_idx cingulate_chan_idx insula_chan_idx EC_chan_idx];
        end
        if strcmp('OFC_FRO_TEMP_CING_INS_EC_sig',reg)
            sig_chan_log_vec = logical(subj_ALL_NC_specificty{sub_counter}(:,2)<0.05);
            chan_idx = chan_idx(ismember(chan_idx, subj_ALL_NC_specificty{sub_counter}(sig_chan_log_vec,1)));
        elseif strcmp('OFC_FRO_TEMP_sig',reg)
            sig_chan_log_vec = logical(subj_NC_specificty{sub_counter}(:,2)<0.05);
            chan_idx = chan_idx(ismember(chan_idx, subj_NC_specificty{sub_counter}(sig_chan_log_vec,1)));
        end
        if isempty(chan_idx)
            continue
        end
        chan_artifact_thresh = chan_artifact_thresh(chan_idx',:);
        if strcmp('onset',lock)
            dur=1;
            strt_idx  = round(pre_stim*fs);
            end_idx   = round((pre_stim+dur)*fs);
        elseif strcmp('response', lock)
            strt_idx = 0.5*fs;
            end_idx  = 1.5*fs;
        end
        
        % loop thru conds
        for cond_num = 2:3
            % get cond data
            [cond1,cond2,cond3,cond4,cond5,cond6] = GetCondData(subj, exp_type, lock, DS, fn_nm, ref);
            
            % index data with an epoch and desired chans
            
            if  cond_num==1
                cond  = 'Repeat';
                cond1 = cond1(:,strt_idx:end_idx,chan_idx);
                trial_data = cond1;
                disp('cond1-Repeat')
                
            elseif cond_num==2
                cond  = 'PatternCompletion';
                cond2 = cond2(:,strt_idx:end_idx,chan_idx);
                trial_data = cond2;
                disp('cond2-pattern comp')
                
            elseif cond_num==3
                cond  = 'PatternSeperation';
                cond3 = cond3(:,strt_idx:end_idx,chan_idx);
                trial_data = cond3;
                disp('cond3-pattern sep')
                
            elseif cond_num==4
                cond  = 'New';
                cond4 = cond4(:,strt_idx:end_idx,chan_idx);
                trial_data = cond4;
                disp('cond4-New')
            end
            

            % remove trials with artifacts
            bad_trials = [];
            for elec = 1:size(trial_data, 3)
                for trials = 1:size(trial_data,1)
                    trial_temp = trial_data(trials, :,elec);
                    
                    % if trail has artifact = nan
                    idx_above_thresh = (trial_temp)>chan_artifact_thresh(elec,2);
                    idx_below_thresh = (trial_temp)<chan_artifact_thresh(elec,1);
                    total_idx = idx_above_thresh+idx_below_thresh;
                    
                    if sum(total_idx)>0 % if trial is artifactual
                        trial_data(trials,:,elec) = nan;
                    end
                end
            end
            
            %% add theta artificially to data to test pqwlch
            
%             dt = 1/fs;                   % seconds per sample
%             t = (0:dt:dur)';     % seconds
%             %%Sine wave:
%             Fc = 5;                     % hertz
%             x = cos(2*pi*Fc*t);
% 
%             for elec = 1:size(trial_data, 3)
%                 for trials = 1:size(trial_data,1)
%                   trial_data(trials, :,elec) = trial_data(trials, :,elec)+ x';
%                     
%                 end
%             end
%             
            % get indiv patient psd, average across trials
            chan_psd = nan(length(freq),size(trial_data,3));
            for iChan = 1:size(trial_data,3)
                [Pxx,F] = pwelch(trial_data(:, :,iChan)',length(trial_data(1, :,iChan)),0,flip(freq'),fs);% input data = timeXtrials                             %  input (data, win size,overlap, NFFT, fs)
                chan_psd (:,iChan) = nanmean(Pxx,2);
            end
            
             % get indiv patient psd, pool trials
            trial_chan_psd_raw = cell(1,size(trial_data,3));
            trial_chan_psd_db  = cell(1,size(trial_data,3));

            for iChan = 1:size(trial_data,3)
                [Pxx,F] = pwelch(trial_data(:, :,iChan)',length(trial_data(1, :,iChan)),0,flip(freq'),fs);% input data = timeXtrials                             %  input (data, win size,overlap, NFFT, fs)
                trial_chan_psd_raw {iChan} = Pxx;
                trial_chan_psd_db  {iChan} = 10*log10(Pxx);
            end
            
            if iLock ==1 && cond_num==2
                chan_psd_allsubj_onset_cond2{sub_counter}    = chan_psd;
                chan_psd_allsubj_onset_cond2_dB{sub_counter} = 10*log10(chan_psd);
                
                trial_chan_psd_allsubj_raw_cond2{sub_counter} = cat(2,trial_chan_psd_raw{:});
                trial_chan_psd_allsubj_db_cond2 {sub_counter} = cat(2,trial_chan_psd_db{:});

            elseif iLock ==2 && cond_num==2
                chan_psd_allsubj_response_cond2{sub_counter}    = chan_psd;
                chan_psd_allsubj_response_cond2_dB{sub_counter} =  10*log10(chan_psd);
                
            elseif iLock ==1 && cond_num==3
                chan_psd_allsubj_onset_cond3{sub_counter}       = chan_psd;
                chan_psd_allsubj_onset_cond3_dB{sub_counter} =  10*log10(chan_psd);
                trial_chan_psd_allsubj_raw_cond3{sub_counter} = cat(2,trial_chan_psd_raw{:});
                trial_chan_psd_allsubj_db_cond3 {sub_counter} = cat(2,trial_chan_psd_db{:});
                
            elseif iLock ==2 && cond_num==3
                chan_psd_allsubj_response_cond3{sub_counter}    = chan_psd;
                chan_psd_allsubj_response_cond3_dB{sub_counter} =  10*log10(chan_psd);
            end
            
            clear trial_data
        end
        
        
    end
end

% power
chan_psd_pooled_allsubj_onset2 = cat(2,chan_psd_allsubj_onset_cond2{:});
chan_psd_pooled_allsubj_resp2  = cat(2,chan_psd_allsubj_response_cond2{:});
chan_psd_pooled_allsubj_onset3 = cat(2,chan_psd_allsubj_onset_cond3{:});
chan_psd_pooled_allsubj_resp3  = cat(2,chan_psd_allsubj_response_cond3{:});

% decible
chan_psd_pooled_allsubj_onset2_dB = cat(2,chan_psd_allsubj_onset_cond2_dB{:});
chan_psd_pooled_allsubj_resp2_dB  = cat(2,chan_psd_allsubj_response_cond2_dB{:});
chan_psd_pooled_allsubj_onset3_dB = cat(2,chan_psd_allsubj_onset_cond3_dB{:});
chan_psd_pooled_allsubj_resp3_dB  = cat(2,chan_psd_allsubj_response_cond3_dB{:});

% power diff - used for paper
diff_mtx_cond3 = chan_psd_pooled_allsubj_onset3 - chan_psd_pooled_allsubj_resp3;
diff_mtx_cond2 = chan_psd_pooled_allsubj_onset2 - chan_psd_pooled_allsubj_resp2;

% decible diff
diff_mtx_cond3_dB = chan_psd_pooled_allsubj_onset3_dB - chan_psd_pooled_allsubj_resp3_dB;
diff_mtx_cond2_dB = chan_psd_pooled_allsubj_onset2_dB - chan_psd_pooled_allsubj_resp2_dB;

trial_chan_psd_allsubj_raw_cond2_pooled = cat(2,trial_chan_psd_allsubj_raw_cond2{:});
trial_chan_psd_allsubj_raw_cond3_pooled = cat(2,trial_chan_psd_allsubj_raw_cond3{:});
trial_chan_psd_allsubj_db_cond2_pooled = cat(2,trial_chan_psd_allsubj_db_cond2{:});
trial_chan_psd_allsubj_db_cond3_pooled = cat(2,trial_chan_psd_allsubj_db_cond3{:});

% *** difference bewteen conditions with chan define for desired figures
% cond2=chan_psd_pooled_allsubj_onset2';
% cond3=chan_psd_pooled_allsubj_onset3';


% cond2=trial_chan_psd_allsubj_db_cond2_pooled';
% cond3=trial_chan_psd_allsubj_db_cond3_pooled';

% used for paper
cond2=diff_mtx_cond2';
cond3=diff_mtx_cond3';


desired_freq_lo = 3;
desired_freq_hi = 200;
freqLogVec=(F>desired_freq_lo)&(F<desired_freq_hi);
F_new = F(freqLogVec);

figure;
hold on;
stdshade(cond3(:,freqLogVec),.1,'m',F(freqLogVec,:),[] ,[], []);
h1 = plot(F(freqLogVec,:),nanmean(cond3(:,freqLogVec),1), 'm', 'LineWidth', 2);
stdshade(cond2(:,freqLogVec),.1,'b',F(freqLogVec,:),[] ,[], []);
h2 = plot(F(freqLogVec,:),nanmean(cond2(:,freqLogVec),1), 'b', 'LineWidth', 2);
xticks(round(F_new(1:10:end)))
%title([reg ' subj' subj],'Interpreter', 'none')
xlabel('freq');
%ylabel('onset psd - response psd')
xlim([F_new(1) F_new(end)])
xlim([2.9 30])

set(gca, 'FontSize', 16, 'FontWeight', 'bold')
print('-clipboard','-dbitmap')
%% 
figure;
for i = 1:34

hold on;
stdshade(cond3(i,:),.1,'m',F,[] ,[], []);
h1 = plot(F,nanmean(cond3(i,:),1), 'm', 'LineWidth', 2);
stdshade(cond2(i,:),.1,'b',F,[] ,[], []);
h2 = plot(F,nanmean(cond2(i,:),1), 'b', 'LineWidth', 2);
xlabel('freq');
%ylabel('onset psd - response psd')
xlim([2.9 30])
pause
close all
end
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
%%
%[zmap,zmapthresh,zmapthresh_for_plot] = permutation_testing_vector(cond3(freqLogVec,:)',cond2(freqLogVec,:)', 1000);
ylim([-20 50])
y=ylim
if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
end
scale=y(2)
%plot(F(freqLogVec,:),scale*zmapthresh_for_plot, 'LineWidth', 2, 'Color', 'k')

%legend([h1 h2],{'lure+','lure-'})
print('-clipboard','-dbitmap')
cd('/mnt/yassamri/iEEG/sandra/GroupFigures')

%%
cond_num=2

if cond_num==2
    cond2=chan_psd_pooled_allsubj_onset2_dB;
cond3=chan_psd_pooled_allsubj_resp2_dB;
    title_val = [reg ' lure-']
else
    cond2=chan_psd_pooled_allsubj_onset3_dB;
cond3=chan_psd_pooled_allsubj_resp3_dB;
    title_val = [reg ' lure+']
end



desired_freq_lo = 3;
desired_freq_hi = 20;
freqLogVec=(F>desired_freq_lo)&(F<desired_freq_hi);
F_new = F(freqLogVec);
figure;
hold on;
stdshade(cond2(freqLogVec,:)',.1,'m',F(freqLogVec,:),[] ,[], []);
h1 = plot(F(freqLogVec,:),nanmean(cond2(freqLogVec,:)',1), 'r', 'LineWidth', 2);
stdshade(cond3(freqLogVec,:)',.1,'b',F(freqLogVec,:),[] ,[], []);
h2 = plot(F(freqLogVec,:),nanmean(cond3(freqLogVec,:)',1), 'g', 'LineWidth', 2);
xticks(round(F_new(1:10:end)))
title(title_val)
xlabel('freq');ylabel('onset psd - response psd')
xlim([F_new(1) F_new(end)])
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

%stats
[zmap,zmapthresh,zmapthresh_for_plot] = permutation_testing_vector(cond3(freqLogVec,:)',cond3(freqLogVec,:)', 1000);
y=ylim
scale=y(2)
plot(F(freqLogVec,:),scale*zmapthresh_for_plot, 'LineWidth', 2, 'Color', 'k')

legend([h1 h2],{'onset','response'})
print('-clipboard','-dbitmap')







