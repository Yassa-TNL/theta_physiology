clear all; close all; clc

fn_ext      = 'cue_resp' %
phase       = 'retrieval'  %'encoding' 'retrieval'
lock        = 'onset'   %'onset'    'response'
reg1_name   = 'TEMP' % OFC, FRO, TEMP
reg2_name   = 'HC';
DS          = 'yes';
ref         = 'LM';
clinical    = 'yes';
lag         = '';
subj_list = {'39' '44'  '57'  '63' '66' '84' '85' '87'};
addpath('/tmp/yassamri/iEEG/sandra/analysis_pipeline_final')

if strcmp('encoding',phase)
    cond_nums = 1:2;
    time1_list = 0:0.01:1.5;
    time2_list = 0.5:0.01:2;
elseif strcmp('retrieval',phase)&&strcmp('onset',lock)
    cond_nums = 2:3;
    time1_list = 0
    time2_list = 1
elseif strcmp('retrieval',phase)&&strcmp('response',lock)
    cond_nums = 2:3;
    time1_list = 0.5;
    time2_list = 1.5;
end


for iTime = 1:length(time1_list)
    
    for cond_num    = cond_nums
        time1       = time1_list(iTime); % response is 0.5 to 1.5 to get 1 sec be4 resp
        time2       = time2_list(iTime);
        if strcmp('encoding',phase)
            exp_type    = phase; 
        else
            exp_type    = 'tuning_correct'; 
        end

        
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
            cd(['/tmp/yassamri/iEEG/sandra/subj_' subj_list{sub_counter}])
            
            % load baseline
            if strcmp('39', subj) || strcmp('44', subj) || strcmp('57', subj) || strcmp('63', subj) || strcmp('66', subj) || strcmp('83', subj)
                load(['normalization_entire_recording_ref_' ref '_baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_' subj '_select_chan_3_fs_500.mat'])
            elseif strcmp('84', subj) || strcmp('85', subj) || strcmp('87', subj)
                load(['normalization_entire_recording_ref_' ref '_baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_' subj '_select_chan_3' fn_nm '_fs_500.mat'])
            end
            
            if strcmp('yes',DS)
                load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3' fn_nm '_fs_500.mat'])
            elseif strcmp('',DS)
                load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3_NotDownsampled.mat'])
            end
            
            % get cond data
            [cond1,cond2,cond3,cond4,cond5,cond6] = GetCondData(subj, exp_type, lock, DS, fn_nm, ref);
            
            % get chans
            [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
                ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj);
            
            if strcmp('OFC',reg1_name)
                reg1 = OFC_chan_idx;
            elseif strcmp('FRO',reg1_name)
                reg1 = fro_chan_idx;
            elseif strcmp('TEMP',reg1_name)
                reg1 = temp_chan_idx;
            elseif strcmp('CING',reg1_name)
                reg1 = cingulate_chan_idx;
            elseif strcmp('INS',reg1_name)
                reg1 = insula_chan_idx;
            elseif strcmp('EC',reg1_name)
                reg1 = EC_chan_idx;
            elseif strcmp('HC',reg1_name)
                reg1 = HC_chan_idx;
            elseif strcmp('CA1',reg1_name)
                reg1 = CA1_chan_idx;
            end
            
            if strcmp('HC',reg2_name)
                reg2 = HC_chan_idx;
            elseif strcmp('CA3',reg2_name)
                reg2 = CA3_chan_idx;
            elseif strcmp('CA1',reg2_name)
                reg2 = CA1_chan_idx;
            elseif strcmp('OFC',reg2_name)
                reg2 = OFC_chan_idx;
            elseif strcmp('FRO',reg2_name)
                reg2 = fro_chan_idx;
            elseif strcmp('TEMP',reg2_name)
                reg2 = temp_chan_idx;
            elseif strcmp('CING',reg2_name)
                reg2 = cingulate_chan_idx;
            elseif strcmp('INS',reg2_name)
                reg2 = insula_chan_idx;
            elseif strcmp('EC',reg2_name)
                reg2 = EC_chan_idx;
            end

            
            allreg = [reg1 reg2];
            length_reg1 = length(reg1);length_reg2 = length(reg2);
            if isempty(reg1) || isempty(reg2)
                continue
            end
            % index data with an epoch and desired chans
            % choose lock
            if strcmp('tuning_correct',exp_type)
                if strcmp('onset',lock)
                    strt_idx  = round((pre_stim+time1)*fs);
                    end_idx   = round((pre_stim+time2)*fs);
                elseif strcmp('response', lock)
                    strt_idx = time1*fs;
                    end_idx  = time2*fs;
                end
                if  cond_num==1
                    cond  = 'Repeat';
                    cond1 = cond1(:,strt_idx:end_idx,allreg);
                    trial_data = cond1;
                    disp('cond1-Repeat')
                    
                elseif cond_num==2
                    cond  = 'PatternCompletion';
                    cond2 = cond2(:,strt_idx:end_idx,allreg);
                    trial_data = cond2;
                    disp('cond2-pattern comp')
                    
                elseif cond_num==3
                    cond  = 'PatternSeperation';
                    cond3 = cond3(:,strt_idx:end_idx,allreg);
                    trial_data = cond3;
                    disp('cond3-pattern sep')
                    
                elseif cond_num==4
                    cond  = 'New';
                    cond4 = cond4(:,strt_idx:end_idx,allreg);
                    trial_data = cond4;
                    disp('cond4-New')
                end
                
            elseif  strcmp('encoding',exp_type)
                strt_idx  = round((pre_stim+time1)*fs);
                end_idx   = round((pre_stim+time2)*fs);
                if cond_num==1
                    cond  = 'PatternSeperation';
                    cond1 = cond1(:,strt_idx:end_idx,allreg);
                    trial_data = cond1;
                    disp('cond1-pattern sep')
                elseif cond_num==2
                    cond  = 'PatternCompletion';
                    cond2 = cond2(:,strt_idx:end_idx,allreg);
                    trial_data = cond2;
                    disp('cond2-pattern comp')
                end
            end
            
            % remove trials with artifacts
            chan_artifact_thresh = chan_artifact_thresh(allreg',:);
            bad_trials = [];
            for elec = 1:size(trial_data, 3)
                for trials = 1:size(trial_data,1)
                    trial_temp = trial_data(trials, :,elec);
                    
                    % if trail has artifact = nan
                    idx_above_thresh = (trial_temp)>chan_artifact_thresh(elec,2);
                    idx_below_thresh = (trial_temp)<chan_artifact_thresh(elec,1);
                    total_idx = idx_above_thresh+idx_below_thresh;
                    if sum(total_idx)>0 % if trial is artifact free
                        bad_trials = [bad_trials trials];
                    end
                end
            end
            bad_trials = unique(bad_trials);
            trial_data = trial_data(~ismember(1:size(trial_data,1),bad_trials),:,:);
            
            % Save subject;condition; reg1-reg2 data
            cd('/tmp/yassamri/iEEG/sandra/DataForDarmon')
            save([reg1_name '_' reg2_name '_' 'subj' num2str(sub_counter) 'Cond' num2str(cond_num) ], ...
            'reg1_name', 'reg2_name',  'length_reg1','length_reg2', 'trial_data')
        end
    end
end
            