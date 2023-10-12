% description: script which will generate the null distribution for a given
% condition, subject, and region. Will shuffle the phase data of the sender
% region. Will save ch1 -> ch2, ch2 -> ch1, and normalized ch->ch2

clear all; close all; clc
format shortg
cond_nums = 1; % condition to generate null distribut
% encoding lure+  = 1,  encoding lure-  = 2
% retrieval lure+ = 3 , retrieval lure-  = 2

%do one area and subj at a time bc takes time.
reg1_name   = 'TEMP'
subj_list = {'87'}% '39' '44'  '57'  '63' '66' '84' '85' '87'
phase       = 'encoding'  %'encoding' 'retrieval'

% relatively fixed params
nPerms      =1000;
fn_ext      = 'cue_resp' %
lock        = 'onset'   %'onset'    'response'
freq_range  = 'deltatheta'
fpass_list  = [4;5]';
reg2_name   = 'HC';
sig_elecs   = ''%'OFC_FRO_TEMP_CING_INS_EC'
addpath('/tmp/yassamri/iEEG/sandra/analysis_pipeline_final')
cd('/tmp/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/LM_reref/tuning_correct_onset')
if strcmp('OFC_FRO_TEMP_CING_INS_EC',sig_elecs)
load('subj_ALL_NC_specificty')
elseif strcmp('OFC_FRO_TEMP',sig_elecs)
load('subj_NC_specificty.mat')
elseif strcmp('HC_NC_Gamma',sig_elecs)
load('subj_ALL_NC_HC_Gamma_specificty.mat')
end

if strcmp('encoding',phase)
    %cond_nums = 1:2;
     %       time1_list = 0:0.01:1.5;
      %      time2_list = 0.5:0.01:2;
    %         time1_list = round([0:0.05:1.5],2);
    %         time2_list = round([0.5:0.05:2],2);
        time1_list = 0;
        time2_list = 2;

elseif strcmp('retrieval',phase)&&strcmp('onset',lock)
 %   cond_nums = 2:3;
%     time1_list = 0:0.01:0.5;
%     time2_list = 0.5:0.01:1;
    
    time1_list = 0;
    time2_list = 1;
    
    
%     time1_list = 1
%     time2_list = 1.5

    %         time1_list = round([0:0.05:0.5],2) % 90% overlap
    %         time2_list = round([0.5:0.05:1],2)
    
    %         time1_list = round([0:0.25:0.5],2)
    %         time2_list = round([0.5:0.25:1],2)
    
elseif strcmp('retrieval',phase)&&strcmp('response',lock)
    %cond_nums = 2:3;
    
    
    %     time1_list = round([0.5:0.25:1],2)
    %     time2_list = round([1:0.25:1.5],2)
    %
    %     time1_list = round([0.5:0.125:1.25],2)
    %     time2_list = round([0.75:0.125:1.5],2)
    time1_list = 0.5;
    time2_list = 1.5;
end
store_idx_strt = [];
store_idx_ed   = [];
for iTime = 1:length(time1_list)
    
    for cond_num    = cond_nums
        time1       = time1_list(iTime); % response is 0.5 to 1.5 to get 1 sec be4 resp
        time2       = time2_list(iTime);
        if strcmp('encoding',phase)
            exp_type    = phase; %'tuning_correct' 'encoding'
        else
            exp_type    = 'tuning_correct'; %'tuning_correct' 'encoding'
        end
        DS          = 'yes';
        ref         = 'LM';
        clinical    = 'yes';
        lag         = '';
        
      %  subj_list = {'39' '44'  '57'  '63' '66' '84' '85' '87'}
        
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
            
            if strcmp('OFC_FRO_TEMP_CING_INS_EC',sig_elecs)
                sig_chan_log_vec = logical(subj_ALL_NC_specificty{sub_counter}(:,2)<0.05);
                reg1 = reg1(ismember(reg1, subj_ALL_NC_specificty{sub_counter}(sig_chan_log_vec,1)));
            elseif strcmp('OFC_FRO_TEMP',sig_elecs)
                sig_chan_log_vec = logical(subj_NC_specificty{sub_counter}(:,2)<0.05);
                reg1 = reg1(ismember(reg1, subj_NC_specificty{sub_counter}(sig_chan_log_vec,1)));
                
            elseif strcmp('HC_NC_Gamma',sig_elecs)
                dprime = subj_HC_NC_specificty{sub_counter}(:,3);
                pvalue =subj_HC_NC_specificty{sub_counter}(:,2);
                sig_lureplus_disc_chans = pvalue<0.05 & dprime>0;
                reg1 = reg1(ismember(reg1, subj_HC_NC_specificty{sub_counter}(sig_lureplus_disc_chans,1)));
                %reg2 = reg2(ismember(reg2, subj_HC_NC_specificty{sub_counter}(sig_lureplus_disc_chans,1)));
            end
            
            allreg = [reg1 reg2];
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
            
            for iFpass = 1:size(fpass_list,1)
                fpass = [fpass_list(iFpass,:)];
                
                %% Get phase MTX
                phase_Data = nan(size(trial_data,2), size(trial_data,1), size(trial_data,3)); % time X trial X chan
                for iChan = 1:size(trial_data,3)
                    % get chan data, then band pass it, then hilbert, then get angle
                    phase_Data(:,:,iChan) = angle(hilbert(bandpass(trial_data (:,:,iChan)', fpass,fs)));
                end
                
                

                %% calculcate from a --> b & b --> a
                if strcmp('', lag)
                    delta_val = round(.1*fs);
                    ch1_to_ch2 =nan(nPerms,length(reg1) * length(length(reg1)+1:length(allreg)));
                    ch2_to_ch1 =nan(nPerms,length(reg1) * length(length(reg1)+1:length(allreg)));
                    PTE_ch1_to_ch2_norm =nan(1,length(reg1) * length(length(reg1)+1:length(allreg)));
                    
                elseif strcmp('yes', lag)
                    delta_val = round((0.005:.005:.2)*fs);
                    ch1_to_ch2 =nan(length(reg1) * length(length(reg1)+1:length(allreg)), length(delta_val));
                    ch2_to_ch1 =nan(length(reg1) * length(length(reg1)+1:length(allreg)), length(delta_val));
                    PTE_ch1_to_ch2_norm = nan(length(reg1) * length(length(reg1)+1:length(allreg)), length(delta_val));
                end
                
                d_cntr = 0;
                for delta =  delta_val
                    d_cntr=d_cntr+1;
                    cntr = 0;
                    
                    for iElec1 = 1:length(reg1) % Reg1 elecs
                        for iElec2 = length(reg1)+1:length(allreg) % Reg2 elecs
                            cntr = cntr+1;
                            ch1_phase = phase_Data(:,:,iElec1);  % input process
                            ch2_phase = phase_Data(:,:,iElec2);  % output process

                            for iPerm = 1:nPerms

                            %% a -> b
                            % organizing the data: assume we have only two channels each of which has n trials.
                            ch1_phaseShuffle = nan(size(ch1_phase));
                            
                            % shuffle one of the data
                            for iTrial=1:size(ch1_phase,2)
                                PermIdx = randperm(size(ch1_phase,1));
                                ch1_phaseShuffle(:,iTrial)=ch1_phase(PermIdx,iTrial);
                            end
                            
                            
                            % estimating the probability mass functions involved in TE definition
                            ch1_bin_width = 3.5*circ_std(ch1_phaseShuffle(:))/length(ch1_phaseShuffle(:))^(1/3); %  used in the original PTE paper.
                            ch2_bin_width = 3.5*circ_std(ch2_phase(:))/length(ch2_phase(:))^(1/3); %  used in the original PTE paper.
                            
                            % estimating pmf of output phase(t - delta)
                            ch2_pst = ch2_phase(1:end-delta,:);
                            ch2_phase_pmf = histcnd(ch2_pst(:), ...
                                linspace(-pi,pi,round(2*pi/ch2_bin_width)))/length(ch2_pst(:));
                            
                            % estimating joint pmf of output phase past(1:end-delta) phase and
                            % output current(1+delta:end) phase
                            ch2_phase_pairs = [reshape(ch2_phase(1:end-delta,:),numel(ch2_phase(1:end-delta,:)),1), ...
                                reshape(ch2_phase(1+delta:end,:),numel(ch2_phase(1+delta:end,:)),1)];
                            ch2_phase_pairs_pmf   = histcnd(ch2_phase_pairs(:,1),ch2_phase_pairs(:,2), linspace(-pi,pi,round(2*pi/ch2_bin_width)),linspace(-pi,pi,round(2*pi/ch2_bin_width)));
                            ch2_phase_pairs_pmf   = ch2_phase_pairs_pmf/size(ch2_phase_pairs,1);
                            
                            % estimating joint pmf of input and output current
                            ch2_ch1_crnt_phase_pairs = [reshape(ch2_phase(1:end-delta,:),numel(ch2_phase(1:end-delta,:)),1), ...
                                reshape(ch1_phaseShuffle(1:end-delta,:),numel(ch1_phaseShuffle(1:end-delta,:)),1)];
                            ch1_ch2_phase_pairs_pmf = histcnd(ch2_ch1_crnt_phase_pairs(:,2),ch2_ch1_crnt_phase_pairs(:,1), ...
                                linspace(-pi,pi,round(2*pi/ch1_bin_width)), linspace(-pi,pi,round(2*pi/ch2_bin_width)))/length(ch2_ch1_crnt_phase_pairs(:,2));
                            
                            % estimating joint pmf of input phase(past) and output phase(past) and output(current)
                            x1pst  = reshape(ch1_phaseShuffle(1:end-delta,:),numel(ch1_phaseShuffle(1:end-delta,:)),1);
                            x2pst  = reshape(ch2_phase(1:end-delta,:),numel(ch2_phase(1:end-delta,:)),1);
                            x2crnt = reshape(ch2_phase(1+delta:end,:),numel(ch2_phase(1+delta:end,:)),1);
                            ch1_ch2_phase_triplets_pmf = histcnd(x1pst, x2pst, x2crnt, ...
                                linspace(-pi,pi,round(2*pi/ch1_bin_width)),linspace(-pi,pi,round(2*pi/ch2_bin_width)),linspace(-pi,pi,round(2*pi/ch2_bin_width)));
                            ch1_ch2_phase_triplets_pmf = ch1_ch2_phase_triplets_pmf/length(x1pst);
                            
                            % Now that we have all the required pmf's, we can calculate the entropies involved in the formula of PTE
                            H_ch2_crnt     = -nansum(nansum(ch2_phase_pmf .* log2(ch2_phase_pmf)));
                            H_ch2_pst_crnt = -nansum(nansum(ch2_phase_pairs_pmf .* log2(ch2_phase_pairs_pmf)));
                            H_ch2_ch1_crnt = -nansum(nansum(ch1_ch2_phase_pairs_pmf .* log2(ch1_ch2_phase_pairs_pmf)));
                            H_ch2_pst_ch1_pst_ch2_crnt = -nansum(nansum(nansum(ch1_ch2_phase_triplets_pmf .* log2(ch1_ch2_phase_triplets_pmf))));
                            PTE_ch1_to_ch2a = H_ch2_pst_crnt + H_ch2_ch1_crnt - H_ch2_crnt - H_ch2_pst_ch1_pst_ch2_crnt;
                            
                            %% from b to a
                            % organizing the data: assume we have only two channels each of which has n trials.
                            ch1_phase = phase_Data(:,:,iElec2); % input process
                            ch1_phaseShuffle = nan(size(ch1_phase));
                            
                            % shuffle one of the data
                            for iTrial=1:size(ch1_phase,2)
                                PermIdx = randperm(size(ch1_phase,1));
                                ch1_phaseShuffle(:,iTrial)=ch1_phase(PermIdx,iTrial);
                            end
                            
                            ch2_phase = phase_Data(:,:,iElec1);  % output process
                            
                            ch1_bin_width = 3.5*circ_std(ch1_phaseShuffle(:))/length(ch1_phaseShuffle(:))^(1/3); %  used in the original PTE paper.
                            ch2_bin_width = 3.5*circ_std(ch2_phase(:))/length(ch2_phase(:))^(1/3); %  used in the original PTE paper.
                            
                            % estimating pmf of output phase(t - delta)
                            ch2_pst = ch2_phase(1:end-delta,:);
                            ch2_phase_pmf = histcnd(ch2_pst(:), ...
                                linspace(-pi,pi,round(2*pi/ch2_bin_width)))/length(ch2_pst(:));
                            
                            % estimating joint pmf of output phase past(1:end-delta) phase and
                            % output current(1+delta:end) phasw
                            ch2_phase_pairs = [reshape(ch2_phase(1:end-delta,:),numel(ch2_phase(1:end-delta,:)),1), ...
                                reshape(ch2_phase(1+delta:end,:),numel(ch2_phase(1+delta:end,:)),1)];
                            ch2_phase_pairs_pmf   = histcnd(ch2_phase_pairs(:,1),ch2_phase_pairs(:,2), linspace(-pi,pi,round(2*pi/ch2_bin_width)),linspace(-pi,pi,round(2*pi/ch2_bin_width)));
                            ch2_phase_pairs_pmf   = ch2_phase_pairs_pmf/size(ch2_phase_pairs,1);
                            
                            % estimating joint pmf of input and output current
                            ch2_ch1_crnt_phase_pairs = [reshape(ch2_phase(1:end-delta,:),numel(ch2_phase(1:end-delta,:)),1), ...
                                reshape(ch1_phaseShuffle(1:end-delta,:),numel(ch1_phaseShuffle(1:end-delta,:)),1)];
                            ch1_ch2_phase_pairs_pmf = histcnd(ch2_ch1_crnt_phase_pairs(:,2),ch2_ch1_crnt_phase_pairs(:,1), ...
                                linspace(-pi,pi,round(2*pi/ch1_bin_width)), linspace(-pi,pi,round(2*pi/ch2_bin_width)))/length(ch2_ch1_crnt_phase_pairs(:,2));
                            
                            % estimating joint pmf of input phase(past) and output phase(past) and output(current)
                            x1pst  = reshape(ch1_phaseShuffle(1:end-delta,:),numel(ch1_phaseShuffle(1:end-delta,:)),1);
                            x2pst  = reshape(ch2_phase(1:end-delta,:),numel(ch2_phase(1:end-delta,:)),1);
                            x2crnt = reshape(ch2_phase(1+delta:end,:),numel(ch2_phase(1+delta:end,:)),1);
                            ch1_ch2_phase_triplets_pmf = histcnd(x1pst, x2pst, x2crnt, ...
                                linspace(-pi,pi,round(2*pi/ch1_bin_width)),linspace(-pi,pi,round(2*pi/ch2_bin_width)),linspace(-pi,pi,round(2*pi/ch2_bin_width)));
                            ch1_ch2_phase_triplets_pmf = ch1_ch2_phase_triplets_pmf/length(x1pst);
                            
                            % Now that we have all the required pmf's, we can calculate the entropies involved in the formula of PTE
                            H_ch2_crnt     = -nansum(nansum(ch2_phase_pmf .* log2(ch2_phase_pmf)));
                            H_ch2_pst_crnt = -nansum(nansum(ch2_phase_pairs_pmf .* log2(ch2_phase_pairs_pmf)));
                            H_ch2_ch1_crnt = -nansum(nansum(ch1_ch2_phase_pairs_pmf .* log2(ch1_ch2_phase_pairs_pmf)));
                            H_ch2_pst_ch1_pst_ch2_crnt = -nansum(nansum(nansum(ch1_ch2_phase_triplets_pmf .* log2(ch1_ch2_phase_triplets_pmf))));
                            PTE_ch1_to_ch2b = H_ch2_pst_crnt + H_ch2_ch1_crnt - H_ch2_crnt - H_ch2_pst_ch1_pst_ch2_crnt;
                            
                            if strcmp('', lag)
                                ch1_to_ch2 (iPerm,cntr) = PTE_ch1_to_ch2a ;
                                ch2_to_ch1 (iPerm,cntr) = PTE_ch1_to_ch2b ;
                                PTE_ch1_to_ch2_norm (iPerm,cntr) = ch1_to_ch2 (iPerm,cntr) / (ch1_to_ch2 (iPerm,cntr) + ch2_to_ch1 (iPerm,cntr));
%                             else
%                                 ch1_to_ch2 (cntr, d_cntr) = PTE_ch1_to_ch2a ;
%                                 ch2_to_ch1 (cntr, d_cntr) = PTE_ch1_to_ch2b ;
%                                 PTE_ch1_to_ch2_norm (cntr, d_cntr) = ch1_to_ch2 (cntr,d_cntr) / (ch1_to_ch2 (cntr,d_cntr) + ch2_to_ch1 (cntr,d_cntr));
                             end
                        end
                    end
                    
                end
                end
                
                cd (['/tmp/yassamri/iEEG/sandra/PTE_results/' fn_ext '/' freq_range '/' phase '/' lock '/'])
                if ~isdir([num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1) '_' num2str(time2) 'sec'])
                    mkdir([num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1) '_' num2str(time2) 'sec'])
                end
                
                cd([num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1) '_' num2str(time2) 'sec'])
                save([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond' num2str(cond_num) sig_elecs 'Chance'], 'PTE_ch1_to_ch2_norm', 'ch1_to_ch2', 'ch2_to_ch1', 'delta_val')
            end
        end
    end
    store_idx_strt = [store_idx_strt strt_idx];
    store_idx_ed   = [store_idx_ed end_idx];
end
disp('done')