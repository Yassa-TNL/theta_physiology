clear all;close all;clc
subj_list     = {'87'}
iSubj         = 1
subj          = subj_list{iSubj}
phase_list    = { 'tuning_correct'}; %encoding, tuning_correc
reg_list      = {'OFC' 'FRO' 'TEMP' 'INS' 'CING'  'EC' 'HC' 'CA3'}
lock          = 'onset'   %onset response
ref           = 'LM'
fs            = 500;
baseline      = '_cond_spec_prestim' % '': entire recording, '_prestim'
desired_freq  = 'deltatheta'
desired_freq_lo  = 2
desired_freq_hi  = 6
trial_count_limit = ''
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
minfreq = 3; maxfreq = 200;
[freq]  = get_freq(fs,minfreq, maxfreq);
tickmarks = 1:21:length(freq);
sig_elecs_per_reg = cell(1,length(reg_list));

%%
for iReg = 1:length(reg_list)
    reg = reg_list{iReg};
    
    % get chan length
[OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
    ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean(subj)

if strcmp('OFC',reg); 
    chan_idx = OFC_chan_idx; 
elseif strcmp('FRO',reg); 
    chan_idx = fro_chan_idx;
elseif strcmp('TEMP',reg)
    chan_idx = temp_chan_idx;
elseif strcmp('CING',reg);
    chan_idx = cingulate_chan_idx;
elseif strcmp('INS',reg)
    chan_idx = insula_chan_idx;
elseif strcmp('EC',reg);
    chan_idx = EC_chan_idx;
elseif strcmp('CA3',reg)
    chan_idx = CA3_chan_idx;
elseif strcmp('HC',reg);
    chan_idx = HC_chan_idx;
end

    chan_sig = zeros(1,length(chan_idx));

    for iPhase = 1:length(phase_list)
        exp_type      = phase_list{iPhase}; %encoding, tuning_correct
        if strcmp('onset', lock)
            pre_stim = 0.5;
            if strcmp('encoding', exp_type)
                des_post_stim_dur=2;
                cond_num=2;
            else
                des_post_stim_dur=1;
                cond_num=4;
            end
        end
        cd(['/mnt/yassamri/iEEG/sandra/group_data/IndivSubjData/' ref '_reref/normalization_' baseline(2:end) '/' exp_type '_' lock '/' subj])
        
        if isfile(['subj' subj '_' reg '_spectrograms_cond_spec_prestim_cue_responseive_.mat']) 
         load(['subj' subj '_' reg '_spectrograms_cond_spec_prestim_cue_responseive_.mat'])
        else
            load(['subj' subj '_' reg '_spectrograms_cond_spec_prestim.mat'])
        end
        for iElec = 1:length(chan_idx)
            trace2 = squeeze(nanmean(cond2((freq>desired_freq_lo)&(freq<desired_freq_hi), :, :,iElec),1))'; % trial by time
            if strcmp('encoding',exp_type)
                trace3 = squeeze(nanmean(cond1((freq>desired_freq_lo)&(freq<desired_freq_hi), :, :,iElec),1))';
                trace1 = [];trace4=[];
            else
                trace1 = squeeze(nanmean(cond1((freq>desired_freq_lo)&(freq<desired_freq_hi), :, :,iElec),1))'; % trial by time
                trace3 = squeeze(nanmean(cond3((freq>desired_freq_lo)&(freq<desired_freq_hi), :, :,iElec),1))'; % trial by time
                trace4 = squeeze(nanmean(cond4((freq>desired_freq_lo)&(freq<desired_freq_hi), :, :,iElec),1))'; % trial by time
                %            remove nans
                trace1 = trace1(~isnan(trace1(:,1)),:);
                trace4 = trace4(~isnan(trace4(:,1)),:);
            end
            
            %            remove nans
            trace2 = trace2(~isnan(trace2(:,1)),:);
            trace3 = trace3(~isnan(trace3(:,1)),:);
            
            if strcmp('encoding',exp_type)
                
                conda = nanmean(trace2(:,0.2*fs:pre_stim*fs),2);
                condb = nanmean(trace2(:,(pre_stim*fs)+1:(pre_stim+des_post_stim_dur)*fs),2);
                reps = 1000; p1 = permutation_paired(conda, condb, reps);
                
                clear conda condb
                conda = nanmean(trace3(:,0.2*fs:pre_stim*fs),2);
                condb = nanmean(trace3(:,(pre_stim*fs)+1:(pre_stim+des_post_stim_dur)*fs),2);
                reps = 1000; p2 = permutation_paired(conda, condb, reps);
                
                all_conds = [trace2; trace3];
                conda = nanmean(all_conds(:,0.2*fs:pre_stim*fs),2);
                condb = nanmean(all_conds(:,(pre_stim*fs)+1:(pre_stim+des_post_stim_dur)*fs),2);
                reps = 1000; p5 = permutation_paired(conda, condb, reps);
                
                p3=100;p4=100;
            elseif strcmp('tuning_correct',exp_type)
                 clear conda condb
                conda = nanmean(trace1(:,0.2*fs:pre_stim*fs),2);
                condb = nanmean(trace1(:,(pre_stim*fs)+1:(pre_stim+des_post_stim_dur)*fs),2);
                reps = 1000; p1 = permutation_paired(conda, condb, reps);
                
                 clear conda condb
                conda = nanmean(trace2(:,0.2*fs:pre_stim*fs),2);
                condb = nanmean(trace2(:,(pre_stim*fs)+1:(pre_stim+des_post_stim_dur)*fs),2);
                reps = 1000; p2 = permutation_paired(conda, condb, reps);
                
                clear conda condb
                conda = nanmean(trace3(:,0.2*fs:pre_stim*fs),2);
                condb = nanmean(trace3(:,(pre_stim*fs)+1:(pre_stim+des_post_stim_dur)*fs),2);
                reps = 1000; p3 = permutation_paired(conda, condb, reps);
                
                clear conda condb
                conda = nanmean(trace4(:,0.2*fs:pre_stim*fs),2);
                condb = nanmean(trace4(:,(pre_stim*fs)+1:(pre_stim+des_post_stim_dur)*fs),2);
                reps = 1000; p4 = permutation_paired(conda, condb, reps);
                
                clear conda condb
                all_conds = [trace1; trace2; trace3; trace4];
                conda = nanmean(all_conds(:,0.2*fs:pre_stim*fs),2);
                condb = nanmean(all_conds(:,(pre_stim*fs)+1:(pre_stim+des_post_stim_dur)*fs),2);
                reps = 1000; p5 = permutation_paired(conda, condb, reps);
            end

            % save elec sig idx
            if p1<0.05 || p2<0.05 || p3<0.05 || p4<0.05 || p5<0.05
                chan_sig (iElec) =  1;
            else
                chan_sig (iElec) =  0;
            end
            
        end
    end
    sig_elecs_per_reg {iReg} = chan_sig;
end

cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
save('sig_elecs_per_reg', 'sig_elecs_per_reg')
disp('done')