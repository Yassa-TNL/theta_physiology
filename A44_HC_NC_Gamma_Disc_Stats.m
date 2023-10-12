% plot cluster determined group bar plots
clear all;close all;clc
fn_ext        = '_cue_responseive_yes'
baseline      = [ '_cond_spec_prestim' fn_ext]
lock     = 'onset' % onset % response
exp_type = 'tuning_correct'% 'encoding' 'tuning_correct'
ref      = 'LM' 
plot_4_conds =  ''
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
subj_list = {'39' '44' '57' '63' '66' '84' '85' '87'}
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock ])
fs=500

if strcmp('onset', lock)
    pre_stim = 0.5; post_stim = 2.0;
elseif strcmp('onset2', lock)
    pre_stim = 0.5; post_stim = 2.6;
end
cd(['/mnt/yassamri/iEEG/sandra/group_data/IndivSubjData/LM_reref/normalization_cond_spec_prestim/tuning_correct_onset/39' ])

load(['subj39_OFC_spectrograms_cond_spec_prestim_cue_responseive_yes.mat'])
% load cluster matrix and preview
if strcmp('onset', lock)
    if strcmp('encoding',exp_type)
        time_idx = [.3*fs 2.5*fs]
    else
        time_idx = [.3*fs 1.5*fs]
    end
elseif strcmp('response', lock)
    time_idx = [51 551];
end
minfreq       = 3;
maxfreq       = 200;
[freq] = get_freq(fs,minfreq, maxfreq);
zmapthresh_temp = zeros(length(freq), size(cond1,2));
zmapthresh_temp(freq>35,1.5*fs:2*fs)=1;
figure;imagesc(linspace(-pre_stim, post_stim,size(zmapthresh_temp,2)), 1:freq,zmapthresh_temp)


subj_NC_Gamma_specificty = cell(8,1);
for iSubj = 1:8
    
    [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
        ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj_list{iSubj})
    
    desired_sites = [OFC_chan_idx fro_chan_idx temp_chan_idx cingulate_chan_idx insula_chan_idx  EC_chan_idx ]; 
    NC_ROI_disc_cond2_cond3 = zeros(length (desired_sites),5);
    NC_ROI_disc_cond2_cond3(:,1) = desired_sites';

    elec_cntr     = 0;
    

    for regCounter = 1:6 % num of regions
            cd(['/mnt/yassamri/iEEG/sandra/group_data/IndivSubjData/LM_reref/normalization_cond_spec_prestim/tuning_correct_onset/' subj_list{iSubj} ])

        if regCounter ==1
            % OFC
            load(['subj' subj_list{iSubj} '_OFC_spectrograms' baseline '.mat'])
            chan_val = OFC_chan_idx;
        elseif regCounter==2
            % fro
            load(['subj' subj_list{iSubj} '_FRO_spectrograms' baseline '.mat'])
            chan_val = fro_chan_idx;
        elseif regCounter==3
            % temp
            load(['subj' subj_list{iSubj} '_TEMP_spectrograms' baseline '.mat'])
            chan_val = temp_chan_idx;
        elseif regCounter ==4
            % CING
            load(['subj' subj_list{iSubj} '_CING_spectrograms' baseline '.mat'])
            chan_val = cingulate_chan_idx;
        elseif regCounter==5
            % INS
            load(['subj' subj_list{iSubj} '_INS_spectrograms' baseline '.mat'])
            chan_val = insula_chan_idx;
        elseif regCounter==6
            % EC
            load(['subj' subj_list{iSubj} '_EC_spectrograms' baseline '.mat'])
            chan_val = EC_chan_idx;
%         elseif regCounter==7
%             % HC
%             load(['subj' subj_list{iSubj} '_HC_spectrograms_cond_spec_prestim_cue_responseive_yes.mat'])
%             chan_val = HC_chan_idx;
%         elseif regCounter==8
%             % CA3
%             load(['subj' subj_list{iSubj} '_CA3_spectrograms_cond_spec_prestim_cue_responseive_yes.mat'])
%             chan_val = CA3_chan_idx;
        end
        
    for iElec = 1:length(chan_val)
        elec_cntr = elec_cntr+1;
        % cond2
        current_elec_2 = cond2(:,:,:,iElec);
        %reset = nanmean(nanmean(current_elec_2(:,0.3*fs:0.5*fs,:),3),2);
        %current_elec_2 = current_elec_2
        cond2_trl_val  = nan(1,size(current_elec_2,3));
        for iTrl = 1:size(current_elec_2,3)
        temp = [];
        temp = current_elec_2(:,:,iTrl);
        cond2_trl_val(iTrl) = nanmean(temp(logical(zmapthresh_temp)));
        end
       
        
        % cond3
        current_elec_3 = cond3(:,:,:,iElec);
        %reset = nanmean(nanmean(current_elec_3(:,0.3*fs:0.5*fs,:),3),2);
        %current_elec_3 = current_elec_3 - reset;
        cond3_trl_val  = nan(1,size(current_elec_3,3));
        for iTrl = 1:size(current_elec_3,3)
        temp = [];
        temp = current_elec_3(:,:,iTrl);
        cond3_trl_val(iTrl) = nanmean(temp(logical(zmapthresh_temp)));
        end
      
        
        % remove nans
        cond2_trl_val = cond2_trl_val(~isnan(cond2_trl_val));
        cond3_trl_val = cond3_trl_val(~isnan(cond3_trl_val));
        
        % stats and dprime
        reps = 500;
        adata = cond3_trl_val;
        bdata = cond2_trl_val;
        p = permutation_unpaired_higher(adata, bdata, reps);
        dprime = nanmean(cond3_trl_val) - nanmean(cond2_trl_val)/mean([std(cond3_trl_val,0,2) std(cond2_trl_val,0,2)]);
        
        NC_ROI_disc_cond2_cond3(elec_cntr,2) = p;
        NC_ROI_disc_cond2_cond3(elec_cntr,3) = dprime;
        NC_ROI_disc_cond2_cond3(elec_cntr,4) = length(cond2_trl_val);
        NC_ROI_disc_cond2_cond3(elec_cntr,5) = length(cond3_trl_val);

%         % make bar plot and save
%         f=figure('visible','off');
%         bar_vector_mn1 = [nanmean(cond2_trl_val) nanmean(cond3_trl_val) ];
%         bar_vector_sem1 = [nanstd(cond2_trl_val,0,2)/sqrt(length(cond3_trl_val))...
%             nanstd(cond3_trl_val,0,2)/sqrt(length(cond2_trl_val))];
%         bar([bar_vector_mn1]);hold on
%         errorbar(1:length(bar_vector_mn1),bar_vector_mn1,bar_vector_sem1, 'rx')
%         set(gca, 'XTick', 1:length(bar_vector_mn1), 'XTickLabel', {'lure-', 'Lure+'},'XTickLabelRotation',45)
%         
%         title([subj_list{iSubj} ' reg' num2str(regCounter) ' elec' num2str(iElec) ' p' num2str(p)])
%         cd('/mnt/yassamri/iEEG/sandra/OFC_FRO_TEMP')
%         saveas(gcf,[[subj_list{iSubj} ' reg' num2str(regCounter) ' elec' num2str(iElec) ' p' num2str(p)] '.bmp'])
%         
%         % make traces and save
%         desired_freq_lo  = 4
%         desired_freq_hi  = 5
%         minfreq = 3; maxfreq = 200;
%         [freq]  = get_freq(fs,minfreq, maxfreq);
% 
%         trace2 = squeeze(nanmean(current_elec_2((freq>desired_freq_lo)&(freq<desired_freq_hi), 1:1.5*fs, :),1))';
%         reset = nanmean(nanmean(trace2(:,0.3*fs:0.5*fs),1));
%         trace2=trace2-reset;
% 
%         trace3 = squeeze(nanmean(current_elec_3((freq>desired_freq_lo)&(freq<desired_freq_hi), 1:1.5*fs, :),1))';
%         reset = nanmean(nanmean(trace3(:,0.3*fs:0.5*fs),1));
%         trace3=trace3-reset;
%         
%         g=figure('visible','off');
%         hold on;
%         
%         stdshade(trace3,.1,'m',linspace(-0.5,1.5, size(trace3,2)),[] ,[], []);
%         h1 = plot(linspace(-0.5,1.5, size(trace3,2)),nanmean(trace3,1), 'm', 'LineWidth', 2);
%         
%         stdshade(trace2,.1,'b',linspace(-0.5,1.5, size(trace2,2)),[] ,[], []);
%         h2 = plot(linspace(-0.5,1.5, size(trace2,2)),nanmean(trace2,1), 'b', 'LineWidth', 2);
%         line([0 0],ylim,'color','k')
%         saveas(gcf,[[subj_list{iSubj} ' reg' num2str(regCounter) ' elec' num2str(iElec) ' p' num2str(p)] 'tr.bmp'])


         clear current_elec_2  current_elec_3 temp  
    end
    end
    subj_NC_Gamma_specificty{iSubj}=NC_ROI_disc_cond2_cond3;
end

% save sig chans
cd('/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/LM_reref/tuning_correct_onset')
save('subj_ALL_NC_Gamma_specificty', 'subj_NC_Gamma_specificty') % ALL regions
