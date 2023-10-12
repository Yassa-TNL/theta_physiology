function save_data_for_average_spectrograms_all_trials_elecs(cue_responsive,ref, baseline, exp_type, lock, subj, norm_freq_acrs_chan_cond_1,norm_freq_acrs_chan_cond_2,...
    norm_freq_acrs_chan_cond_3, norm_freq_acrs_chan_cond_4)

% [will save region specific matrices containing all trials and elecs for each sub when run thru A003_condition_wise_spectral_analysis.m]
% output directory: /mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms
if strcmp('yes',cue_responsive)
    [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
        ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj);
else
    [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
        ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean(subj);
end

OFC_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,OFC_chan_idx),4);
OFC_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,OFC_chan_idx),4);
OFC_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,OFC_chan_idx),4);
OFC_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,OFC_chan_idx),4);

fro_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,fro_chan_idx),4);
fro_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,fro_chan_idx),4);
fro_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,fro_chan_idx),4);
fro_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,fro_chan_idx),4);

temp_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,temp_chan_idx),4);
temp_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,temp_chan_idx),4);
temp_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,temp_chan_idx),4);
temp_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,temp_chan_idx),4);

ins_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,insula_chan_idx),4);
ins_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,insula_chan_idx),4);
ins_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,insula_chan_idx),4);
ins_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,insula_chan_idx),4);

cing_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,cingulate_chan_idx),4);
cing_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,cingulate_chan_idx),4);
cing_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,cingulate_chan_idx),4);
cing_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,cingulate_chan_idx),4);

ACC_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,ACC_chan_idx),4);
ACC_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,ACC_chan_idx),4);
ACC_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,ACC_chan_idx),4);
ACC_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,ACC_chan_idx),4);

EC_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,EC_chan_idx),4);
EC_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,EC_chan_idx),4);
EC_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,EC_chan_idx),4);
EC_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,EC_chan_idx),4);

CA3_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,CA3_chan_idx),4);
CA3_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,CA3_chan_idx),4);
CA3_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,CA3_chan_idx),4);
CA3_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,CA3_chan_idx),4);

CA1_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,CA1_chan_idx),4);
CA1_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,CA1_chan_idx),4);
CA1_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,CA1_chan_idx),4);
CA1_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,CA1_chan_idx),4);

HC_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,HC_chan_idx),4);
HC_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,HC_chan_idx),4);
HC_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,HC_chan_idx),4);
HC_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,HC_chan_idx),4);

MTL_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,MTL_chan_idx),4);
MTL_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,MTL_chan_idx),4);
MTL_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,MTL_chan_idx),4);
MTL_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,MTL_chan_idx),4);

NC_cond1 = nanmean(norm_freq_acrs_chan_cond_1(:,:,:,NC_chan_idx),4);
NC_cond2 = nanmean(norm_freq_acrs_chan_cond_2(:,:,:,NC_chan_idx),4);
NC_cond3 = nanmean(norm_freq_acrs_chan_cond_3(:,:,:,NC_chan_idx),4);
NC_cond4 = nanmean(norm_freq_acrs_chan_cond_4(:,:,:,NC_chan_idx),4);

cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock  ])
save(['subj' subj 'spectrograms_all_trials_across_channels_' baseline '_cue_responseive' cue_responsive], ...
    'OFC_cond1','OFC_cond2','OFC_cond3','OFC_cond4',...
    'fro_cond1','fro_cond2','fro_cond3','fro_cond4',...
    'temp_cond1','temp_cond2','temp_cond3','temp_cond4',...
    'ins_cond1','ins_cond2','ins_cond3','ins_cond4',...
    'cing_cond1','cing_cond2','cing_cond3', 'cing_cond4',...
    'ACC_cond1','ACC_cond2','ACC_cond3', 'ACC_cond4',...
    'EC_cond1','EC_cond2','EC_cond3', 'EC_cond4',...
    'CA3_cond1','CA3_cond2','CA3_cond3','CA3_cond4',...
    'CA1_cond1','CA1_cond2','CA1_cond3','CA1_cond4', ...
    'HC_cond1','HC_cond2','HC_cond3','HC_cond4',...
    'MTL_cond1','MTL_cond2','MTL_cond3','MTL_cond4', ...
    'NC_cond1','NC_cond2','NC_cond3','NC_cond4')

end

