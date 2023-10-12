function save_individual_subj_region_data(cue_responsive,ref,baseline, exp_type, lock, subj,removeERP, ...
    norm_freq_acrs_chan_cond_1, norm_freq_acrs_chan_cond_2,norm_freq_acrs_chan_cond_3,norm_freq_acrs_chan_cond_4)

% [will save an  spec per chan for each sub when run thru A003_condition_wise_spectral_analysis.m]
% output directory: /mnt/yassamri/iEEG/sandra/group_data/IndivSubjData
if strcmp('yes',cue_responsive)
    [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
        ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj);
else
    [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
        ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean(subj);
end
cd(['/mnt/yassamri/iEEG/sandra/group_data/IndivSubjData/' ref '_reref/normalization_' baseline '/' ...
    exp_type '_' lock '/' subj])

fnExtension = [baseline '_cue_responseive_' cue_responsive 'removeERP' removeERP];

cond1 = norm_freq_acrs_chan_cond_1(:,:,:,OFC_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,OFC_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,OFC_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,OFC_chan_idx);
save(['subj' subj '_OFC_spectrograms_' fnExtension],'cond1','cond2','cond3','cond4','-v7.3')

clear cond1 cond2 cond3 cond4
cond1 = norm_freq_acrs_chan_cond_1(:,:,:,fro_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,fro_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,fro_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,fro_chan_idx);
save(['subj' subj '_FRO_spectrograms_' fnExtension],'cond1','cond2','cond3','cond4','-v7.3')

clear cond1 cond2 cond3 cond4
cond1 = norm_freq_acrs_chan_cond_1(:,:,:,temp_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,temp_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,temp_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,temp_chan_idx);
save(['subj' subj '_TEMP_spectrograms_' fnExtension], 'cond1','cond2','cond3','cond4','-v7.3')

clear cond1 cond2 cond3 cond4
cond1 = norm_freq_acrs_chan_cond_1(:,:,:,insula_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,insula_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,insula_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,insula_chan_idx);
save(['subj' subj '_INS_spectrograms_' fnExtension], 'cond1','cond2','cond3','cond4','-v7.3')

clear cond1 cond2 cond3 cond4
cond1 = norm_freq_acrs_chan_cond_1(:,:,:,cingulate_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,cingulate_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,cingulate_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,cingulate_chan_idx);
save(['subj' subj '_CING_spectrograms_' fnExtension], 'cond1','cond2','cond3','cond4','-v7.3')

clear cond1 cond2 cond3 cond4
cond1 = norm_freq_acrs_chan_cond_1(:,:,:,EC_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,EC_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,EC_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,EC_chan_idx);
save(['subj' subj '_EC_spectrograms_' fnExtension], 'cond1','cond2','cond3','cond4','-v7.3')

clear cond1 cond2 cond3 cond4
cond1 = norm_freq_acrs_chan_cond_1(:,:,:,CA3_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,CA3_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,CA3_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,CA3_chan_idx);
save(['subj' subj '_CA3_spectrograms_' fnExtension], 'cond1','cond2','cond3','cond4','-v7.3')

clear cond1 cond2 cond3 cond4
cond1 = norm_freq_acrs_chan_cond_1(:,:,:,CA1_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,CA1_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,CA1_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,CA1_chan_idx);
save(['subj' subj '_CA1_spectrograms_' fnExtension], 'cond1','cond2','cond3','cond4','-v7.3')

clear cond1 cond2 cond3 cond4
cond1 = norm_freq_acrs_chan_cond_1(:,:,:,HC_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,HC_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,HC_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,HC_chan_idx);
save(['subj' subj '_HC_spectrograms_' fnExtension],'cond1','cond2','cond3','cond4','-v7.3')

clear cond1 cond2 cond3 cond4
cond1 = norm_freq_acrs_chan_cond_1(:,:,:,NC_chan_idx);
cond2 = norm_freq_acrs_chan_cond_2(:,:,:,NC_chan_idx);
cond3 = norm_freq_acrs_chan_cond_3(:,:,:,NC_chan_idx);
cond4 = norm_freq_acrs_chan_cond_4(:,:,:,NC_chan_idx);
save(['subj' subj '_NC_spectrograms_' fnExtension], 'cond1','cond2','cond3','cond4','-v7.3')

end