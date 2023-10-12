clear all; close all;clc
subj_list = {'39' '44'  '57'  '63' '66' '84' '85' '87'};%
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
group_roi_mni = cell(1,8)
i=0
roi_name = 'OFC_FRO_TEMP_CING_INS_EC'%'OFC_FRO_TEMP_CING_INS_EC' %'OFC_FRO_TEMP'
for iSubj = 1:length(subj_list)
    
    subj=subj_list{iSubj};
    
    % get cure resp elecs
    [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
        ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj)
    if strcmp('OFC_FRO_TEMP_CING_INS_EC',roi_name)
    rois = [OFC_chan_idx fro_chan_idx temp_chan_idx cingulate_chan_idx insula_chan_idx EC_chan_idx];
    elseif strcmp('OFC_FRO_TEMP',roi_name)
    rois = [OFC_chan_idx fro_chan_idx temp_chan_idx];    
    elseif strcmp('NC_HC_Gamma',roi_name)
    rois=[OFC_chan_idx fro_chan_idx temp_chan_idx cingulate_chan_idx insula_chan_idx  EC_chan_idx ];    
    end
    elec_pos_chan_roi = nan(length(rois),3);
 
 % load chan_label 
     cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
     if strcmp('84',subj) || strcmp('85',subj)|| strcmp('87',subj)
     load(['normalization_cond_spec_prestim_ref_LM_baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_' subj '_select_chan_3_clinical_fs_500.mat']) 
     else
     load(['normalization_cond_spec_prestim_ref_LM_baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_' subj '_select_chan_3_fs_500.mat']) 
     end
    
     % fix some error in chan label
     if strcmp('84',subj)
         chan_label{103} = 'ASI2';
     end
     
     % load mni info
     cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/FT_Pipeline/Electrodes'])
     load(['IR' subj '_elec_mni_frv.mat'])

    % define position and mni variable
    chan_label_mni =   elec_mni_frv.label(:);
    elec_pos       = elec_mni_frv.elecpos;

    % align chan_labels, and get mni for only chns of int
    if ~strcmp('84',subj) && ~strcmp('85',subj) && ~strcmp('87',subj)
        for chan = 1:length(chan_label)
            chan_label{chan} = chan_label{chan}(4:end-3);
        end
    end

    chan_label_roi = chan_label(rois);

        for iChan_SG = 1:length(chan_label_roi)
            for iChan_mni = 1:length(chan_label_mni)
                if strcmp(chan_label_mni{iChan_mni}, chan_label_roi{iChan_SG})
                    elec_pos_chan_roi(iChan_SG,:) =elec_pos(iChan_mni,:);
                end
            end
        end

      group_roi_mni{iSubj}=elec_pos_chan_roi;


end

cd('/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/LM_reref/tuning_correct_onset')
if strcmp('OFC_FRO_TEMP_CING_INS_EC',roi_name)
    load('subj_ALL_NC_specificty')
    subj_site_specificty = cat(1,subj_ALL_NC_specificty{:});
elseif strcmp('OFC_FRO_TEMP',roi_name)
    load('subj_NC_specificty.mat')
    subj_site_specificty = cat(1,subj_NC_specificty{:});
elseif strcmp('NC_HC_Gamma',roi_name)
    load('subj_ALL_NC_HC_Gamma_specificty.mat')
    subj_site_specificty=cat(1,subj_NC_Gamma_specificty{:});
elseif strcmp('NC_Gamma',roi_name)
    load('subj_ALL_NC_Gamma_specificty.mat')
    subj_site_specificty=cat(1,subj_NC_Gamma_specificty{:});
end

group_roi_mni = cat(1,group_roi_mni{:})

% plot sig chans d prime
dprime = subj_site_specificty(:,3)
pvalue =subj_site_specificty(:,2)
sig_lureplus_disc_chans = pvalue<0.05 & dprime>0
imagesc(dprime(sig_lureplus_disc_chans));colorbar;colormap jet
%%
group_roi_mni_sig = group_roi_mni(sig_lureplus_disc_chans,:)
group_roi_mni_sig = group_roi_mni_sig(~isnan(group_roi_mni_sig(:,1)),:)

filler = cell(size(group_roi_mni_sig,1),1)
for i = 1:size(filler,1)
  filler{i} = [roi_name num2str(i) ]
end

% to do manually: make a text file out of data_table and put in e_loc folder
Dat_Table =  table(filler,group_roi_mni_sig(:,1),group_roi_mni_sig(:,2),group_roi_mni_sig(:,3));
magenta_region_label  = zeros(size(group_roi_mni_sig));% cue resp
magenta_region_label(:,[1 3]) = 255;

