clear all; close all;clc
subj_chan_name = cell(1600,1);
X = nan(1600,1);
Y = nan(1600,1);
Z = nan(1600,1);
anat_selec_index_group = nan(1600,1);
cue_resp_index_group = nan(1600,1);

subj_list = {'39' '44'  '57'  '63' '66' '84' '85' '87'};
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')

i=0
for iSubj = 1:length(subj_list)

     subj=subj_list{iSubj};
     
     % get anatomical elecs
[OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
    ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean(subj);
 anatomical_selec =  [OFC_chan_idx    fro_chan_idx   temp_chan_idx     insula_chan_idx    cingulate_chan_idx   HC_chan_idx EC_chan_idx  ];
      
 clear OFC_chan_idx    fro_chan_idx   temp_chan_idx     insula_chan_idx    cingulate_chan_idx   HC_chan_idx EC_chan_idx
 
 % get cure resp elecs
[OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
    ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj)
 cue_resp_elec =  [OFC_chan_idx    fro_chan_idx   temp_chan_idx     insula_chan_idx    cingulate_chan_idx   HC_chan_idx EC_chan_idx  ];
 
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
        chan_label_mni = elec_mni_frv.label(:);
        elec_pos       = elec_mni_frv.elecpos;
    
    % align chan_labels, and get mni for only chns of int
    if ~strcmp('84',subj) && ~strcmp('85',subj) && ~strcmp('87',subj)
     for chan = 1:length(chan_label)
        chan_label{chan} = chan_label{chan}(4:end-3);
     end
    end
    
    if ~strcmp('84',subj) && ~strcmp('85',subj) && ~strcmp('87',subj)
                anat_selec_index = zeros(1,length(chan_label));
                cue_resp_index   = zeros(1,length(chan_label));

        for iChan_SG = 1:length(chan_label)
            for iChan_mni = 1:length(chan_label_mni)
                if strcmp(chan_label_mni{iChan_mni}, chan_label{iChan_SG})
                    desired_chan_idx (iChan_SG) = iChan_mni;
                    if ismember (iChan_SG, anatomical_selec)
                    anat_selec_index (iChan_SG) = 1;
                    end
                    
                    if ismember (iChan_SG, cue_resp_elec)
                        cue_resp_index (iChan_SG) = 1;
                    end
                    
                end
            end
        end
        
        if any(desired_chan_idx==0)
            desired_chan_idx = desired_chan_idx(desired_chan_idx~=0);
            chan_label = chan_label(desired_chan_idx);
            chan_label_reduc_for_mni = chan_label;
            save('chan_label_reduc_for_mni','chan_label_reduc_for_mni')
        end
     elec_pos_reordered= elec_pos(desired_chan_idx,:);
    for iElec = 1:length(elec_pos_reordered)
        i = i+1;
        subj_chan_name{i} =[subj chan_label{iElec}];
        X(i) =elec_pos_reordered(iElec,1);
        Y(i) =elec_pos_reordered(iElec,2);
        Z(i) =elec_pos_reordered(iElec,3);
        anat_selec_index_group(i) = anat_selec_index (iElec);
       cue_resp_index_group(i) = cue_resp_index (iElec);
    end
    
    
    elseif strcmp('84',subj) || strcmp('85',subj) || strcmp('87',subj)
        anat_selec_index = zeros(1,length(chan_label_mni));
        cue_resp_index = zeros(1,length(chan_label_mni));
        
        for iChan_mni = 1:length(chan_label_mni)
            for iChan_SG = 1:length(chan_label)
                if strcmp(chan_label_mni{iChan_mni}, chan_label{iChan_SG})
                    desired_chan_idx (iChan_mni) = iChan_SG;
                    
                    if ismember (iChan_SG, anatomical_selec)
                    anat_selec_index (iChan_mni) = 1;
                    end
                    
                    if ismember (iChan_SG, cue_resp_elec)
                        cue_resp_index (iChan_mni) = 1;
                    end
                end
            end
        end
        
        if length(desired_chan_idx)== length(chan_label_mni) % if all mni chan label are found in my data, then save them as is
            for iElec = 1:length(chan_label_mni)
                i = i+1;
                subj_chan_name{i} =[subj chan_label_mni{iElec}];
                X(i) =elec_pos(iElec,1);
                Y(i) =elec_pos(iElec,2);
                Z(i) =elec_pos(iElec,3);
               anat_selec_index_group(i) = anat_selec_index (iElec);
               cue_resp_index_group(i) = cue_resp_index (iElec);
            end    
        end        
    end
    
        clearvars -except iSubj i subj_list Z X Y subj_chan_name anat_selec_index_group cue_resp_index_group

end
% append a number to indicate order for subject name
for iName = 1:i
subj_chan_name{iName} = [num2str(iName) '_' subj_chan_name{iName} ];
end
Dat_Table = table(subj_chan_name(1:i,:),X(1:i,:),Y(1:i,:),Z(1:i,:),anat_selec_index_group(1:i), cue_resp_index_group(1:i));


%% make black region label
region_label = zeros(size(Dat_Table,1),3);% all coverage
region_label(logical(anat_selec_index_group(1:i)==1&cue_resp_index_group(1:i)==0),2) = 255;% anatomy 
region_label(logical(cue_resp_index_group(1:i)),1) = 255;% cue resp
region_label(logical(cue_resp_index_group(1:i)),3) = 255;

