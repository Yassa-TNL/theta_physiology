function [cond1a,cond2a,cond3a,cond4a] = get_power_subj_elecs(subj_list,exp_type,reg, ref, baseline,lock)
% Function output: condition specific spectrogram power from every elec pooled
% from all subjs for a given region
cd('/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/LM_reref/tuning_correct_onset')

% load OFC FRO TEMP chan specificity
load('subj_NC_specificty')
load('subj_ALL_NC_specificty')
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock ])

cond1a      = cell(1,length(subj_list));
cond2a      = cell(1,length(subj_list));
cond3a      = cell(1,length(subj_list));
cond4a      = cell(1,length(subj_list));
for subj = 1:length(subj_list)
    load(['subj' subj_list{subj} 'spectrograms' baseline '.mat'])
    if strcmp('OFC',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = OFC_cond2;
            cond3a{subj} = OFC_cond1;
        else
            cond1a{subj} = OFC_cond1;
            cond2a{subj} = OFC_cond2;
            cond3a{subj} = OFC_cond3;
            cond4a{subj} = OFC_cond4;
        end
        
    elseif strcmp('FRO',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = fro_cond2;
            cond3a{subj} = fro_cond1;
        else
            cond1a{subj} = fro_cond1;
            cond2a{subj} = fro_cond2;
            cond3a{subj} = fro_cond3;
            cond4a{subj} = fro_cond4;
        end
    elseif strcmp('EC',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = EC_cond2;
            cond3a{subj} = EC_cond1;
        else
            cond1a{subj} = EC_cond1;
            cond2a{subj} = EC_cond2;
            cond3a{subj} = EC_cond3;
            cond4a{subj} = EC_cond4;
        end
    elseif strcmp('TEMP',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = temp_cond2;
            cond3a{subj} = temp_cond1;
        else
            cond1a{subj} = temp_cond1;
            cond2a{subj} = temp_cond2;
            cond3a{subj} = temp_cond3;
            cond4a{subj} = temp_cond4;
        end
    elseif strcmp('CING',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = cing_cond2;
            cond3a{subj} = cing_cond1;
        else
            cond1a{subj} = cing_cond1;
            cond2a{subj} = cing_cond2;
            cond3a{subj} = cing_cond3;
            cond4a{subj} = cing_cond4;
        end
        
    elseif strcmp('MTL',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = MTL_cond2;
            cond3a{subj} = MTL_cond1;
        else
            cond1a{subj} = MTL_cond1;
            cond2a{subj} = MTL_cond2;
            cond3a{subj} = MTL_cond3;
            cond4a{subj} = MTL_cond4;
        end
    elseif strcmp('CA1',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = CA1_cond2;
            cond3a{subj} = CA1_cond1;
        else
            cond1a{subj} = CA1_cond1;
            cond2a{subj} = CA1_cond2;
            cond3a{subj} = CA1_cond3;
            cond4a{subj} = CA1_cond4;
        end
        
    elseif strcmp('CA3',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = CA3_cond2;
            cond3a{subj} = CA3_cond1;
        else
            cond1a{subj} = CA3_cond1;
            cond2a{subj} = CA3_cond2;
            cond3a{subj} = CA3_cond3;
            cond4a{subj} = CA3_cond4;
        end
    elseif strcmp('HC',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = HC_cond2;
            cond3a{subj} = HC_cond1;
        else
            cond1a{subj} = HC_cond1;
            cond2a{subj} = HC_cond2;
            cond3a{subj} = HC_cond3;
            cond4a{subj} = HC_cond4;
        end
        % compare anat diff in ephy resp for lure+ cond
    elseif strcmp('CA3_CA1_lure+',reg) 
            cond3a{subj} = CA3_cond3;
            cond2a{subj} = CA1_cond3; 
       elseif strcmp('CA3_CA1_lure-',reg)
            cond3a{subj} = CA3_cond2;
            cond2a{subj} = CA1_cond2;      
    elseif strcmp('ins',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = ins_cond2;
            cond3a{subj} = ins_cond1;
        else
            cond1a{subj} = ins_cond1;
            cond2a{subj} = ins_cond2;
            cond3a{subj} = ins_cond3;
            cond4a{subj} = ins_cond4;
        end
    elseif strcmp('NC',reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = NC_cond2;
            cond3a{subj} = NC_cond1;
        else
            cond1a{subj} = NC_cond1;
            cond2a{subj} = NC_cond2;
            cond3a{subj} = NC_cond3;
            cond4a{subj} = NC_cond4;
        end
    elseif strcmp('OFC_FRO_TEMP',reg)
        % get OFC_FRO_TEMP
        if strcmp('encoding', exp_type)
            cond2a{subj} = cat(3,OFC_cond2,fro_cond2,temp_cond2)
            cond3a{subj} = cat(3,OFC_cond1,fro_cond1,temp_cond1);
        else
            cond1a{subj} = cat(3,OFC_cond1,fro_cond1,temp_cond1);
            cond2a{subj} = cat(3,OFC_cond2,fro_cond2,temp_cond2);
            cond3a{subj} = cat(3,OFC_cond3,fro_cond3,temp_cond3);
            cond4a{subj} = cat(3,OFC_cond4,fro_cond4,temp_cond4);
        end
    elseif strcmp('OFC_FRO_TEMP_significant',reg)
        sig_chan_log_vec = logical(subj_NC_specificty{subj}(:,2)<0.05);
        if strcmp('encoding', exp_type)
            temp2 = cat(3,OFC_cond2,fro_cond2,temp_cond2);
            temp3 = cat(3,OFC_cond1,fro_cond1,temp_cond1);
            
            cond2a{subj} = temp2(:,:,sig_chan_log_vec);
            cond3a{subj} = temp3(:,:,sig_chan_log_vec);
            
        else
            % get all chans
            temp1 = cat(3,OFC_cond1,fro_cond1,temp_cond1);
            temp2 = cat(3,OFC_cond2,fro_cond2,temp_cond2);
            temp3 = cat(3,OFC_cond3,fro_cond3,temp_cond3);
            temp4 = cat(3,OFC_cond4,fro_cond4,temp_cond4);
            
            % index significant one
            cond1a{subj} = temp1(:,:,sig_chan_log_vec);
            cond2a{subj} = temp2(:,:,sig_chan_log_vec);
            cond3a{subj} = temp3(:,:,sig_chan_log_vec);
            cond4a{subj} = temp4(:,:,sig_chan_log_vec);
        end
        
    elseif strcmp('OFC_FRO_TEMP_CING_INS_EC_significant',reg)
        sig_chan_log_vec = logical(subj_ALL_NC_specificty{subj}(:,2)<0.05);
        if strcmp('encoding', exp_type)
            temp2 = cat(3,OFC_cond2,fro_cond2,temp_cond2,cing_cond2, ins_cond2, EC_cond2);
            temp3 = cat(3,OFC_cond1,fro_cond1,temp_cond1,cing_cond1, ins_cond1, EC_cond1);
            
            cond2a{subj} = temp2(:,:,sig_chan_log_vec);
            cond3a{subj} = temp3(:,:,sig_chan_log_vec);
            
        else
            % get all chans
            temp1 = cat(3,OFC_cond1,fro_cond1,temp_cond1,cing_cond1, ins_cond1, EC_cond1);
            temp2 = cat(3,OFC_cond2,fro_cond2,temp_cond2,cing_cond2, ins_cond2, EC_cond2);
            temp3 = cat(3,OFC_cond3,fro_cond3,temp_cond3,cing_cond3, ins_cond3, EC_cond3);
            temp4 = cat(3,OFC_cond4,fro_cond4,temp_cond4,cing_cond4, ins_cond4, EC_cond4);
            
            % index significant one
            cond1a{subj} = temp1(:,:,sig_chan_log_vec);
            cond2a{subj} = temp2(:,:,sig_chan_log_vec);
            cond3a{subj} = temp3(:,:,sig_chan_log_vec);
            cond4a{subj} = temp4(:,:,sig_chan_log_vec);
        end
        
    end
    
end
end


