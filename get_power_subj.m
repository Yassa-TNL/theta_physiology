function [cond1,cond2,cond3,cond4] = get_power_subj(subj_list,exp_type,reg, ref, lock)
% Function output: mean spectrogram power for a given region, using a
% single matrix per subj for all conds
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock ])

for iSubj = 1:length(subj_list)
    % load subj data
    load(['subj' subj_list{iSubj} 'spectrograms.mat'])
    
    if strcmp('OFC',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(OFC_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(OFC_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(OFC_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(OFC_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(OFC_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(OFC_cond4,3); % new
        end
        
    elseif strcmp('FRO',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(fro_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(fro_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(fro_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(fro_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(fro_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(fro_cond4,3); % new
        end
        
    elseif strcmp('TEMP',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(temp_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(temp_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(temp_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(temp_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(temp_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(temp_cond4,3); % new
        end
    elseif strcmp('CING',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(cing_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(cing_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(cing_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(cing_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(cing_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(cing_cond4,3); % new
        end
        
    elseif strcmp('ins',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(ins_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(ins_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(ins_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(ins_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(ins_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(ins_cond4,3); % new
        end
    elseif strcmp('EC',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(EC_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(EC_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(EC_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(EC_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(EC_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(EC_cond4,3); % new
            cond3 = [];
            cond4 = [];
        end
    elseif strcmp('HC',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(HC_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(HC_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(HC_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(HC_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(HC_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(HC_cond4,3); % new
        end
    elseif strcmp('CA1',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(CA1_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(CA1_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(CA1_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(CA1_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(CA1_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(CA1_cond4,3); % new
        end
    elseif strcmp('CA3',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(CA3_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(CA3_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(CA3_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(CA3_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(CA3_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(CA3_cond4,3); % new
        end
        
    elseif strcmp('NC',reg)
        if strcmp('encoding', exp_type)
            cond1(:,:,iSubj) = nanmean(NC_cond1,3); % --> lure +
            cond2(:,:,iSubj) = nanmean(NC_cond2,3); % --> lure -
            cond3 = [];
            cond4 = [];
        else
            cond1(:,:,iSubj) = nanmean(NC_cond1,3); % rep
            cond2(:,:,iSubj) = nanmean(NC_cond2,3); % lure -
            cond3(:,:,iSubj) = nanmean(NC_cond3,3); % lure +
            cond4(:,:,iSubj) = nanmean(NC_cond4,3); % new
        end
        
end
end

