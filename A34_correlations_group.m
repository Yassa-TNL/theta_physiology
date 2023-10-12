x_reg = 'OFC'
y_reg = 'CA3'
exp_type = 'tuning_correct'
lock     = 'onset'

if strcmp('OFC',x_reg)
    subj_list  = {'39' '57' '44' '63' '66' '84'}
elseif strcmp('FRO',x_reg)
    subj_list  = {'39' '57' '44' '63' '66' '83' '84'}    
elseif strcmp('TEMP',x_reg)
    subj_list  = {'39' '57' '44' '63' '66' '83' '84'}
elseif strcmp('CING',x_reg)
    subj_list  = {'39' '57' '44' '63' '66' }    
elseif strcmp('MTL',x_reg)
    subj_list  = {'39' '57' '44' '63' '66'  '84'}    
elseif strcmp('CA3',x_reg)
    subj_list  = {'39' '57' '44'  '66' }
elseif strcmp('ins',x_reg)
    subj_list  = {'39' '57'}
end

if strcmp('response', lock)
    cd(['/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_response'])
    time_range = [1 1301];
elseif strcmp('onset', lock)
    cd(['/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_onset_' exp_type])
    if strcmp('encoding', exp_type)
        time_range = [1 2101];
        pre_stim   =.5;
        post_stim  = 2;
    else
        time_range = [301 1301];
        pre_stim   = 0;
        post_stim  = 1;
    end
end
total_subj = length(subj_list);
cd('/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_onset_tuning_correct')

%% x axis data
cond2a      = cell(1,total_subj);
cond3a      = cell(1,total_subj);
for subj = 1:total_subj
    load(['subj' subj_list{subj} 'spectrograms.mat'])
    if strcmp('OFC',x_reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = OFC_cond2;
            cond3a{subj} = OFC_cond1;
        else
            cond2a{subj} = OFC_cond2;
            cond3a{subj} = OFC_cond3;
        end
        
    elseif strcmp('FRO',x_reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = fro_cond2;
            cond3a{subj} = fro_cond1;
        else
            cond2a{subj} = fro_cond2;
            cond3a{subj} = fro_cond3;
        end
        
    elseif strcmp('TEMP',x_reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = temp_cond2;
            cond3a{subj} = temp_cond1;        
        else
            cond2a{subj} = temp_cond2;
            cond3a{subj} = temp_cond3;
        end
    elseif strcmp('CING',x_reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = cing_cond2;
            cond3a{subj} = cing_cond1;            
        else
            cond2a{subj} = cing_cond2;
            cond3a{subj} = cing_cond3;
        end
        
         
    elseif strcmp('ins',x_reg)
        if strcmp('encoding', exp_type)
            cond2a{subj} = ins_cosnd2;
            cond3a{subj} = ins_cond1;            
        else
            cond2a{subj} = ins_cond2;
            cond3a{subj} = ins_cond3;
        end
        
    end
end


cond2_x = cat(3,cond2a{:});
cond3_x = cat(3,cond3a{:});

%% y axis
cond2a      = cell(1,total_subj);
cond3a      = cell(1,total_subj);
for subj = 1:total_subj
    load(['subj' subj_list{subj} 'spectrograms.mat'])
    if strcmp('MTL',y_reg)
        
        cond2a{subj} = MTL_cond2;
        cond3a{subj} = MTL_cond3;

    elseif strcmp('CA3',y_reg)
        cond2a{subj} = CA3_cond2;
        cond3a{subj} = CA3_cond3;
    end
end

cond2_y = cat(3,cond2a{:});
cond3_y = cat(3,cond3a{:});
%% 

 cd(['/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_' lock '_' exp_type '/cluster_matrices'])
 load('CA3_onset_group_cluster_delta.mat')
 y_axis_clust = zmapthresh;