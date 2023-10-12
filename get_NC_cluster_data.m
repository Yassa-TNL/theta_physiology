function [cond1a,cond2a,cond3a,cond4a] = get_NC_cluster_data(subj_list, ref, baseline, exp_type, lock)
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock ])

regList = { 'OFC' 'FRO' 'TEMP' 'CING' 'INS' 'EC'};
cntr = 0;
for iReg = 1:length(regList)
    
    for subj = 1:length(subj_list)
        load(['subj' subj_list{subj} 'spectrograms' baseline '.mat'])
        
        if strcmp('encoding',exp_type)
            if  strcmp('OFC',regList{iReg})
                data2 = OFC_cond2;
                data3 = OFC_cond1;
                
            elseif strcmp('FRO',regList{iReg})
                data2 = fro_cond2;
                data3 = fro_cond1;
                
            elseif strcmp('TEMP',regList{iReg})
                data2 = temp_cond2;
                data3 = temp_cond1;
                
            elseif strcmp('CING',regList{iReg})
                data2 = cing_cond2;
                data3 = cing_cond1;
                
            elseif strcmp('INS',regList{iReg})
                data2 = ins_cond2;
                data3 = ins_cond1;
                
            elseif strcmp('EC',regList{iReg})
                data2 = EC_cond2;
                data3 = EC_cond1;
                
            end
            
            if ~isempty(data2)
                cntr = cntr+1;
                cond2a{cntr} = nanmean(data2,3);
                cond3a{cntr} = nanmean(data3,3);
                cond1a=[];cond4a=[];
            end
        else strcmp('tuning_correct', exp_type)
            
            
            
             if  strcmp('OFC',regList{iReg})
                data1 = OFC_cond1;
                data2 = OFC_cond2;
                data3 = OFC_cond3;
                data4 = OFC_cond4;
            elseif strcmp('FRO',regList{iReg})
                data1 = fro_cond1;
                data2 = fro_cond2;
                data3 = fro_cond3;
                data4 = fro_cond4;                
            elseif strcmp('TEMP',regList{iReg})
                data1 = temp_cond1;
                data2 = temp_cond2;
                data3 = temp_cond3;
                data4 = temp_cond4; 
                
            elseif strcmp('CING',regList{iReg})
                data1 = cing_cond1;
                data2 = cing_cond2;
                data3 = cing_cond3;
                data4 = cing_cond4; 
                
            elseif strcmp('INS',regList{iReg})
                data1 = ins_cond1;
                data2 = ins_cond2;
                data3 = ins_cond3;
                data4 = ins_cond4; 
                
            elseif strcmp('EC',regList{iReg})
                data1 = EC_cond1;
                data2 = EC_cond2;
                data3 = EC_cond3;
                data4 = EC_cond4; 
                
            end
            
            if ~isempty(data2)
                cntr = cntr+1;
                cond1a{cntr} = nanmean(data1,3);
                cond2a{cntr} = nanmean(data2,3);
                cond3a{cntr} = nanmean(data3,3);
                cond4a{cntr} = nanmean(data4,3);
            end
            
            
            
            
            
        end
    end
end


end

