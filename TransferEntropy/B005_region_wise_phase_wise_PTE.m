clear all;close all;clc
figure
hold on
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')

% define params
group_plots = 'yes'
freq_range  = 'deltatheta';
fpass       = [4 5];
times1      = {'0' '0' '0.5'};%'0' '0' '0.5'
times2      = {'2' '1' '1.5'};%'2' '1' '1.5'
bidirec     = 0.5;
pattern      = 1

phase_list  = {'encoding' 'retrieval' 'retrieval'};% 'retrieval'
lock_list   = {'onset'    'onset'     'response'    }; %'response'

% define regions
reg2_list = {'HC'};

if pattern ==1
    reg1_list = {'OFC' 'TEMP' 'FRO' };  
    colorList = {'r-o', 'c-o', 'g-o'};
    color     = { 'r','c', 'g'};
elseif pattern == 2
    reg1_list = {'CING' 'INS' 'EC' };
    colorList = {'m-o', 'k-o', 'b-o',};
    color     = { 'm','k', 'b'};
end


% init
phase_data = cell(1,length(phase_list));
mtx_sz     = 40;

for iReg1 = 1:length(reg1_list)

    % get NC site
    reg1_name = reg1_list{iReg1};
    reg2_name = reg2_list{1};
    
    if strcmp('yes',group_plots)
        subj_list = {'39' '44' '57'  '63' '66' '84' '85' '87'};
    else
        subj_list{1} = subj;
    end
    
    for iPhase = 1:length(phase_list)
        time1 = times1{iPhase};
        time2 = times2{iPhase};
        phase = phase_list{iPhase};
        lock  = lock_list{iPhase};
        
        cd (['/mnt/yassamri/iEEG/sandra/PTE_results/cue_resp/' freq_range '/' phase '/' lock '/' num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1) '_' num2str(time2) 'sec'])

        % get all data from each subj
        i = 0;
        reg1_reg2_condb = cell(1,mtx_sz);
        for sub_counter = 1:length(subj_list)
            i = i+1;
            if strcmp('retrieval',phase)
               if isfile([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3.mat' ])
                % load lure +
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3' ])
                reg1_reg2_condb{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
               end
                
            else
                 if isfile([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1.mat' ])
                % load lure +
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' ])
                reg1_reg2_condb{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                 end
            end
        end
        
        % lure+
        reg1_reg2_condb    = cat(1,reg1_reg2_condb{:});
        
        phase_data{iPhase} = reg1_reg2_condb-bidirec;

    end
    
    conda = phase_data{1};
    condb = phase_data{2};
    condc = phase_data{3};
    
    % save source data
   phaseDataPerRegion{iReg1} = phase_data;

    %% plot only first two conditions

plot([nanmean(conda) nanmean(condb) ],colorList{iReg1}, 'MarkerFaceColor', color{iReg1}, 'LineWidth', 1)
errorbar(1:2,[nanmean(conda) nanmean(condb) ],[nanstd(conda)/sqrt(length(conda)) nanstd(condb)/sqrt(length(condb))], colorList{iReg1},'LineWidth',1)
set(gca, 'XTick', 1:2, 'XTickLabel', {'encode lure+' 'ret onst lure+' },'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
    % Uncomment this to plot all three phases
    plot([nanmean(conda) nanmean(condb) nanmean(condc)],colorList{iReg1}, 'MarkerFaceColor', color{iReg1}, 'LineWidth', 1)
    errorbar(1:3,[nanmean(conda) nanmean(condb) nanmean(condc)],[nanstd(conda)/sqrt(length(conda)) nanstd(condb)/sqrt(length(condb)) nanstd(condc)/sqrt(length(condc))], colorList{iReg1},'LineWidth',1)
    set(gca, 'XTick', 1:3, 'XTickLabel', {'encode lure+' 'ret onst lure+' 'ret resp lure+'},'XTickLabelRotation',45)
    set(gca, 'FontSize', 16, 'FontWeight', 'bold')
    
    %% stats for phase plots
    % change second input to cond a vs. cond c to compare retreival onset
    % vs. resonse with encoding, respectively
    [p,ef] = permutation_paired(condc, conda, 1000); 
    p_val(iReg1) =  p;
    effect(iReg1) =  ef;

end
% after getting all reg
[p_fdr, p_masked] = fdr( p_val, 0.05);



