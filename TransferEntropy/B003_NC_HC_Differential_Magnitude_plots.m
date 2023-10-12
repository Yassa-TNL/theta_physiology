
clear all;%close all;clc
fn_ext           = 'cue_resp' %

group_plots      = 'yes'
subj             = ''
cond             = 'lure+'
freq_range       = 'deltatheta'
fpass            = [4 5]
time1            = '0'
time2            = '1'
plot_type        = 'magnitude' % main phase % magnitude
lock             = 'onset'
%lock            = 'response'
phase            = 'encoding'
phase           = 'retrieval'
bidirec          = 0.5

if strcmp('39', subj)  || strcmp('57', subj) || strcmp('84', subj) || strcmp('85', subj) ||  strcmp('yes',group_plots)
    reg1_list        = {'OFC' 'FRO' 'TEMP' 'CING' 'INS' 'EC'};
elseif (strcmp('44', subj)  || strcmp('63', subj)) &&  strcmp('',group_plots)
    reg1_list        = {'OFC' 'FRO' 'TEMP' 'CING'  'EC'};
elseif strcmp('66', subj) &&  strcmp('',group_plots)
    reg1_list        = {'OFC' 'FRO' 'TEMP' 'CING' };
elseif strcmp('87', subj) &&  strcmp('',group_plots)
    reg1_list        = {'OFC' 'FRO'  'TEMP' 'CING' 'INS'};
end
    reg1_list        = {'OFC' 'FRO'  'TEMP' 'CING' 'INS' 'EC'};

reg2_list = {'HC'};

addpath('/tmp/yassamri/iEEG/sandra/analysis_pipeline_final')
cd (['/tmp/yassamri/iEEG/sandra/PTE_results/' fn_ext '/' freq_range '/' phase '/' lock '/' num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1) '_' num2str(time2) 'sec'])
reg_spec_values_cond3 = cell(1,length(reg1_list));
reg_spec_values_cond2 = cell(1,length(reg1_list));
for iReg1 = 1:length(reg1_list)
    
    % get NC site
    reg1_name = reg1_list{iReg1};
    reg2_name = reg2_list{1};
    if strcmp('yes',group_plots)
        % get subj list
            subj_list = {'39' '44' '57'  '63' '66' '84' '85' '87'}
    else
        subj_list{1} = subj;
    end
    
    % prepare this for all conds
    i = 0
    % cond 2: lure -
    reg1_reg2_conda = cell(1,length(subj_list));
    % cond 3: lure+
    reg1_reg2_condb = cell(1,length(subj_list));
    % cond 1: repeat
    reg1_reg2_condc = cell(1,length(subj_list));
    % cond 4: new
    reg1_reg2_condd = cell(1,length(subj_list));
    
    % get all data from each subj
    for sub_counter = 1:length(subj_list)
        i = i+1;
        if strcmp('retrieval',phase)
            % load lure -
            if isfile([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2.mat' ])
            load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2' ])
            
                
            reg1_reg2_conda{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
            else 
                
                continue
            end
            % load lure +
            load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3' ])
            reg1_reg2_condb{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
            
            % load repeat
            load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' ])
            reg1_reg2_condc{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
            
            % load new
            load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond4' ])
            reg1_reg2_condd{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
        else
            if isfile([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2.mat' ])
            % load lure -
            load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2' ])
            reg1_reg2_conda{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
            
            % load lure +
            load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' ])
            reg1_reg2_condb{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
            else
                continue
            end
            end
    end
    reg_spec_values_cond3{iReg1} = cat(1,reg1_reg2_condb{:});
    reg_spec_values_cond2{iReg1} = cat(1,reg1_reg2_conda{:});
    reg_spec_values_cond1{iReg1} = cat(1,reg1_reg2_condc{:});
    reg_spec_values_cond4{iReg1} = cat(1,reg1_reg2_condd{:});
end


% plot
mn = -.06;
mx = .06;

if   strcmp('lure+',cond)
    data_mtx = reg_spec_values_cond3;
elseif strcmp('lure-',cond)
    data_mtx = reg_spec_values_cond2;
elseif strcmp('new',cond)
    data_mtx = reg_spec_values_cond4;
elseif strcmp('repeat',cond)
    data_mtx = reg_spec_values_cond1;
end
    
% plot
if strcmp('encoding',phase)
        figure
    hold on
     if strcmp('yes',group_plots) || strcmp('39', subj)  || strcmp('57', subj) || strcmp('84', subj) || strcmp('85', subj)
         

    men = [nanmean(data_mtx{2})-bidirec nanmean(data_mtx{1})-bidirec nanmean(data_mtx{4})-bidirec nanmean(data_mtx{5})-bidirec...
        nanmean(data_mtx{6})-bidirec nanmean(data_mtx{3})-bidirec]
    se  = [nanstd(data_mtx{2})/sqrt(length(data_mtx{2})) nanstd(data_mtx{1})/sqrt(length(data_mtx{1}))  nanstd(data_mtx{4})/sqrt(length(data_mtx{4})) ...
       nanstd(data_mtx{5})/sqrt(length(data_mtx{5}))  nanstd(data_mtx{6})/sqrt(length(data_mtx{6}))  nanstd(data_mtx{3})/sqrt(length(data_mtx{3}))]
    bar(men)
    errorbar(1:6,men, se, 'rx')
    set(gca, 'XTick', 1:6, 'XTickLabel', {'FRO ' 'OFC '  ...
         'CING' 'INS' 'EC' 'TEMP'},'XTickLabelRotation',45)
     elseif strcmp('',group_plots)&& (strcmp('44', subj)  || strcmp('63', subj))
    men = [nanmean(data_mtx{2})-bidirec nanmean(data_mtx{1})-bidirec nanmean(data_mtx{4})-bidirec ...
        nanmean(data_mtx{5})-bidirec nanmean(data_mtx{3})-bidirec]
    se  = [nanstd(data_mtx{2})/sqrt(length(data_mtx{2})) nanstd(data_mtx{1})/sqrt(length(data_mtx{1}))  nanstd(data_mtx{4})/sqrt(length(data_mtx{4})) ...
       nanstd(data_mtx{5})/sqrt(length(data_mtx{5}))  nanstd(data_mtx{3})/sqrt(length(data_mtx{3}))]
    bar(men)
    errorbar(1:5,men, se, 'rx')
    set(gca, 'XTick', 1:5, 'XTickLabel', {'FRO ' 'OFC '  ...
         'CING'  'EC' 'TEMP'},'XTickLabelRotation',45) 
     elseif strcmp('',group_plots)&& ( strcmp('66', subj))
      men = [nanmean(data_mtx{2})-bidirec nanmean(data_mtx{1})-bidirec nanmean(data_mtx{4})-bidirec ...
         nanmean(data_mtx{3})-bidirec]
    se  = [nanstd(data_mtx{2})/sqrt(length(data_mtx{2})) nanstd(data_mtx{1})/sqrt(length(data_mtx{1}))  nanstd(data_mtx{4})/sqrt(length(data_mtx{4})) ...
        nanstd(data_mtx{3})/sqrt(length(data_mtx{3}))]
    bar(men)
    errorbar(1:4,men, se, 'rx')
    set(gca, 'XTick', 1:4, 'XTickLabel', {'FRO ' 'OFC '  ...
         'CING'   'TEMP'},'XTickLabelRotation',45)   
     elseif strcmp('',group_plots)&& ( strcmp('87', subj))
     men = [nanmean(data_mtx{2})-bidirec nanmean(data_mtx{1})-bidirec nanmean(data_mtx{4})-bidirec nanmean(data_mtx{5})-bidirec ...
         nanmean(data_mtx{3})-bidirec]
    se  = [nanstd(data_mtx{2})/sqrt(length(data_mtx{2})) nanstd(data_mtx{1})/sqrt(length(data_mtx{1}))  nanstd(data_mtx{4})/sqrt(length(data_mtx{4})) ...
        nanstd(data_mtx{5})/sqrt(length(data_mtx{5})) nanstd(data_mtx{3})/sqrt(length(data_mtx{3}))]
    bar(men)
    errorbar(1:5,men, se, 'rx')
    set(gca, 'XTick', 1:5, 'XTickLabel', {'FRO ' 'OFC '  ...
         'CING' 'INS'  'TEMP'},'XTickLabelRotation',45)          
   
     end
     

elseif strcmp('retrieval',phase)&& strcmp('onset',lock)
    figure
    hold on
     if strcmp('yes',group_plots) || strcmp('39', subj)  || strcmp('57', subj) || strcmp('84', subj) || strcmp('85', subj)
         

    men = [nanmean(data_mtx{3})-bidirec nanmean(data_mtx{5})-bidirec nanmean(data_mtx{4})-bidirec...
        nanmean(data_mtx{1})-bidirec  nanmean(data_mtx{6})-bidirec nanmean(data_mtx{2})-bidirec  ]
    se  = [nanstd(data_mtx{3})/sqrt(length(data_mtx{3})) nanstd(data_mtx{5})/sqrt(length(data_mtx{5})) ...
        nanstd(data_mtx{4})/sqrt(length(data_mtx{4})) nanstd(data_mtx{1})/sqrt(length(data_mtx{1}))    ...
         nanstd(data_mtx{6})/sqrt(length(data_mtx{6})) nanstd(data_mtx{2})/sqrt(length(data_mtx{2})) ]
    bar(men)
    errorbar(1:6,men, se, 'rx')
    set(gca, 'XTick', 1:6, 'XTickLabel', {'TEMP ' 'INS ' 'CING' 'OFC' 'EC' 'FRO'},'XTickLabelRotation',45)
     end
elseif strcmp('retrieval',phase) && strcmp('response',lock)
    if strcmp('gamma',freq_range)
    figure
    hold on
    men = [nanmean(data_mtx{2}) nanmean(data_mtx{5})  nanmean(data_mtx{6}) ...
         nanmean(data_mtx{3}) nanmean(data_mtx{1}) nanmean(data_mtx{4})]-bidirec 
        
    se  = [nanstd(data_mtx{2})/sqrt(length(data_mtx{2})) nanstd(data_mtx{5})/sqrt(length(data_mtx{5})) nanstd(data_mtx{6})/sqrt(length(data_mtx{6}))...
         nanstd(data_mtx{3})/sqrt(length(data_mtx{3})) nanstd(data_mtx{1})/sqrt(length(data_mtx{1}))  nanstd(data_mtx{4})/sqrt(length(data_mtx{4}))]
    bar(men)
    errorbar(1:6,men, se, 'rx')
    set(gca, 'XTick', 1:6, 'XTickLabel', {'FRO' 'INS' 'EC' 'TEMP' 'OFC '  'CING'  },'XTickLabelRotation',45)
    ylim([-0.002 0.003])
    else
        figure
        hold on
        men = [nanmean(data_mtx{2}) nanmean(data_mtx{5})  nanmean(data_mtx{6}) ...
            nanmean(data_mtx{3}) nanmean(data_mtx{1}) nanmean(data_mtx{4})]-bidirec
        
        se  = [nanstd(data_mtx{2})/sqrt(length(data_mtx{2})) nanstd(data_mtx{5})/sqrt(length(data_mtx{5})) nanstd(data_mtx{6})/sqrt(length(data_mtx{6}))...
            nanstd(data_mtx{3})/sqrt(length(data_mtx{3})) nanstd(data_mtx{1})/sqrt(length(data_mtx{1}))  nanstd(data_mtx{4})/sqrt(length(data_mtx{4}))]
        bar(men)
        errorbar(1:6,men, se, 'rx')
        set(gca, 'XTick', 1:6, 'XTickLabel', {'FRO' 'INS' 'EC' 'TEMP' 'OFC '  'CING'  },'XTickLabelRotation',45)
          %  ylim([-0.009 0.008])

    end
end

set(gca, 'FontSize', 16, 'FontWeight', 'bold')
%y=ylim
%yticks([y(1) 0 y(2)])
title(cond)
cd('/tmp/yassamri/iEEG/sandra/GroupFigures')