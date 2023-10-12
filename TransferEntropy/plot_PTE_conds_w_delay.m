% plot condition traces as a function of delay

clear all
sig_elecs   = '';%'OFC_FRO_TEMP_CING_INS_EC'
fn_ext      = 'cue_resp'; %
freq_range  = 'deltatheta';
fpass  = [4;5]';
time1 = 0;
time2 = 1;
lock  = 'onset';
bidirec = 0.5;

%lock = 'response'
phase = 'encoding'
phase = 'retrieval';

reg1_list = {'OFC' 'FRO' 'TEMP'  }; %' 'CING' 'INS'
reg2_list = {'HC'};

addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
cd (['/mnt/yassamri/iEEG/sandra/PTE_results/' fn_ext '/' freq_range '/' phase '/' lock '/'])
cd([num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1) '_' num2str(time2) 'sec'])


mtx_sz = 40;
% cond 2
reg1_reg2_conda = cell(1,mtx_sz);
% cond 3
reg1_reg2_condb = cell(1,mtx_sz);

subj_list = {'39' '44' '57'  '63' '66' '84' '85' '87' }

i = 0;
for iReg1 = 1:length(reg1_list)
    
    % get NC site
    reg1_name = reg1_list{iReg1};
    reg2_name = reg2_list{1};
    
    % get all data from each subj
    for sub_counter = 1:length(subj_list)
        fn1 = [reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2' sig_elecs '_FOTofLag.mat' ];
        if isfile(fn1)
            i = i+1;
            
            % load cond2
            load(fn1)
            reg1_reg2_conda{i} = PTE_ch1_to_ch2_norm;
            clear PTE_ch1_to_ch2_norm
            % load cond3
            if strcmp('encoding',phase)
                fn2 = [reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' sig_elecs '_FOTofLag.mat' ];
            else
                fn2 = [reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3' sig_elecs '_FOTofLag.mat' ];
            end
            load(fn2)
            
            reg1_reg2_condb{i} = PTE_ch1_to_ch2_norm;
            clear PTE_ch1_to_ch2_norm
        end
    end
    
end

reg1_reg2_conda = cat(1,reg1_reg2_conda{:});
reg1_reg2_condb = cat(1,reg1_reg2_condb{:});
%%
figure;
hold on

x_stps =(0.05:.025:.5)*1000; %to get msec          

stdshade(reg1_reg2_condb-bidirec,.1,'m',x_stps, [],[])
h1 = plot(x_stps, nanmean(reg1_reg2_condb,1)-bidirec, 'm', 'LineWidth', 3)

stdshade(reg1_reg2_conda-bidirec,.1,'b',x_stps, [],[])
h2 = plot(x_stps, nanmean(reg1_reg2_conda,1)-bidirec, 'b', 'LineWidth', 3)

xlim([x_stps(1) x_stps(end)])
%ylim([-.005 .005])
line([x_stps(1) x_stps(end)],[0 0], 'Color', 'k', 'LineWidth',3)
legend([h1, h2],{'lure +' 'lure -' })
xlabel(['delay (msec)']);
ylabel(['PTE'])
title([phase ' ' lock ' NC vs. '  reg2_list{1}]) 
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
