clear all;close all;clc
group_plots      = 'yes'
iSubj            = 1% change this for indiv subj line color
subj             = '39'
phase            = 'retrieval'
freq_range       = 'deltatheta'
fn_ext           = 'cue_resp' %
plot_4_conds     = ''
locks            = {'onset'  }
times1           = {'0' '0.5'}
times2           = {'1' '1.5'}
fpass_list       = [4; 5]'


reg1_list   = {'OFC' 'FRO' 'TEMP'};
reg2_list   = {'HC'};
if strcmp('yes',group_plots)
    subj_list = {'39' '44' '57'  '63' '66' '84' '85' '87' }
else
    subj_list{1} = subj;
end
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')

mtx_sz = 40

% NC -> HC
reg1_reg2_onset_cond2 = cell(1,mtx_sz);
reg1_reg2_onset_cond3 = cell(1,mtx_sz);

% HC -> NC
reg2_reg1_onset_cond2 = cell(1,mtx_sz);
reg2_reg1_onset_cond3 = cell(1,mtx_sz);

i1 = 0;
for iFpass = 1:size(fpass_list,1)
    fpass = [fpass_list(iFpass,:)];
        lock  = locks{1}
        time1 = times1{1}
        time2 = times2{1}
        cd (['/mnt/yassamri/iEEG/sandra/PTE_results/' fn_ext '/' freq_range '/' phase '/' lock '/' num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1) '_' num2str(time2) 'sec'])
        
        for iReg1 = 1:length(reg1_list)
            
            % get NC site
            reg1_name = reg1_list{iReg1};
            % get HC site
            reg2_name = reg2_list{1};
            
            
            % get all data from each subj
            for sub_counter = 1:length(subj_list)
                
                if isfile([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2.mat' ])
                   
                    % load lure -
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2' ])
                    i1 = i1+1;
                    % NC -> HC
                    reg1_reg2_onset_cond2{i1} = ch1_to_ch2;
                    % HC -> NC
                    reg2_reg1_onset_cond2{i1} = ch2_to_ch1;
                    clear ch1_to_ch2 ch2_to_ch1
                    
                    % load lure +
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3' ])
                    % NC -> HC
                    reg1_reg2_onset_cond3{i1} = ch1_to_ch2;
                    % HC -> NC
                    reg2_reg1_onset_cond3{i1} = ch2_to_ch1;                    
                    clear ch1_to_ch2 ch2_to_ch1
                    
                end      
            end
        end
    end
%%

reg1_reg2_onset_cond3_pooled = cat(2,reg1_reg2_onset_cond3{:});
reg2_reg1_onset_cond3_pooled = cat(2,reg2_reg1_onset_cond3{:});

reg1_reg2_onset_cond2_pooled = cat(2,reg1_reg2_onset_cond2{:});
reg2_reg1_onset_cond2_pooled = cat(2,reg2_reg1_onset_cond2{:});
%%
figure
subplot(1,2,1)
h1 = histogram(reg1_reg2_onset_cond3_pooled,'normalization', 'probability','FaceColor', 'm', 'BinWidth', 0.01 )
hold on
h2=histogram(reg1_reg2_onset_cond2_pooled,'normalization', 'probability','FaceColor', 'b','BinWidth', 0.01  )
ylabel('counts')
xlabel('NC -> HC PTE')
legend({'lure+', 'lure-'})
set (gca, 'FontSize', 14, 'FontWeight', 'bold')

subplot(1,2,2)
histogram(reg2_reg1_onset_cond3_pooled,'normalization', 'probability', 'FaceColor', 'm','BinWidth', 0.01 )
hold on
histogram(reg2_reg1_onset_cond2_pooled,'normalization', 'probability','FaceColor', 'b','BinWidth', 0.01 )
ylabel('counts')
xlabel('HC -> NC PTE')
legend({'lure+', 'lure-'})
suptitle([reg1_name ' vs. ' reg2_name])
set (gca, 'FontSize', 14, 'FontWeight', 'bold')
print('-clipboard','-dbitmap')
