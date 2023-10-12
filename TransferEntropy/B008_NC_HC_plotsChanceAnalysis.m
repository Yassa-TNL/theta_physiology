
% description: script which will access the null distributions and observed
% data across subjects and regions then plot.
% note: as it stands, script can only be run for chance analysis for lure+
% condition NC-->HC (non-normalized direc) during retrievial and HC--> NC during encoding. 
% Would need to update script for the three other conditions, other directions, and normalized direction.

clear all;close all;clc

lock             = 'onset';
phase            = 'encoding';
%phase            = 'retrieval';
reg1_list   = {'OFC' 'FRO' 'TEMP'};
reg2_list   = {'HC'};
nPerms           = 1000;
sig_chans        ='';
group_plots      = 'yes';
subj             = '85';
freq_range       = 'deltatheta';
fpass            = [4 5];
fn_ext           = 'cue_resp'; %
plot_4_conds     = '';
times1           = {'0' '0' '0.5'};
times2           = {'2' '1' '1.5'};
condCounterList  =2; % default is only compare lure+ vs. lure- unles plot 4 conds
if strcmp('encoding',phase)
    time1 = times1{1};
    time2 = times2{1};
    lock = 'onset';
elseif strcmp('retrieval',phase) && strcmp('onset', lock)
    time1 = times1{2};
    time2 = times2{2};
elseif strcmp('retrieval',phase) && strcmp('response', lock)
    time1 = times1{3};
    time2 = times2{3};
end

addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
cd (['/mnt/yassamri/iEEG/sandra/PTE_results/' fn_ext '/' freq_range '/' phase '/' lock '/' num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1) '_' num2str(time2) 'sec'])

mtx_sz = 40;
% cond 2: normalized lure -
reg1_reg2_conda = cell(1,mtx_sz);
% cond 3: normalized lure+
reg1_reg2_condb = cell(1,mtx_sz);
% cond 3: unidirec lure+
ch1_to_ch2_condb = cell(1,mtx_sz); % for retrieval
ch2_to_ch1_condb = cell(1,mtx_sz); % for encoding

% cond1: normalized repeat
reg1_reg2_condc = cell(1,mtx_sz);
% cond 4: normalized new
reg1_reg2_condd = cell(1,mtx_sz);

% unidirec null
ch1_to_ch2_Mean_NullDist_condb = cell(1,mtx_sz);
ch2_to_ch1_Mean_NullDist_condb = cell(1,mtx_sz);

% norm direc null
ch2_to_ch1_norm_Mean_nullDist_condb = cell(1,mtx_sz);

i = 0;
for iReg1 = 1:length(reg1_list)
    
    % get NC site
    reg1_name = reg1_list{iReg1};
    reg2_name = reg2_list{1};
    
    if strcmp('yes',group_plots)
        subj_list = {'39' '44' '57' '63' '66'  '84' '85' '87' };
    else
        subj_list{1} = subj;
    end
    % get all data from each subj
    for sub_counter = 1:length(subj_list)
        
        % get null data
        % load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3Chance' sig_chans '.mat'])
        if strcmp('retrieval',phase)
            if isfile([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3Chance' sig_chans '.mat'])
                i = i+1;
                
                % load chance for lure+ per subject and area
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3Chance' sig_chans ])
                
                ch1_to_ch2_Mean_NullDist_condb{i} = ch1_to_ch2';clear ch1_to_ch2 % NC ==> HC
                ch2_to_ch1_Mean_NullDist_condb{i} = ch2_to_ch1';clear ch2_to_ch1 % NC ==> HC
                ch2_to_ch1_norm_Mean_nullDist_condb{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm % NC ==> HC
                
                % load true lure +
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3' sig_chans ])
                reg1_reg2_condb{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                ch1_to_ch2_condb{i} = ch1_to_ch2';clear ch1_to_ch2
                
                
                % can add same section as above here for lure- if want to
                % do stats
                
                % update this section for 2 other conds if want to do stats
                if strcmp('yes',plot_4_conds)
                    condCounterList = [1 2 3];
                    
                    % load repeat
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' sig_chans ])
                    reg1_reg2_condc{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                    
                    % load new
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond4' sig_chans ])
                    reg1_reg2_condd{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                end
            end
        else
            if isfile([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1Chance' sig_chans '.mat'])
                i = i+1;
                
                % load chance for lure+ per subject and area
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1Chance' sig_chans ])
                
                ch1_to_ch2_Mean_NullDist_condb{i} = ch1_to_ch2';clear ch1_to_ch2 % NC ==> HC
                ch2_to_ch1_Mean_NullDist_condb{i} = ch2_to_ch1';clear ch2_to_ch1 % NC ==> HC
                ch2_to_ch1_norm_Mean_nullDist_condb{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm % NC ==> HC
                
                % load lure +
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' sig_chans ])
                reg1_reg2_condb{i}  = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                ch2_to_ch1_condb{i} = ch2_to_ch1';clear ch2_to_ch1
                
                % load lure -: update if you want to do stats
                
            end
        end
    end
end
% TRUE obsv: pool observations across subjects - pairs of regions for data
reg1_reg2_conda = cat(1,reg1_reg2_conda{:}); % lure-
reg1_reg2_condb = cat(1,reg1_reg2_condb{:}); % lure+
reg1_reg2_condc = cat(1,reg1_reg2_condc{:}); % repeat
reg1_reg2_condd = cat(1,reg1_reg2_condd{:}); % new

% Null obsv: non-norm PTE - unidirec
ch1_to_ch2_condb = cat(1,ch1_to_ch2_condb{:}); % retrieval lure+
ch2_to_ch1_condb = cat(1,ch2_to_ch1_condb{:}); % encoding lure+

% pool null dist across subjects - 1000 X pairs of regions for data
ch1_to_ch2_Mean_NullDist_condb_Cat = cat(1,ch1_to_ch2_Mean_NullDist_condb{:});
ch2_to_ch1_Mean_NullDist_condb_Cat = cat(1,ch2_to_ch1_Mean_NullDist_condb{:});

if strcmp('retrieval',phase) % plot for retrieval lure+  
    %get a single mean for true observation and a mean for null obs, but retain the 1000 data points for that mean
    ch1_to_ch2_Mean_NullDist_condb_Cat_Final = nanmean(ch1_to_ch2_Mean_NullDist_condb_Cat,1);
    
    % see where the observed lies along this dist
    pval = sum(ch1_to_ch2_Mean_NullDist_condb_Cat_Final>nanmean(ch1_to_ch2_condb))/nPerms;
    figure;histogram(ch1_to_ch2_Mean_NullDist_condb_Cat_Final, 'normalization', 'probability')
    hold on; y=ylim;
    line( [nanmean(ch1_to_ch2_condb) nanmean(ch1_to_ch2_condb)],[0 y(2)], 'color', 'r')
    
else % plot for encoding lure+
    ch2_to_ch1_Mean_NullDist_condb_Cat_Final = nanmean(ch2_to_ch1_Mean_NullDist_condb_Cat,1);
    
    % see where the observed lies along this dist
    pval = sum(ch2_to_ch1_Mean_NullDist_condb_Cat_Final>nanmean(ch2_to_ch1_condb))/nPerms;
    figure;histogram(ch2_to_ch1_Mean_NullDist_condb_Cat_Final, 'normalization', 'probability')
    hold on; y=ylim;
    line( [nanmean(ch2_to_ch1_condb) nanmean(ch2_to_ch1_condb)],[0 y(2)], 'color', 'r')
end