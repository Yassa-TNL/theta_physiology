% Script: generates figures of pwoer traces as a functin of a desired freq.
% Runs perm testing stats of traces, using channels as observations.
% Generates also group matrices, using subj as observations (last cell).

clear all;close all;clc
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
exp_type      = 'encoding' %encoding, tuning_correct
reg           = 'OFC_FRO_TEMP_CING_INS_EC_significant' %'OFC_FRO_TEMP_significant' 'OFC_FRO_TEMP_CING_INS_EC_significant' 'HC' 'NC'
lock          = 'onset'   %onset response
lock_ext      = '' % to get 0-2 seconds
ref           = 'LM'
fs            = 500;
fn_ext        = '_cue_responseiveyes'
baseline      = [ '_cond_spec_prestim' fn_ext] % '': entire recording, '_prestim'

desired_freq  = ''
desired_freq_lo  = 4
desired_freq_hi  = 5
plots_4_conds    =''

if strcmp('encoding', exp_type)
    cond_num=2;    titles = {'lure+','lure-'}
    xtick_vals = [-0.2 0 2];
else
    cond_num = 4; titles = {'repeat', 'lure-', 'lure+', 'new'}
end

subj_list = {'39' '44' '57' '63' '66' '84' '85' '87'}  ;

if strcmp('response', lock)
    pre_stim           = 1.5;
    post_stim          = .5;
    pre_stim_minus_edge  = 1.5-200/fs;
    post_stim_minus_edge = 0.5-200/fs;
    time_range = [0.1*fs 1.2*fs];
    current_pre_stim   = -1;
    current_post_stim  = 0.1;
    xtick_vals         = [-1  -0.5 0 0.1];
elseif strcmp('onset', lock)
    pre_stim = 0.5; post_stim = 2.0;
    if strcmp('encoding', exp_type)
        time_range = [.3*fs (pre_stim+post_stim)*fs];
        current_pre_stim  = -.2;
        current_post_stim = 2;
        xtick_vals = [current_pre_stim 0 1 2];
    else
        time_range = [.3*fs 1.5*fs];
        current_pre_stim  =-.2;
        current_post_stim = 1; % was 1;
        xtick_vals = [-0.2 0 1];
    end
% elseif strcmp('onset2', lock)
%     pre_stim = 0.5; post_stim = 2.6;
%     time_range = [.3*fs (pre_stim+2)*fs]; % was [301 1301];
%     current_pre_stim =-.2;
%     current_post_stim= 2;
%     xtick_vals = [current_pre_stim 0 1 2];    
end
if strcmp('yes',lock_ext) && strcmp('onset',lock)
    time_range(2)=(pre_stim+2)*fs;
    current_pre_stim =-.2;
    current_post_stim= 2;
    xtick_vals = [-0.5 0 0.5 1 1.5 2];
end
minfreq = 3; maxfreq = 200;
[freq]  = get_freq(fs,minfreq, maxfreq);
% tickmarks = 1:21:length(freq);

[cond1a,cond2a,cond3a,cond4a] = get_power_subj_elecs(subj_list,exp_type,reg,ref,baseline,lock);


if strcmp('tuning_correct',exp_type)
    cond1 = cat(3,cond1a{:});
    cond4 = cat(3,cond4a{:});
    
    trace1 = squeeze(nanmean(cond1((freq>desired_freq_lo)&(freq<desired_freq_hi), time_range(1):time_range(2), :)))';
    trace4 = squeeze(nanmean(cond4((freq>desired_freq_lo)&(freq<desired_freq_hi), time_range(1):time_range(2), :)))';
end
cond2 = cat(3,cond2a{:});
cond3 = cat(3,cond3a{:});

% Matrices for traces
trace2 = squeeze(nanmean(cond2((freq>desired_freq_lo)&(freq<desired_freq_hi), time_range(1):time_range(2), :),1))';
trace3 = squeeze(nanmean(cond3((freq>desired_freq_lo)&(freq<desired_freq_hi), time_range(1):time_range(2), :),1))';

if strcmp('onset',lock)
reset = nanmean(nanmean(trace2(:,1:-1*current_pre_stim*fs),1)); 
trace2=trace2-reset;
reset = nanmean(nanmean(trace3(:,1:-1*current_pre_stim*fs),1)); 
trace3=trace3-reset;
end

% stats pre vs. post
% condb = nanmean(trace2(:,(-1*current_pre_stim*fs)+1:(-1*current_pre_stim+1)*fs),2);
% conda = nanmean(trace3(:,(-1*current_pre_stim*fs)+1:(-1*current_pre_stim+1)*fs),2);
% reps = 1000; p1 = permutation_paired(conda, condb, reps);

% run stats on 0 - 2 seconds
%
[zmap,zmapthresh,zmapthresh_for_plot] = permutation_testing_vector(trace2,trace3, 1000);

if nansum(zmapthresh_for_plot)>0
   disp('sig!!') 
end

% plot traces and stats
figure;
hold on;
stdshade(trace3,.1,'m',linspace(current_pre_stim,current_post_stim, size(trace3,2)),[] ,[], []);
plot(linspace(current_pre_stim,current_post_stim, size(trace3,2)),nanmean(trace3,1), 'm', 'LineWidth', 2);
stdshade(trace2,.1,'b',linspace(current_pre_stim,current_post_stim, size(trace2,2)),[] ,[], []);
plot(linspace(current_pre_stim,current_post_stim, size(trace2,2)),nanmean(trace2,1), 'b', 'LineWidth', 2);
if strcmp('tuning_correct',exp_type)&&strcmp('yes',plots_4_conds)
    reset = nanmean(nanmean(trace1(:,1:-1*current_pre_stim*fs),1));
    trace1=trace1-reset;
    stdshade(trace1,.1,'r',linspace(current_pre_stim,current_post_stim, size(trace1,2)),[] ,[], []);
    h3 = plot(linspace(current_pre_stim,current_post_stim, size(trace3,2)),nanmean(trace1,1), 'r', 'LineWidth', 2);
    reset = nanmean(nanmean(trace4(:,1:-1*current_pre_stim*fs),1));
    trace4=trace4-reset;
    stdshade(trace4,.1,'k',linspace(current_pre_stim,current_post_stim, size(trace4,2)),[] ,[], []);
    h4 = plot(linspace(current_pre_stim,current_post_stim, size(trace2,2)),nanmean(trace4,1), 'k', 'LineWidth', 2);
end
scale = max([nanmean(trace3,1) nanmean(trace2,1)])+0.05;
hold on
% if strcmp('response',lock)
% plot(linspace(current_pre_stim,current_post_stim, size(trace2,2)),scale*zmapthresh_for_plot, 'LineWidth', 2, 'Color', 'k')
% else
% plot(linspace(-1*current_pre_stim,current_post_stim, size(trace2,2)),scale*zmapthresh_for_plot, 'LineWidth', 2, 'Color', 'k')
% end
xlim([current_pre_stim current_post_stim]);

%xlabel('Time (s)')
%title(reg )
%line([0 0],ylim,'color','k')
%ylim([-0.5 0.5])
y=ylim;
if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
end
print('-clipboard','-dbitmap')
cd('/mnt/yassamri/iEEG/sandra/GroupFigures')
%% gamma bar plots
cond1_grp_val=nanmean(trace1(:,1.2*fs:1.7*fs),2);
cond2_grp_val=nanmean(trace2(:,1.2*fs:1.7*fs),2);
cond3_grp_val=nanmean(trace3(:,1.2*fs:1.7*fs),2);
cond4_grp_val=nanmean(trace4(:,1.2*fs:1.7*fs),2);
bar_vector_mn = [nanmean(cond1_grp_val)  nanmean(cond2_grp_val)  ...
    nanmean(cond3_grp_val)  nanmean(cond4_grp_val)  ]
bar_vector_sem = [nanstd(cond1_grp_val) nanstd(cond2_grp_val)...
                  nanstd(cond3_grp_val) nanstd(cond4_grp_val)]
figure;hold on
bar([bar_vector_mn])         
errorbar(1:length(bar_vector_mn),bar_vector_mn,bar_vector_sem, 'kx')
set(gca, 'XTick', 1:length(bar_vector_mn), 'XTickLabel', {'repeat+' ,'lure+', 'lure-' ,'new'},'XTickLabelRotation',45)
y=ylim;
if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
end
%ylim([-0.1 0.045]);
ylabel([reg ' power'])

    % retrieval stats
reps  = 1000;
adata = cond3_grp_val; 
bdata = cond1_grp_val;
p1     = permutation_paired(adata, bdata, reps)

adata = cond3_grp_val; 
bdata = cond4_grp_val;
p3    = permutation_paired(adata, bdata, reps)

alpha = .05
pvals = [p1 p3];
[p_fdr, p_masked] = fdr( pvals, alpha)
[pvals p_fdr]
%% gamma vs delta traces
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/LM_reref/tuning_correct_onset2'])

% gamma3 = trace3; gamma2 = trace2;save([reg '_gamma'],'gamma2', 'gamma3')
% deltatheta3 = trace3; deltatheta2 = trace2;save([reg '_deltatheta'],'deltatheta2', 'deltatheta3')
load([reg '_deltatheta.mat'])
load([reg '_gamma.mat'])

figure;hold on;
stdshade(deltatheta - nanmean(nanmean(deltatheta(:,1:(fs*-current_pre_stim)))),.1,'c',linspace(current_pre_stim,current_post_stim, size(trace3,2)),[] ,[], []);
h1 = plot(linspace(current_pre_stim,current_post_stim, size(trace3,2)),nanmean(deltatheta,1)-nanmean(nanmean(deltatheta(:,1:(fs*-current_pre_stim)))), 'c', 'LineWidth', 2);

stdshade(gamma,.1,'r',linspace(current_pre_stim,current_post_stim, size(trace3,2)),[] ,[], []);
h2 = plot(linspace(current_pre_stim,current_post_stim, size(trace3,2)),nanmean(gamma,1), 'r', 'LineWidth', 2);

xlim([current_pre_stim current_post_stim]);
xticks(xtick_vals)
ylim([-0.15 0.2])
y = ylim;
line([0 0],ylim,'color','k')
yticks([y(1) 0 y(2)])
print('-clipboard','-dbitmap')



%% bar plots first and second second for retrieval, group then idiv subj

% group
figure;hold on
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/LM_reref/tuning_correct_onset2'])
load([reg '_deltatheta.mat'])
load([reg '_gamma.mat'])
time_range1 = [(-1*current_pre_stim)*fs (-1*current_pre_stim+1)*fs]; % first second
time_idx = time_range1
cond1a_vals = nanmean(deltatheta3(:,time_idx(1):time_idx(2)),2)-nanmean(deltatheta2(:,time_idx(1):time_idx(2)),2);
cond2a_vals = nanmean(gamma3(:,time_idx(1):time_idx(2)),2)-nanmean(gamma2(:,time_idx(1):time_idx(2)),2);


time_range2 = [(-1*current_pre_stim+1)*fs (-1*current_pre_stim+2)*fs]; % second second
time_idx = time_range2
cond1b_vals = nanmean(deltatheta3(:,time_idx(1):time_idx(2)),2)-nanmean(deltatheta2(:,time_idx(1):time_idx(2)),2);
cond2b_vals = nanmean(gamma3(:,time_idx(1):time_idx(2)),2)-nanmean(gamma2(:,time_idx(1):time_idx(2)),2);

bar_vector_mn = [nanmean(cond1a_vals) nanmean(cond2a_vals) nanmean(cond1b_vals) nanmean(cond2b_vals) ]
bar_vector_sem = [nanstd(cond1a_vals)/sqrt(length(cond1a_vals)) nanstd(cond2a_vals)/sqrt(length(cond2a_vals))...
    nanstd(cond1b_vals)/sqrt(length(cond1b_vals)) nanstd(cond2b_vals)/sqrt(length(cond2b_vals))]

cd('/mnt/yassamri/iEEG/sandra/GroupFigures')

bar([bar_vector_mn])         
errorbar(1:length(bar_vector_mn),bar_vector_mn,bar_vector_sem, 'kx')
set(gca, 'XTick', 1:length(bar_vector_mn), 'XTickLabel', {'3-5 Hz' ,'>45 Hz', '3-5 Hz' ,'>45 Hz'},'XTickLabelRotation',45)
%ylim([-0.1 0.045]);
ylabel([reg ' power'])

% stats
reps  = 1000;
adata = cond1a_vals; 
bdata = cond2a_vals;
p1     = permutation_paired(adata, bdata, reps)

reps  = 1000;
adata = cond1a_vals; 
bdata = cond1b_vals;
p1     = permutation_paired(adata, bdata, reps)

reps  = 1000;
adata = cond1b_vals; 
bdata = cond2b_vals;
p1     = permutation_paired(adata, bdata, reps)

reps  = 1000;
adata = cond1b_vals; 
bdata = cond2b_vals;
p1     = permutation_paired(adata, bdata, reps)

clear cond1a cond2a cond1b cond2b
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock  ])
for iSubj = 1:length(subj_list)
    load(['subj' subj_list{iSubj} 'spectrograms_all_trials_across_channels_entire_recording'])
    
    if  strcmp('OFC',reg)
        data1 = OFC_cond1;
        data2 = OFC_cond2;
        data3 = OFC_cond3;
        data4 = OFC_cond4;
        
    elseif strcmp('FRO',reg)
        data1 = fro_cond1;
        data2 = fro_cond2;
        data3 = fro_cond3;
        data4 = fro_cond4;
        
    elseif strcmp('TEMP',reg)
        data1 = temp_cond1;
        data2 = temp_cond2;
        data3 = temp_cond3;
        data4 = temp_cond4;
        
    elseif strcmp('CING',reg)
        data1 = cing_cond1;
        data2 = cing_cond2;
        data3 = cing_cond3;
        data4 = cing_cond4;
        
    elseif strcmp('ins',reg)
        data1 = ins_cond1;
        data2 = ins_cond2;
        data3 = ins_cond3;
        data4 = ins_cond4;
        
    elseif strcmp('EC',reg)
        data1 = EC_cond1;
        data2 = EC_cond2;
        data3 = EC_cond3;
        data4 = EC_cond4;
        
    elseif strcmp('HC',reg)
        data1 = HC_cond1;
        data2 = HC_cond2;
        data3 = HC_cond3;
        data4 = HC_cond4;
        
    elseif strcmp('CA1',reg)
        data1 = CA1_cond1;
        data2 = CA1_cond2;
        data3 = CA1_cond3;
        data4 = CA1_cond4;
        
    elseif strcmp('CA3',reg)
        data1 = CA3_cond1;
        data2 = CA3_cond2;
        data3 = CA3_cond3;
        data4 = CA3_cond4;
        
    elseif strcmp('NC',reg)
        data1 = NC_cond1;
        data2 = NC_cond2;
        data3 = NC_cond3;
        data4 = NC_cond4;
    end

    % lure+ delta first sec, gamma first, dela second, gamma second
    time_range_trl = [pre_stim*fs+1:(pre_stim+1)*fs];
    cond1a(iSubj) = nanmean(squeeze(nanmean(nanmean(data3(freq>3 & freq<5, time_range_trl, :),2),1)))-...
        nanmean(squeeze(nanmean(nanmean(data2(freq>3 & freq<5, time_range_trl, :),2),1)));
    cond2a(iSubj) =  nanmean(squeeze(nanmean(nanmean(data3(freq>45, time_range_trl, :),2),1)))-...
        nanmean(squeeze(nanmean(nanmean(data2(freq>45, time_range_trl, :),2),1)));
    clear time_range_trl
    
    time_range_trl = [(pre_stim+1)*fs+1:(pre_stim+2)*fs];
    cond1b(iSubj) =  nanmean(squeeze(nanmean(nanmean(data3(freq>3 & freq<5, time_range_trl, :),2),1)))-...
        nanmean(squeeze(nanmean(nanmean(data2(freq>3 & freq<5, time_range_trl, :),2),1)));
    cond2b(iSubj) =  nanmean(squeeze(nanmean(nanmean(data3(freq>45, time_range_trl, :),2),1)))-...
        nanmean(squeeze(nanmean(nanmean(data2(freq>45, time_range_trl, :),2),1)));
    
end
          
colorList = {'ro', 'co', 'go', 'yo', 'mo', 'ko', 'bo', '-ko'};
color     = {'r', 'c', 'g', 'y', 'm', 'k', 'b', 'k'};
for iSubj = 1:length(subj_list)
    plot( [cond1a(iSubj) cond2a(iSubj) cond1b(iSubj) cond2b(iSubj)], colorList{iSubj}, 'MarkerFaceColor', color{iSubj}, 'LineWidth', 1)
end



%% 4 conds second second at desired freq, group then indiv
close all
time_range2 = [(-1*current_pre_stim+1)*fs (-1*current_pre_stim+2)*fs]; % second second
time_idx = time_range2
cond1_vals = nanmean(trace1(:,time_idx(1):time_idx(2)),2);
cond2_vals = nanmean(trace2(:,time_idx(1):time_idx(2)),2);
cond3_vals = nanmean(trace3(:,time_idx(1):time_idx(2)),2);
cond4_vals = nanmean(trace4(:,time_idx(1):time_idx(2)),2);

bar_vector_mn = [nanmean(cond1_vals) nanmean(cond2_vals) nanmean(cond3_vals) nanmean(cond4_vals) ]
bar_vector_sem = [nanstd(cond1_vals)/sqrt(length(cond1_vals)) nanstd(cond2_vals)/sqrt(length(cond2_vals))...
    nanstd(cond3_vals)/sqrt(length(cond3_vals)) nanstd(cond4_vals)/sqrt(length(cond4_vals))]
  
figure
hold on
bar([bar_vector_mn],'EdgeColor', 'none')         
errorbar(1:length(bar_vector_mn),bar_vector_mn,bar_vector_sem, 'kx')
set(gca, 'XTick', [])

% run stats
reps  = 1000;
adata = cond1_vals; 
bdata = cond3_vals;
p1     = permutation_paired(adata, bdata, reps)

adata = cond2_vals; 
bdata = cond3_vals;
p2    = permutation_paired(adata, bdata, reps)

adata = cond4_vals; 
bdata = cond3_vals;
p3    = permutation_paired(adata, bdata, reps)

alpha = .05
pvals = [p1 p2 p3];
[p_fdr, p_masked] = fdr( pvals, alpha)
pvals
p_fdr

cond1_subj = nan(1,length(subj_list));
cond2_subj = nan(1,length(subj_list));
cond3_subj = nan(1,length(subj_list));
cond4_subj = nan(1,length(subj_list));

cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock  ])
for iSubj = 1:length(subj_list)
    load(['subj' subj_list{iSubj} 'spectrograms_all_trials_across_channels_entire_recording'])
    
    if  strcmp('OFC',reg)
        data1 = OFC_cond1;
        data2 = OFC_cond2;
        data3 = OFC_cond3;
        data4 = OFC_cond4;
        
    elseif strcmp('FRO',reg)
        data1 = fro_cond1;
        data2 = fro_cond2;
        data3 = fro_cond3;
        data4 = fro_cond4;
        
    elseif strcmp('TEMP',reg)
        data1 = temp_cond1;
        data2 = temp_cond2;
        data3 = temp_cond3;
        data4 = temp_cond4;
        
    elseif strcmp('CING',reg)
        data1 = cing_cond1;
        data2 = cing_cond2;
        data3 = cing_cond3;
        data4 = cing_cond4;
        
    elseif strcmp('ins',reg)
        data1 = ins_cond1;
        data2 = ins_cond2;
        data3 = ins_cond3;
        data4 = ins_cond4;
        
    elseif strcmp('EC',reg)
        data1 = EC_cond1;
        data2 = EC_cond2;
        data3 = EC_cond3;
        data4 = EC_cond4;
        
    elseif strcmp('HC',reg)
        data1 = HC_cond1;
        data2 = HC_cond2;
        data3 = HC_cond3;
        data4 = HC_cond4;
        
    elseif strcmp('CA1',reg)
        data1 = CA1_cond1;
        data2 = CA1_cond2;
        data3 = CA1_cond3;
        data4 = CA1_cond4;
        
    elseif strcmp('CA3',reg)
        data1 = CA3_cond1;
        data2 = CA3_cond2;
        data3 = CA3_cond3;
        data4 = CA3_cond4;
        
    elseif strcmp('NC',reg)
        data1 = NC_cond1;
        data2 = NC_cond2;
        data3 = NC_cond3;
        data4 = NC_cond4;
    end

    % lure+ delta first sec, gamma first, dela second, gamma second
  
    time_range_trl = [(pre_stim+1)*fs+1:(pre_stim+2)*fs];
    cond1_subj(iSubj) =  nanmean(squeeze(nanmean(nanmean(data1(freq>desired_freq_lo & freq<desired_freq_hi, time_range_trl, :),2),1)));
    cond2_subj(iSubj) =  nanmean(squeeze(nanmean(nanmean(data2(freq>desired_freq_lo & freq<desired_freq_hi, time_range_trl, :),2),1)));
    cond3_subj(iSubj) =  nanmean(squeeze(nanmean(nanmean(data3(freq>desired_freq_lo & freq<desired_freq_hi, time_range_trl, :),2),1)));
    cond4_subj(iSubj) =  nanmean(squeeze(nanmean(nanmean(data4(freq>desired_freq_lo & freq<desired_freq_hi, time_range_trl, :),2),1)));
    
end
          
colorList = {[1 0 1], [0 1 1], [1 0 0], [0 1 0], [0 0 1], [0.5 0.1 0.2], [0.1 0.5 0.3], [0.3 0.1 0.5]};
colors = {[1 0 1], [0 1 1], [1 0 0], [0 1 0], [0 0 1], [0.5 0.1 0.2], [0.1 0.5 0.3], [0.3 0.1 0.5]};
for iSubj = 1:length(subj_list)
    plot( [cond1_subj(iSubj) cond2_subj(iSubj) cond3_subj(iSubj) cond4_subj(iSubj)],'o','MarkerFaceColor',[colorList{iSubj}],'MarkerEdgeColor',[colorList{iSubj}])
end
y=ylim;
if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
end
cd('/mnt/yassamri/iEEG/sandra/GroupFigures')
%%
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/cluster_matrices'])
load('HC_onset_encoding_groupcluster_delta.mat')
zmapthresh_for_plot  = nan(size(zmapthresh))
zmapthresh_for_plot(logical(sum(zmapthresh,1)))=1

%% NC vs. HC
figure;hold on; stdshade(trace3,.1,'c',linspace(current_pre_stim,current_post_stim, size(trace3,2)),[] ,[], []);
h1 = plot(linspace(current_pre_stim,current_post_stim, size(trace3,2)),nanmean(trace3,1), 'c', 'LineWidth', 2);

stdshade(trace3,.1,'b',linspace(current_pre_stim,current_post_stim, size(trace3,2)),[] ,[], []);
h2 = plot(linspace(current_pre_stim,current_post_stim, size(trace3,2)),nanmean(trace3,1), 'g', 'LineWidth', 2);

%% save matrices for later bar plots
if strcmp('gamma',desired_freq)
     zmapthresh_for_plot(end-current_post_stim*fs:end ) = nan;
     zmapthresh_for_plot(1:end-current_post_stim*fs) = 1;
end

% retrieval onset:  HC
if strcmp('gamma',desired_freq)
    zmapthresh_for_plot = nan(size(zmapthresh))
    zmapthresh_for_plot(1:(abs(current_pre_stim)+0.5)*fs ) = nan;
    zmapthresh_for_plot((abs(current_pre_stim)+0.5)*fs:end) = 1;
end

% retrieval onset:  NC and CA3
if strcmp('gamma',desired_freq)
    zmapthresh_for_plot(1:(abs(current_pre_stim)+0.5)*fs ) = nan;
end

cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/cluster_matrices'])
if strcmp('response',lock)
    zmapthresh_for_plot(end-current_post_stim*fs:end ) = nan; % NC resposne delta
end

%save([reg '_' lock '_' exp_type '_' ref '_''zmapthresh_for_plotting_mean_after_traces'], 'zmapthresh_for_plot') CA3 save as is

%% make bar plots for response 
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/cluster_matrices'])
load([reg '_' lock '_' exp_type '_' ref '_' desired_freq '_zmapthresh_for_plotting_mean_after_traces'])

cond1_vals = nanmean(trace1(:,~isnan(zmapthresh_for_plot)),2);
cond2_vals = nanmean(trace2(:,~isnan(zmapthresh_for_plot)),2);
cond3_vals = nanmean(trace3(:,~isnan(zmapthresh_for_plot)),2);
cond4_vals = nanmean(trace4(:,~isnan(zmapthresh_for_plot)),2);

bar_vector_mn = [nanmean(cond1_vals) nanmean(cond2_vals) nanmean(cond3_vals) nanmean(cond4_vals)]
bar_vector_sem = [nanmean(cond1_vals)/sqrt(length(cond1_vals)) nanmean(cond2_vals)/sqrt(length(cond2_vals))...
             nanmean(cond3_vals)/sqrt(length(cond3_vals)) nanmean(cond4_vals)/sqrt(length(cond4_vals))]

cd('/mnt/yassamri/iEEG/sandra/GroupFigures')

figure
hold on
bar([bar_vector_mn])         
errorbar(1:length(bar_vector_mn),bar_vector_mn,bar_vector_sem, 'rx')
set(gca, 'XTick', 1:length(bar_vector_mn), 'XTickLabel', titles,'XTickLabelRotation',45)
y = ylim;
ylabel([reg ' power'])
set(gca, 'FontSize', 16, 'FontWeight', 'bold') 

% run stats
reps  = 1000;
adata = cond1_vals; 
bdata = cond3_vals;
p1     = permutation_paired(adata, bdata, reps)

adata = cond2_vals; 
bdata = cond3_vals;
p2    = permutation_paired(adata, bdata, reps)

adata = cond4_vals; 
bdata = cond3_vals;
p3    = permutation_paired(adata, bdata, reps)

alpha = .05
pvals = [p1 p2 p3];
[p_fdr, p_masked] = fdr( pvals, alpha)
pvals
p_fdr;


%%  get power values, obs = subj
% init
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock ])
load('subj39spectrograms.mat')
cond1 = nan(size(HC_cond2,1), size(HC_cond2,2), length(subj_list));
cond2 = nan(size(HC_cond2,1), size(HC_cond2,2), length(subj_list));
cond3 = nan(size(HC_cond2,1), size(HC_cond2,2), length(subj_list));
cond4 = nan(size(HC_cond2,1), size(HC_cond2,2), length(subj_list));

% redefined matrices (freqXtimeXsubj)
[cond1,cond2,cond3,cond4] = get_power_subj(subj_list,exp_type,reg, ref, lock);




