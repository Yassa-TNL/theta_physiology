clear all;close all;clc
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
exp_type      = 'encoding' %encoding, tuning_correct
reg           = 'CA3'
lock          = 'onset'   %onset response
ref           = 'LM'
fs            = 500;
if strcmp('encoding', exp_type);
    cond_num=2;    titles = {'lure+','lure-'}
else cond_num = 4; titles = {'repeat', 'lure-', 'lure+', 'new'}
end

[subj_list]   = get_subj_list(reg);

% init
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock ])
load('subj39spectrograms.mat')
cond1 = nan(size(HC_cond2,1), size(HC_cond2,2), length(subj_list));
cond2 = nan(size(HC_cond2,1), size(HC_cond2,2), length(subj_list));
cond3 = nan(size(HC_cond2,1), size(HC_cond2,2), length(subj_list));
cond4 = nan(size(HC_cond2,1), size(HC_cond2,2), length(subj_list));

if strcmp('response', lock)
    data_length = 0;
elseif strcmp('onset', lock)
    data_length = size(fro_cond1,2);
    pre_stim = 0.5; post_stim = 2.0; 
end

% get power values per subj
[cond1,cond2,cond3,cond4] = get_power_subj(subj_list,exp_type,reg, ref, lock);

cond1_grp = nanmean(cond1,3);
cond2_grp = nanmean(cond1,3);
cond3_grp = nanmean(cond1,3);
cond4_grp = nanmean(cond1,3);

% generate cluster mtx
if strcmp('encoding', exp_type)
    time_range = [.3*fs (pre_stim+post_stim)*fs];
    current_pre_stim  = -.2;
    current_post_stim = 2;
else strcmp('tuning_correct', exp_type) && strcmp('onset', lock)
    time_range = [.3*fs 1.5*fs]; 
    current_pre_stim =-.2;
    current_post_stim= 1; 
end
minfreq = 3; maxfreq = 200;
[freq] = get_freq(fs,minfreq, maxfreq);
tickmarks = 1:21:length(freq);
%%
cluster = zeros(size(nanmean(cond1,3)));
range   ='delta';
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/cluster_matrices'])
load([reg '_' lock '_' exp_type '_groupcluster_' range '.mat']);
cluster(:,time_range(1):time_range(2),:) = zmapthresh;
figure;imagesc(cluster)
% group plots w outline
mx = .25
mn = -.25

figure
suptitle(reg)


for i = 1:cond_num
    if     i ==1; cond = cond1; elseif i ==2; cond = cond2;
    elseif i ==3; cond = cond3; elseif i ==4; cond = cond4; end
    subplot (cond_num,1,i)
    contourf(linspace(-pre_stim,post_stim,size(cond,2)), 1:length(freq), mean(cond,3), 150,'linecolor','none');
    if (strcmp('encoding', exp_type) && i==1) || (strcmp('retrieval', exp_type) && i==3)
        hold on;
        contour(linspace(-pre_stim,post_stim,size(cond,2)), 1:length(freq), cluster,[-1 1],'linecolor','k','LineWidth',3)
    end
    set(gca,'YDir','reverse', 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    colormap jet
    caxis([mn mx])
    xlabel('time (s)')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title(titles{i})
    set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
end

%% plot w/ contour
mx = .2
mn = -.2
cond        = cond1(:,time_range(1):time_range(2),:); % chage to cond1 for encoding
cluster_adj = cluster(:,time_range(1):time_range(2),:);
figure
contourf(linspace(current_pre_stim,current_post_stim,size(cond,2)), 1:length(freq), mean(cond,3), 150,'linecolor','none')
hold on
contour(linspace(current_pre_stim,current_post_stim,size(cond,2)), 1:length(freq), cluster_adj,[-1 1],'linecolor','k','LineWidth',3)
set(gca,'YDir','reverse', 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
caxis([mn mx])
colorbar('Ticks',[mn mx])
colormap jet
xlabel('time (s)')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
xticks([0 1 ])
colorbar('Ticks',[mn mx],...
  'TickLabels',[mn mx])
print('-clipboard','-dbitmap')

%% make bar plots and run stats (single value per subj)

for iSubj = 1:length(subj_list)
     if strcmp('encoding', exp_type)
        temp2 = [];
        temp2 = (cond1(:,:,iSubj));
        cond2_grp_val(iSubj) = [nanmean(temp2(logical(cluster)) )];
        
        temp3 = [];
        temp3 = (cond2(:,:,iSubj));
        cond3_grp_val(iSubj) = [nanmean(temp3(logical(cluster)) )];
     elseif strcmp('tuning_correct', exp_type)
        temp1 = [];
        temp1 = (cond1(:,:,iSubj));
        cond1_grp_val(iSubj) = [nanmean(temp1(logical(cluster)) )];
        
        temp2 = [];
        temp2 = (cond2(:,:,iSubj));
        cond2_grp_val(iSubj) = [nanmean(temp2(logical(cluster)) )];
        
        temp3 = [];
        temp3 = (cond3(:,:,iSubj));
        cond3_grp_val(iSubj) = [nanmean(temp3(logical(cluster)) )];
        
        temp4 = [];
        temp4 = (cond4(:,:,iSubj));
        cond4_grp_val(iSubj) = [nanmean(temp4(logical(cluster)) )];        
     end
end

if strcmp('encoding', exp_type)
    bar_vector_mn1  = [nanmean(cond2_grp_val) nanmean(cond3_grp_val) ];
    bar_vector_std1 = [ nanstd(cond2_grp_val,0,2) nanstd(cond3_grp_val,0,2)];
    labels =  {'lure +' 'lure -'};
    
    % encoding stats
    reps = 1000;
    adata = cond3_grp_val;
    bdata = cond2_grp_val;
    p2    = permutation_paired(adata, bdata, reps)
    
else
    bar_vector_mn1 = [nanmean(cond1_grp_val) nanmean(cond2_grp_val) nanmean(cond3_grp_val) nanmean(cond4_grp_val)];
    bar_vector_std1 = [nanstd(cond1_grp_val,0,2) nanstd(cond2_grp_val,0,2) nanstd(cond3_grp_val,0,2) nanstd(cond4_grp_val,0,2)];
    labels = {'repeat' 'lure -' 'lure +' 'new'};
    
end

% plot
figure
hold on
bar([bar_vector_mn1])         
errorbar(1:length(bar_vector_mn1),bar_vector_mn1,bar_vector_std1./sqrt(length(subj_list)), 'rx')
set(gca, 'XTick', 1:length(bar_vector_mn1), 'XTickLabel', labels,'XTickLabelRotation',45)
ylabel([reg ' power'])
set(gca, 'FontSize', 16, 'FontWeight', 'bold') 
title([reg ' ' range ])
% print('-clipboard','-dbitmap')
% print('-clipboard','-depsc')



