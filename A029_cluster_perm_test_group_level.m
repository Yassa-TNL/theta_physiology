clear all; clc; close all
rmpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/'))
reset = 1% reset to zero 0=no, 1 = yes
exp_type      = 'encoding' %encoding, tuning_correct
smoothgamma   = ''
fn_ext        = ''%'_cue_responseiveyes'%'_cue_responseiveyes'
removeERP     = ''
if strcmp('' ,removeERP)
    baseline      = ['_cond_spec_prestim' fn_ext] %'_entire_recording', '_cond_spec_prestim' fn_ext 
else
    baseline      = [ '_cond_spec_prestim' fn_ext 'removeERP' removeERP] %'_entire_recording'
end
reg           = 'HC' %'CA3' 'CA3_CA1_lure-' 'OFC_FRO_TEMP_significant' 'OFC_FRO_TEMP_CING_INS_EC_significant'
lock          = 'onset'   %onset response
lock_ext      = '' %empty or yes to get two second window
clustercorect = 'maxsize' %maxsum maxsize ws
ref           = 'LM'
fs            = 500;
minfreq       = 3;
maxfreq       = 200;
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')

subj_list = {'39' '44' '57' '63' '66' '84' '85' '87'};
%subj_list = {'39' '66' '84' '87'};% HC path
%subj_list = {'44' '57' '63' '85'};% TEMP path


total_subj = length(subj_list);

% get wavelet params
[freq] = get_freq(fs,minfreq, maxfreq);
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock ])

if strcmp('response', lock)
    pre_stim           = 1.5;
    post_stim          = .5;
    pre_stim_minus_edge  = 1.5-200/fs;
    post_stim_minus_edge = 0.5-200/fs;
    time_range = [0.1*fs 1.2*fs];
    current_pre_stim  = -1;
    current_post_stim = 0.1;
    x_tick_vec = [-1 -0.5 0] ;
elseif strcmp('onset', lock)
    pre_stim = 0.5; post_stim = 2.0;
    if strcmp('encoding', exp_type)
        time_range = [.3*fs (pre_stim+post_stim)*fs];
        current_pre_stim  = -.2;
        current_post_stim = 2;
        x_tick_vec = [0  1 current_post_stim];
    else
        current_pre_stim = -.2;
        current_post_stim= 1; % was 1;
        time_range = [(pre_stim-abs(current_pre_stim)) pre_stim+current_post_stim]*fs;
        x_tick_vec = [0 0.5 current_post_stim];
    end
end

% look at retrieval from 0-2 sec
if strcmp('yes',lock_ext)&& strcmp('onset', lock)
    time_range(2)=(pre_stim+2)*fs;
    current_pre_stim =-.2;
    current_post_stim= 2;
    x_tick_vec = [0 0.5 1 1.5 current_post_stim];
end
%
[cond1a,cond2a,cond3a,cond4a] = get_power_subj_elecs(subj_list,exp_type,reg,ref,baseline,lock);
cond2 = cat(3,cond2a{:});
cond3 = cat(3,cond3a{:});

% reset only for onset
if reset==1
if strcmp('onset',lock)
resetVec = nanmean(nanmean(cond2(:,0.3*fs:pre_stim*fs,:),3),2);
cond2 = cond2 - resetVec;
resetVec = nanmean(nanmean(cond3(:,0.3*fs:pre_stim*fs,:),3),2);
cond3 = cond3 - resetVec;
end
end
% get the time ranges over which CBPT will be run
cond2=cond2(:,time_range(1):time_range(2),:);
cond3=cond3(:,time_range(1):time_range(2),:);

% remove any page that is all 0
for iElec = 1:size(cond2,3)
    if sum(all(cond2(:,:,iElec), 1)==0)>0
        cond2(:,:,iElec) =nan;
        disp('detect 0 chan - cond2')
    end

end
for iElec = 1:size(cond3,3)
    if sum(all(cond3(:,:,iElec), 1)==0)>0
        cond3(:,:,iElec) =nan;
        disp('detect 0 chan - cond2')
    end
end


both_conds = cat(3,cond3, cond2);

% plot to preview
mx = 4
mn = -mx

tickmarks = 1:18:length(freq);
figure;
subplot(2,1,1);imagesc(linspace(current_pre_stim,current_post_stim, size(cond2,2)), 1:length(freq),nanmean(cond2,3)./nanstd(cond2,0,3));hold on; colormap jet; title('lure --> -')
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Arial');colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
caxis([mn mx])
subplot(2,1,2);imagesc(linspace(current_pre_stim,current_post_stim, size(cond3,2)), 1:length(freq),nanmean(cond3,3)./nanstd(cond3,0,3));hold on; colormap jet; title('lure --> +')
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),'FontSize', 11, 'FontWeight', 'bold', 'FontName', 'Arial');colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
caxis([mn mx])
suptitle(reg)

% smoothen average gamma power for each subject
if strcmp('yes',smoothgamma)
    win= round(.135*fs)
    for subjn = 1:size(both_conds,3)
        for gamma_freq = find(freq>50)
            both_conds(gamma_freq,:,subjn) = conv(both_conds(gamma_freq,:,subjn), ones(1,win)/win,'same');
        end
    end
end
real_condition_mapping = [-ones(1,size(cond3,3)) ones(1,size(cond2,3))];
num_frex               = length(freq);
nTimepoints            = size(both_conds,2);
voxel_pval             = 0.05;
mcc_voxel_pval         = 0.05;
mcc_cluster_pval       = 0.05;
n_permutes             = 1000;

% compute actual t-test of difference (using unequal N and std)
tnum   = squeeze(nanmean(both_conds(:,:,real_condition_mapping==-1),3) - nanmean(both_conds(:,:,real_condition_mapping==1),3));
tdenom = sqrt((nanstd(both_conds(:,:,real_condition_mapping==-1),0,3).^2)./sum(real_condition_mapping==-1) + (nanstd(both_conds(:,:,real_condition_mapping==1),0,3).^2)./sum(real_condition_mapping==1));
real_t = tnum./tdenom;

% initialize null hypothesis matrices
permuted_tvals  = zeros(n_permutes,num_frex,nTimepoints);
max_pixel_pvals = zeros(n_permutes,2);
max_clust_info  = zeros(n_permutes,1);
max_clust_t_stat= zeros(n_permutes,2);

% generate pixel-specific null hypothesis parameter distributions
for permi = 1:n_permutes
    
    % shuffle labels
    fake_condition_mapping = real_condition_mapping;
    shuf_cntr = [];
    for ii = 1:size(cond3,3)
        shuf_cntr = [shuf_cntr randi([0 1],1)];
    end
    shuffle = logical([shuf_cntr  shuf_cntr]);
    fake_condition_mapping(shuffle) = -1*real_condition_mapping(shuffle);
    
    % compute t-map of null hypothesis
    tnum   = squeeze(nanmean(both_conds(:,:,fake_condition_mapping==-1),3)-nanmean(both_conds(:,:,fake_condition_mapping==1),3));
    tdenom = sqrt((nanstd(both_conds(:,:,fake_condition_mapping==-1),0,3).^2)./sum(fake_condition_mapping==-1) + (nanstd(both_conds(:,:,fake_condition_mapping==1),0,3).^2)./sum(fake_condition_mapping==1) );
    tmap   = tnum./tdenom;
    
    % save all permuted values
    permuted_tvals(permi,:,:) = tmap;
    
    % save maximum pixel values
    max_pixel_pvals(permi,:) = [min(tmap(:)) max(tmap(:))];
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % clusters obtained by parametrically thresholding the t-maps
    tmap(abs(tmap)<tinv(1-voxel_pval,size(both_conds,3)-2))=0;
    
    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(tmap);
    clust_size = [];
    for iClust = 1:size(clustinfo.PixelIdxList,2)
        clust_size = [clust_size size(clustinfo.PixelIdxList{iClust},1)];
    end
    max_clust_info(permi) = max([0 clust_size]); % zero accounts for empty maps
    
    % loop thru clusters, find cluster w/ max t-statistic and max
    % weighted sum
    clust_sum = [];
    clust_weighted_sum = [];
    for clus_cntr = 1:length(clustinfo.PixelIdxList)
        clust_sum = [clust_sum sum(tmap(clustinfo.PixelIdxList{clus_cntr}))];
        clust_weighted_sum =[clust_weighted_sum sum(tmap(clustinfo.PixelIdxList{clus_cntr}))/size(clustinfo.PixelIdxList{clus_cntr},1)];
        
    end
    max_clust_t_stat(permi,:) = [min(clust_sum) max(clust_sum)];
    max_clust_t_stat_weighted(permi,:) = [min(clust_weighted_sum) max(clust_weighted_sum)];
    disp(permi)
end

% plot
tickmarks = 1:20:length(freq);

% now compute Z-map
zmap = (real_t-squeeze(nanmean(permuted_tvals,1)))./squeeze(nanstd(permuted_tvals));

% apply uncorrected threshold
figure
hold on
mn = -4
mx =  4
if strcmp('encoding', exp_type)
    contourf(linspace(current_pre_stim,current_post_stim, size(zmap,2)), 1:length(freq),zmap,40,'linecolor','none')
else
    contourf(linspace(current_pre_stim,current_post_stim, size(zmap,2)), 1:length(freq),zmap,40,'linecolor','none')
end
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar

zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=false; % going from Z value to p value % norminv(1-voxel_pval) is the standard def unit that corrsep to p<.05
zmapthresh=logical(zmapthresh);
hold on
if strcmp('encoding', exp_type)
    contour(linspace(current_pre_stim,current_post_stim, size(zmap,2)), 1:length(freq),zmapthresh,3,'linecolor','k')
else
    contour(linspace(current_pre_stim,current_post_stim, size(zmap,2)), 1:length(freq),zmapthresh,3,'linecolor','k')
end
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
colormap jet
caxis([mn mx])
xticks([x_tick_vec])
colorbar('Ticks',[mn mx],...
  'TickLabels',[mn mx])
%title(reg)
xlabel('Time (s)'), ylabel('Frequency (Hz)')
print('-clipboard','-dbitmap')
%% get p-map play
p_twoSided = (1-normcdf(abs(zmap)))*2;

figure;
imagesc(linspace(current_pre_stim,current_post_stim, size(zmap,2)), 1:length(freq),...
    -1*p_twoSided)
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
hold on
xlim([current_pre_stim current_post_stim])
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
colormap jet
caxis([ -.05 0.0001 ])
colorbar
xticks([x_tick_vec])
xlabel('Time (s)'), ylabel('Frequency (Hz)')
print('-clipboard','-dbitmap')

%% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
clustinfo  = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);

if strcmp('maxsize',clustercorect)
    % % MC based on size
    clust_size_threshold   = prctile(max_clust_info,100-mcc_cluster_pval*100);
    % identify clusters to remove
    whichclusters2remove   = find(clust_info<clust_size_threshold);
    
elseif strcmp('ws',clustercorect)
    %% based on weightsum
    clust_WSUM_upper_threshold      = prctile(max_clust_t_stat_weighted(:,2),100-mcc_cluster_pval*100);
    clust_WSUM_lower_threshold = prctile(max_clust_t_stat_weighted(:,1),mcc_cluster_pval*100);
    clust_WSUM=[];
    for clus_cntr = 1:length(clust_info)
        clust_WSUM(clus_cntr) = sum(zmapthresh(clustinfo.PixelIdxList{clus_cntr}))/size(clustinfo.PixelIdxList{clus_cntr},1);
    end
    sig_pos_clust = find(clust_WSUM>clust_WSUM_upper_threshold);
    sig_neg_clust = find(clust_WSUM<clust_WSUM_lower_threshold);
    sig_clusters = [sig_pos_clust  sig_neg_clust];
    non_sig_clusters = ones(1,length(clust_info));
    non_sig_clusters(sig_clusters) = 0;
    whichclusters2remove   = find(non_sig_clusters);
    
elseif strcmp('maxsum',clustercorect)   
    % MC based on sum of stat
    clust_t_stat_upper_threshold = prctile(max_clust_t_stat(:,2),100-mcc_cluster_pval*100);
    clust_t_stat_lower_threshold = prctile(max_clust_t_stat(:,1),mcc_cluster_pval*100);
    clust_t_sum=[];
    for clus_cntr = 1:length(clust_info)
        clust_t_sum(clus_cntr,:) = sum(zmapthresh(clustinfo.PixelIdxList{clus_cntr}));
    end
    sig_pos_clust = find(clust_t_sum>clust_t_stat_upper_threshold);
    sig_neg_clust = find(clust_t_sum<clust_t_stat_lower_threshold);
    sig_clusters = [sig_pos_clust ; sig_neg_clust];
    
    non_sig_clusters = ones(1,length(clust_info));
    non_sig_clusters(sig_clusters) = 0;
    whichclusters2remove   = find(non_sig_clusters);
end

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

mn = -4
mx = 4
figure
contourf(linspace(current_pre_stim,current_post_stim, size(zmapthresh,2)), 1:length(freq), zmapthresh,40,'linecolor','none')
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
xlabel('Time (s)'), ylabel('Frequency (Hz)')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
caxis([mn mx])
colorbar('Ticks',[mn 0 mx])
colorbar('Ticks',[mn 0 mx],...
  'TickLabels',[mn 0 mx])
xticks([-1 0  current_post_stim ])
print('-clipboard','-dbitmap')
title('Cluster-corrected Z map')
%suptitle(['grp - region: ' reg ])


%% run once
zmapthresh_temp = zmapthresh;

%% refresh
zmapgamma = zmapthresh_temp;

%% encoding onset
zmapgamma (freq>7,:)=0;
zmapgamma (:,1:-1*current_pre_stim*fs)=0;

%GAMMA
zmapgamma (freq<55,:)=0;
% CA1
%first low freq
zmapgamma (:,1:.7*fs)=0;
%first low freq
zmapgamma (:,.72*fs:end)=0;

% CA3
%DELTA

zmapgamma (freq>3.9,:)=0;
%THETA
zmapgamma (freq<3.9,:)=0;
zmapgamma (freq>6.5,:)=0;
zmapgamma (freq>5.5,.9*fs:1.1*fs)=0;
zmapgamma (freq>6,1.1*fs:1.2*fs)=0;
% ALPHA
zmapgamma (freq<6.4,:)=0;

% HC
%delta
zmapgamma (freq>5,:)=0;
zmapgamma (:,1:.45*fs)=0;
zmapgamma (freq>3.2,.45*fs:.55*fs)=0;
zmapgamma (freq>4.6,1.12*fs:1.6*fs)=0;
zmapgamma (freq>4.3,1.3*fs:1.6*fs)=0;

% alpha
zmapgamma (:,1:.75*fs)=0;
zmapgamma (freq<5.5,:)=0;
zmapgamma (freq>12,.75*fs:.85*fs)=0;
zmapgamma (freq<6.5,1.4*fs:1.72*fs)=0;
% alpha

% FRO
% gamma
zmapgamma (freq<29,:)=0;
zmapgamma (freq<44,1:1*fs)=0;
% alpha
zmapgamma (freq>44,:)=0;
zmapgamma (freq>29,1.1*fs:end)=0;

%%
zmapgamma (freq>4,:)=0;
zmapgamma (:,1:0.8*fs)=0;

%%
mn = -4
mx = 4
tickmarks = 1:20:length(freq);
figure
imagesc(linspace(current_pre_stim,current_post_stim, size(zmap,2)),1:length(freq), zmapgamma)
hold on
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
xlabel('Time (s)'), ylabel('Frequency (Hz)')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
colormap jet
title(reg)
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
caxis([mn mx])
colorbar('Ticks',[mn mx],...
  'TickLabels',[mn mx])
xticks([0 1])
print('-clipboard','-dbitmap')
%% retrieval RESPONSE

% CA1
%low freq
zmapgamma (freq>6,:)=0;
zmapgamma (freq<5,:)=0;

%GAMMA
zmapgamma (freq<13,:)=0;

% HC
%low freq
zmapgamma (freq>15,:)=0;



% OFC
%low freq
zmapgamma (freq>17,:)=0;
%GAMMA
zmapgamma (freq<17,:)=0;

% TEMP
%low freq
zmapgamma (freq>8,:)=0;
zmapgamma (freq>5,1:.5*fs)=0;

% CING
%low freq
zmapgamma (freq>7,:)=0;
% gamma
zmapgamma (freq<35,:)=0;

% EC
% gamma
zmapgamma (freq<35,:)=0;
% theta
zmapgamma (freq>35,:)=0;


%% retrieval onset
% NC
%low freq
zmapgamma (freq>7,:)=0;
zmapgamma (:,1:.2*fs)=0;

%delta
zmapgamma (freq>4,:)=0;
%theta
zmapgamma (freq<4,:)=0;
zmapgamma (freq>5.6,:)=0;

%GAMMA hi
zmapgamma (:,1:.3*fs)=0;
zmapgamma (freq<19,:)=0;
zmapgamma (:,1:.88*fs)=0;

%GAMMA lo
zmapgamma (:,1:.3*fs)=0;
zmapgamma (freq<19,:)=0;
zmapgamma (:,.89*fs:end)=0;

% CA1
%low freq
zmapgamma (:,.85*fs:end)=0;
zmapgamma (freq>6.5,:)=0;
%hi freq
zmapgamma (:,1:.5*fs)=0;
zmapgamma (freq<6.5,:)=0;

% HC
%low freq
zmapgamma (freq>6.3,:)=0;
zmapgamma (freq>5.5,.88*fs:end)=0;
%hi freq
zmapgamma (freq<29,:)=0;
zmapgamma (:,.1:.62*fs)=0;
zmapgamma (freq<34,.75*fs:.88*fs)=0;

% CA3
%low freq
zmapamma (freq>5,:)=0;

%hi freq
zmapgamma (freq<15,:)=0;

%OFC
%delta
zmapgamma (freq>5.7,:)=0;
%theta
zmapgamma (freq<5.7,:)=0;

%EC
%delta
zmapgamma (:,1:.5*fs)=0;

%FRO
%deltatheta
zmapgamma (freq>10,:)=0;
zmapgamma (freq>7, 1:.4*fs)=0;
%delta
zmapgamma (freq>4.2,1:.5*fs)=0;
zmapgamma (freq>5, .5*fs:end)=0;
zmapgamma (freq>4.4, .48*fs:.7*fs)=0;
zmapgamma (freq>4, .53*fs:.58*fs)=0;
zmapgamma (freq>4.3, .7*fs:.75*fs)=0;
%theta
zmapgamma (freq<4.2,1:.5*fs)=0;
zmapgamma (freq<5, .5*fs:end)=0;
zmapgamma (freq<4.4, .48*fs:.7*fs)=0;
zmapgamma (freq<4, .53*fs:.58*fs)=0;
zmapgamma (freq<4.3, .7*fs:.75*fs)=0;
%gamma
zmapgamma (freq<17,:)=0;
zmapgamma (:,1:.2*fs)=0;
zmapgamma (freq<20,1:.2*fs)=0;
zmapgamma (freq<24,.25:.34*fs)=0;


%temp
%gamma
zmapgamma (freq<23,:)=0;


%% save specific cluster
zmapthresh=zmapgamma;
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/cluster_matrices'])
%%
save([reg '_' lock '_' exp_type '_groupcluster_ALL'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_deltatheta'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_delta'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_theta1'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_theta2'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_theta'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_alpha'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_beta'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_gamma'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_gammaHI'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_gammaLO'],'zmapthresh')
save([reg '_' lock '_' exp_type '_groupcluster_deactive_deltatheta'],'zmapthresh')
