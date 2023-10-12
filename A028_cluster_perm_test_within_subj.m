% Arrange data for statistical testing
close all;clc
regList       = {'INS'}  % HC CA3 CA1 OFC FRO TEMP CING INS EC
acrsTrls      = 'no';
clustercorect = 'maxsum'; %maxsum maxsize ws


addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
[fro_chan_idx,MTL_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,OFC_chan_idx,CA3_chan_idx,CA1_chan_idx,HC_chan_idx] = get_elecs(subj);
time_range = [501 1500]; %0ms before stim onset to 1300ms after

for counter = 1:length(regList)
clear both_conds permuted_tvals max_clust_t_stat clust_sum zmapthresh zmap cond1 cond2 temp1 temp2
reg = regList{counter}
if strcmp('OFC',reg)
    chan_idx = OFC_chan_idx;
elseif strcmp('FRO',reg)
    chan_idx = fro_chan_idx;
elseif strcmp('TEMP',reg)
    chan_idx = temp_chan_idx;
elseif strcmp('CING',reg)
    chan_idx = cingulate_chan_idx;
elseif strcmp('INS',reg)
    chan_idx = insula_chan_idx;
elseif strcmp('MTL',reg)
    chan_idx = MTL_chan_idx;
elseif strcmp('HC',reg)
    chan_idx = HC_chan_idx;
elseif strcmp('CA1',reg)
    chan_idx = CA1_chan_idx;
elseif strcmp('CA3',reg) 
    chan_idx = CA3_chan_idx;
end

    
% get data
 if strcmp('no',acrsTrls)
    cond1 = squeeze(nanmean(norm_freq_acrs_chan_cond_3(:, time_range(1):time_range(2),:, chan_idx),3));
    cond2 = squeeze(nanmean(norm_freq_acrs_chan_cond_2(:, time_range(1):time_range(2),:, chan_idx),3));
 else
    temp1 = norm_freq_acrs_chan_cond_3(:, time_range(1):time_range(2),:, chan_idx);
    temp2 = norm_freq_acrs_chan_cond_2(:, time_range(1):time_range(2),:, chan_idx);
    
    cond1 = nan(size(temp1,1), size(temp1,2), size(temp1,3)*length(chan_idx));
    cond2 = nan(size(temp2,1), size(temp2,2), size(temp2,3)*length(chan_idx));
    
    for chan = 1:length(chan_idx)
        if chan ==1
            cond1 (:,:, 1:size(temp1,3)) = temp1(:,:,:,chan);
            cond2 (:,:, 1:size(temp2,3)) = temp2(:,:,:,chan);
        else
            cond1 (:,:, (size(temp1,3)*(chan-1))+1:size(temp1,3)*chan)= temp1(:,:,:,chan);
            cond2 (:,:, (size(temp2,3)*(chan-1))+1:size(temp2,3)*chan)= temp2(:,:,:,chan);
        end
    end
end
both_conds = cat(3,cond1, cond2);

% smoothen average gamma power for each subject
win= .135*fs;
for subjn = 1:size(both_conds,3)
    for gamma_freq = find(freq>40)
        both_conds(gamma_freq,:,subjn) = conv(both_conds(gamma_freq,:,subjn), ones(1,win)/win,'same');
    end
end
%%
real_condition_mapping = [-ones(1,size(cond1,3)) ones(1,size(cond2,3))];
num_frex               = length(freq);
nTimepoints            = size(both_conds,2);
voxel_pval             = 0.05;
mcc_cluster_pval       = 0.05;
n_permutes             = min(nchoosek(length(chan_idx)*2, length(chan_idx)),1000);

% compute actual t-test of difference (using unequal N and std)
tnum   = squeeze(nanmean(both_conds(:,:,real_condition_mapping==-1),3) - nanmean(both_conds(:,:,real_condition_mapping==1),3));
tdenom = sqrt((nanstd(both_conds(:,:,real_condition_mapping==-1),0,3).^2)./sum(real_condition_mapping==-1) + (nanstd(both_conds(:,:,real_condition_mapping==1),0,3).^2)./sum(real_condition_mapping==1));
real_t = tnum./tdenom;

% initialize null hypothesis matrices
permuted_tvals  = zeros(n_permutes,num_frex,nTimepoints);
max_pixel_pvals = zeros(n_permutes,2);
max_clust_info  = zeros(n_permutes,1);
max_clust_t_stat= zeros(n_permutes,2);
max_clust_t_stat_weighted = zeros(n_permutes,2);

if length(chan_idx)>0
    
    % generate pixel-specific null hypothesis parameter distributions
    for permi = 1:n_permutes
        fake_cond_idx = randperm(length(real_condition_mapping));
        fake_condition_mapping = real_condition_mapping(fake_cond_idx);
        % compute t-map of null hypothesis
        tnum   = squeeze(nanmean(both_conds(:,:,fake_condition_mapping==-1),3)-nanmean(both_conds(:,:,fake_condition_mapping==1),3));
        tdenom = sqrt((nanstd(both_conds(:,:,fake_condition_mapping==-1),0,3).^2)./sum(fake_condition_mapping==-1) + (nanstd(both_conds(:,:,fake_condition_mapping==1),0,3).^2)./sum(fake_condition_mapping==1) );
        tmap   = tnum./tdenom;
        
        % save all permuted values
        permuted_tvals(permi,:,:) = tmap;
        
        % save maximum pixel values
        max_pixel_pvals(permi,:) = [min(tmap(:)) max(tmap(:))];
        
        % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
        % note that here, clusters were obtained by parametrically thresholding
        % the t-maps
        tmap(abs(tmap)<tinv(1-voxel_pval,size(both_conds,3)-2))=0;
        
        % get number of elements in largest supra-threshold cluster
        clustinfo = bwconncomp(tmap);
        max_clust_info(permi) = max([0 cellfun(@numel,clustinfo.PixelIdxList)]); % zero accounts for empty maps
        
        % loop thru clusters, find cluster w/ max t-statistic and max
        % weighted sum
        clust_sum = [];
        clust_weighted_sum = [];
        for clus_cntr = 1:length(clustinfo.PixelIdxList)
            clust_sum(clus_cntr) = sum(tmap(clustinfo.PixelIdxList{clus_cntr}));
            clust_weighted_sum(clus_cntr) = sum(tmap(clustinfo.PixelIdxList{clus_cntr}))/size(clustinfo.PixelIdxList{clus_cntr},1);

        end
        max_clust_t_stat(permi,:) = [min(clust_sum) max(clust_sum)];
        max_clust_t_stat_weighted(permi,:) = [min(clust_weighted_sum) max(clust_weighted_sum)];

        disp(permi)
    end
end

%%
% now compute Z-map
zmap = (real_t-squeeze(nanmean(permuted_tvals,1)))./squeeze(nanstd(permuted_tvals));

% computer a p matrix making no assumption about distribution
p_mtx = nan(size(real_t));
for row = 1:size(real_t,1)
    for column = 1:size(real_t,2)
        time_freq_null   = permuted_tvals(:,row,column);
        if real_t(row,column)>0
            p_mtx(row,column) =  sum(time_freq_null>real_t(row,column))/n_permutes; % right tail
        elseif  real_t(row,column)<0
            p_mtx(row,column) =  sum(time_freq_null<real_t(row,column))/n_permutes; % left tail
        end
        clear time_freq_null
    end
end

% threshold the p matrix
p_mtx(p_mtx>0.05 ) = nan;
zmapthresh=p_mtx; 


%hold on
%contour(linspace(0,1,size(both_conds,2)), 1:length(freq),zmapthresh,1,'linecolor','k')

% find connected voxels
clustinfo  = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);

if strcmp('maxsize', clustercorect)
    %% MC based on maxsize
    clust_size_threshold   = prctile(max_clust_info,100-mcc_cluster_pval*100);
    % identify clusters to remove
    whichclusters2remove   = find(clust_info<clust_size_threshold);
    % remove clusters
    for i=1:length(whichclusters2remove)
        zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end
    
elseif strcmp('ws', clustercorect)
    %% based on weightsum
    clust_WSUM_upper_threshold      = prctile(max_clust_t_stat_weighted(:,2),100-mcc_cluster_pval*100);
    clust_WSUM_lower_threshold = prctile(max_clust_t_stat_weighted(:,1),mcc_cluster_pval*100);
    clust_WSUM=[];
    for clus_cntr = 1:length(clust_info)
        clust_WSUM(clus_cntr) = sum(zmapthresh(clustinfo.PixelIdxList{clus_cntr}))/size(clustinfo.PixelIdxList{clus_cntr},1);
    end
    sig_pos_clust = find(clust_WSUM>clust_WSUM_upper_threshold);
    sig_neg_clust = find(clust_WSUM<clust_WSUM_lower_threshold);
    sig_clusters = [sig_pos_clust ; sig_neg_clust];
    non_sig_clusters = ones(1,length(clust_info));
    non_sig_clusters(sig_clusters) = 0;
    whichclusters2remove   = find(non_sig_clusters);
    % remove clusters
    for i=1:length(whichclusters2remove)
        zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end
    
elseif strcmp('maxsum', clustercorect)
    
    %% MC based on sum of maxstat
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
    % remove clusters
    for i=1:length(whichclusters2remove)
        zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
    end
    
end

if ~any(any(zmapthresh))
   % plot
mx = 0.05
mn = 0
figure
tickmarks = 1:20:length(freq);
contourf(linspace(0,1,size(both_conds,2)), 1:length(freq),p_mtx,40,'linecolor','none')
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),'clim',[mn mx])
axis square
colorbar
axis square
colorbar
title('Cluster-corrected Z map')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
suptitle(['subj:' num2str(subj) ' - region:' reg])
else
    
    % plot
mx = 0.05
mn = 0
figure
tickmarks = 1:20:length(freq);
subplot(211)
contourf(linspace(0,1,size(both_conds,2)), 1:length(freq),p_mtx,40,'linecolor','none')
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),'clim',[mn mx])
axis square
colorbar
zmapthresh=p_mtx;
% plot after correc for multiple comparisons
subplot(212)
contourf(linspace(0,1,size(both_conds,2)), 1:length(freq),zmapthresh,40,'linecolor','none')
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),'clim',[mn mx])
axis square
colorbar
title('Cluster-corrected Z map')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
suptitle(['subj:' num2str(subj) ' - region:' reg])
end
% save in subj folder
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/single_subj_perm_test'])
save(['perm_test_single_subj_tmapz_sig_cluster_' reg], 'zmapthresh', 'zmap', 'permuted_tvals')

% save in group folder
cd(['/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_' lock '_' exp_type])
save(['subj_' subj '_reg_' reg 'zmapthresh'],'zmapthresh')
end

%%
% a=5.7253
% line([0 1], [174 174])
% line([0 1], [max(find(round(freq)==a)) max(find(round(freq)==a))])
% zmapthresh (freq>120,180:200)=0;