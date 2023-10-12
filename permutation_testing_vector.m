function [zmap,zmapthresh,zmapthresh_for_plot] = permutation_testing_vector(cond3,cond2,n_permutes)
% Runs perm testing on vectors rather than matrices, and channels being the
% observations (first dim). 
% Inputs: two conditions, observationsXtime
% Output: zmap, zmapthresh, and thresh for plot
both_conds = cat(1,cond3, cond2);
real_condition_mapping = [-ones(1,size(cond3,1)) ones(1,size(cond2,1))];
voxel_pval             = 0.05;
mcc_voxel_pval         = 0.05;
mcc_cluster_pval       = 0.05;

% compute actual t-test of difference (using unequal N and std)
tnum   = squeeze(nanmean(both_conds(real_condition_mapping==-1,:),1) - nanmean(both_conds(real_condition_mapping==1,:),1));
tdenom = sqrt((nanstd(both_conds(real_condition_mapping==-1,:),0,1).^2)./sum(real_condition_mapping==-1) + (nanstd(both_conds(real_condition_mapping==1,:),0,1).^2)./sum(real_condition_mapping==1));
real_t = tnum./tdenom;

% initialize null hypothesis matrices
permuted_tvals  = zeros(n_permutes,length(real_t));
max_pixel_pvals = zeros(n_permutes,2);
max_clust_info  = zeros(n_permutes,1);
max_clust_t_stat= zeros(n_permutes,2);

%% generate pixel-specific null hypothesis parameter distributions
for permi = 1:n_permutes
    % perm
     numtot = length(real_condition_mapping); rand_idx=randperm(numtot);
idx_cond1 = rand_idx(1:size(cond3,1)); idx_cond2 =  rand_idx(size(cond3,1)+1:end);
   %fake_condition_mapping = sign(randn(length(real_condition_mapping),1));

    % compute t-map of null hypothesis
    tnum   = squeeze(nanmean(both_conds(idx_cond1,:),1) - nanmean(both_conds(idx_cond2,:),1));
    tdenom = sqrt((nanstd(both_conds(idx_cond1,:),0,1).^2)./length(idx_cond1) + (nanstd(both_conds(idx_cond2,:),0,1).^2)./length(idx_cond2));
    tmap   = tnum./tdenom;     
    % save all permuted values
    permuted_tvals(permi,:) = tmap;
    
    % save maximum pixel values
    max_pixel_pvals(permi,:) = [min(tmap(:)) max(tmap(:))];
    
    % for cluster correction, apply uncorrected threshold and get maximum cluster sizes
    % note that here, clusters were obtained by parametrically thresholding
    % the t-maps
    tmap(abs(tmap)<tinv(1-voxel_pval,size(both_conds,1)-2))=0;

    % get number of elements in largest supra-threshold cluster
    clustinfo = bwconncomp(tmap);
    max_clust_info(permi) = max([0 cellfun(@numel,clustinfo.PixelIdxList)]); % zero accounts for empty maps
    
    % loop thru clusters, find cluster w/ max t-statistic and max
        % weighted sum
        clust_sum = [];
        for clus_cntr = 1:length(clustinfo.PixelIdxList)
            clust_sum(clus_cntr) = sum(tmap(clustinfo.PixelIdxList{clus_cntr}));
        end
        if ~isempty([min(clust_sum) max(clust_sum)]) 
        max_clust_t_stat(permi,:) = [min(clust_sum) max(clust_sum)];
        end

end

%% now compute Z-map
zmap = (real_t-nanmean(permuted_tvals,1))./squeeze(nanstd(permuted_tvals,0,1));
zmapthresh = zmap;

zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=false; % going from Z value to p value % norminv(1-voxel_pval) is the standard def unit that corrsep to p<.05


%%
clustinfo  = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);
%%
% % MC based on size
clust_size_threshold   = prctile(max_clust_info,100-mcc_cluster_pval*100);
% identify clusters to remove
whichclusters2remove   = find(clust_info<clust_size_threshold);
    
    %% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end

%% nans for non sig voxels, 1 for sig voxels 
zmapthresh=logical(zmapthresh);
zmapthresh_for_plot = nan(size(zmapthresh));
zmapthresh_for_plot(zmapthresh~=0) = 1;
end

