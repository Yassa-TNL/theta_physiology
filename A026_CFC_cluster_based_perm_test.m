phase_regions_list = {'CING'} 
power_regions_list = {'HC'}  

    % cortex w/ MTL
if strcmp('OFC', phase_regions_list) | strcmp('FRO', phase_regions_list) | strcmp('TEMP', phase_regions_list) && strcmp('MTL', power_regions_list)
    group = {'39' '57'  '66'  '44' '63' '84'}

elseif strcmp('CING', phase_regions_list) && strcmp('MTL', power_regions_list)
    group = {'39' '57'  '66'  '44' '63'}
    
    % cortex w/ HC and cortex w/ CA3
elseif strcmp('OFC', phase_regions_list) | strcmp('FRO', phase_regions_list) | strcmp('TEMP', phase_regions_list) | strcmp('CING', phase_regions_list) && strcmp('HC', power_regions_list) | strcmp('CA3', power_regions_list)
    group = {'39' '57'  '66'  '44' }
end

cond2a      = cell(1,length(group));
cond3a      = cell(1,length(group));

cond2a_null = cell(1,length(group));
cond3a_null = cell(1,length(group));
for reg1 = 1:length(phase_regions_list) %  loop thru pairs of regions
    phase_region = phase_regions_list{reg1};
    for reg2 = 1:length(power_regions_list)
        power_region = power_regions_list{reg2};
        for ptnt= 1:length(group) % loop thru ptnts
            % go to patient dir
            cd(['/mnt/yassamri/iEEG/sandra/subj_' group{ptnt} '/figures/cfc'])

            % cond2
            load(['pacz_' phase_region '_' power_region '_cond2.mat'])
            cond2a{ptnt} = pacz;

            % cond3
            load(['pacz_' phase_region '_' power_region '_cond3.mat'])
            cond3a{ptnt} = pacz;
        end
    end
end
cond2 = cat(3,cond2a{:});
cond3 = cat(3,cond3a{:});

% 1 vector, average across gamma
cond2 = nanmean(cat(3,cond2a{:}),1);
cond3 = nanmean(cat(3,cond3a{:}),1);

% Runs perm testing on vectors rather than matrices, and channels being the
% observations. 
% Inputs: two conditions that are XxYxZ, X = 1. Y = vector dimension. Z =
% elec or observation dimension
% Output: zmap
both_conds = cat(3,cond3, cond2);
real_condition_mapping = [-ones(1,size(cond3,3)) ones(1,size(cond2,3))];
voxel_pval             = 0.02;
mcc_voxel_pval         = 0.05;
mcc_cluster_pval       = 0.05;
n_permutes             = 1000;
num_frex_pow           = size(cond3,1);
num_frex_phas          = size(cond3,2);
% compute actual t-test of difference (using unequal N and std)
tnum   = squeeze(nanmean(both_conds(:,:,real_condition_mapping==-1),3) - nanmean(both_conds(:,:,real_condition_mapping==1),3));
tdenom = sqrt((nanstd(both_conds(:,:,real_condition_mapping==-1),0,3).^2)./sum(real_condition_mapping==-1) + (nanstd(both_conds(:,:,real_condition_mapping==1),0,3).^2)./sum(real_condition_mapping==1));
real_t = tnum./tdenom;

% initialize null hypothesis matrices
permuted_tvals  = zeros(n_permutes,num_frex_pow,num_frex_phas);
max_pixel_pvals = zeros(n_permutes,2);
max_clust_info  = zeros(n_permutes,1);
max_clust_t_stat= zeros(n_permutes,2);

% generate pixel-specific null hypothesis parameter distributions
for permi = 1:n_permutes
    % choose how many subj to perm
   fake_condition_mapping = sign(randn(length(real_condition_mapping),1));
  
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
        for clus_cntr = 1:length(clustinfo.PixelIdxList)
            clust_sum(clus_cntr) = sum(tmap(clustinfo.PixelIdxList{clus_cntr}));
        end
        if ~isempty([min(clust_sum) max(clust_sum)]) 
        max_clust_t_stat(permi,:) = [min(clust_sum) max(clust_sum)];
        end
    disp(permi)
end

% now compute Z-map
zmap = (real_t-squeeze(nanmean(permuted_tvals,1)))./squeeze(nanstd(permuted_tvals));

% now compute Z-map
power_freqs = 2.^(5.3:0.15:8);
phas_freqs  = 2.^(1.56:0.3:5.15); 
x=floor(phas_freqs*10)/10;
y=floor(power_freqs*10)/10

% apply uncorrected threshold
figure
subplot(211)
contourf(zmap,40,'linecolor','none')
set(gca,'xtick', 1:2:length(phas_freqs),'xticklabel', x(1:2:length(phas_freqs)),'YDir','normal', 'ytick', 1:2:length(power_freqs),'yticklabel', y(1:2:length(y)),'FontSize', 14, 'FontWeight', 'bold');
colorbar
%caxis([mn mx])
colormap jet
zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=false; % going from Z value to p value % norminv(1-voxel_pval) is the standard def unit that corrsep to p<.05
zmapthresh=logical(zmapthresh);
hold on
contour(zmapthresh,1,'linecolor','k')
title('Unthresholded Z map')
xlabel(['Freq for Phase ']);
ylabel(['Freq for Power ' ]);
set(gca, 'FontSize', 12, 'FontWeight', 'bold', 'FontName', 'Arial')

% apply cluster-level corrected threshold
zmapthresh = zmap;
% uncorrected pixel-level threshold
zmapthresh(abs(zmapthresh)<norminv(1-voxel_pval))=0;
clustinfo  = bwconncomp(zmapthresh);
clust_info = cellfun(@numel,clustinfo.PixelIdxList);

% % MC based on size
clust_size_threshold   = prctile(max_clust_info,100-mcc_cluster_pval*100);
% identify clusters to remove
whichclusters2remove   = find(clust_info<clust_size_threshold);

% remove clusters
for i=1:length(whichclusters2remove)
    zmapthresh(clustinfo.PixelIdxList{whichclusters2remove(i)})=0;
end 

subplot(212)
contourf(zmapthresh,40,'linecolor','none')
set(gca,'xticklabel', x(1:2:length(phas_freqs)),'YDir','normal', 'ytick', 1:2:length(power_freqs),'yticklabel', y(1:2:length(y)),'FontSize', 14, 'FontWeight', 'bold');
colorbar
title('Cluster-corrected Z map')
title('Unthresholded Z map')
xlabel(['Freq for Phase ']);
ylabel(['Freq for Power ' ]);
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
