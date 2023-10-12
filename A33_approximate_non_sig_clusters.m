load('subj_39_reg_MTLzmapthresh.mat')
zmap_all(:,:,1) = zmapthresh;
load('subj_44_reg_MTLzmapthresh.mat')
zmap_all(:,:,2) = zmapthresh;
load('subj_57_reg_MTLzmapthresh.mat')
zmap_all(:,:,3) = zmapthresh;
load('subj_63_reg_MTLzmapthresh.mat')
zmap_all(:,:,4) = zmapthresh;
load('subj_66_reg_MTLzmapthresh.mat')
avgmap = nanmean(zmap_all,3);
imagesc(avgmap./(max(max(avgmap))))
colorbar
clear zmapthresh
zmapthresh= avgmap./(max(max(avgmap)));
imagesc(zmapthresh)