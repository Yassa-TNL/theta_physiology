% average across time after getting single value for all subjs and chans
figure;bar([ nanmean(conda_val_time_mn)-bidirec nanmean(condb_val_time_mn)-bidirec])

% average across time for each single value coming from subjs X chans
conda_collapse_time = nanmean(conda_val_time_mn_store,1)
condb_collapse_time = nanmean(condb_val_time_mn_store,1)
conda_collapse_grp = nanmean(conda_collapse_time)
condb_collapse_grp = nanmean(condb_collapse_time)

figure;bar([ conda_collapse_grp condb_collapse_grp]-bidirec)
ylim([-0.01 .004])