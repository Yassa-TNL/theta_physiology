cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/figures/cfc'])
% for a given region pair
perm_distrib = squeeze(nanmean(Perm_PAC, 4)); % average null distribution from all elec pairs
p            = 95;
per_val_map  = prctile(perm_distrib,p,3);
per_val_map_z = (per_val_map - nanmean(perm_distrib,3))./nanstd(perm_distrib,0,3);
%% compare observed z score to significant z score
pacz_mean    = nanmean(pacz,3);
pacz_mean(find(pacz_mean<per_val_map_z)) = 0;
imagesc(pacz_mean)



%% individual subj uncorr stats

for reg1 = 1:length(phase_regions_list) %  loop thru pairs of regions
        phase_region = phase_regions_list{reg1};
        
        for reg2 = 1:length(power_regions_list)
            power_region = power_regions_list{reg2};
                   % go to patient dir
            cd(['/mnt/yassamri/iEEG/sandra/subj_' group{ptnt} '/figures/cfc'])
            

        end
end
