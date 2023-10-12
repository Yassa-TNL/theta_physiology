group = {'39' '57' '44' '63' '66' '84'}%
phase_regions_list = {'OFC', 'FRO', 'TEMP'} %, 'CING'
power_regions_list = {'MTL'} % 'HC'
%%
group = {'39' '57' '44' '63' '66' }%'84'
phase_regions_list = {'CING'} %, 'CING'
power_regions_list = {'MTL'} % 'HC'

%%
group = {'39' '57' '44'  '66' }%'84''63'
phase_regions_list = {'OFC', 'FRO', 'TEMP' 'CING'} %, 'CING'
power_regions_list = {'HC'} % 'HC'

%% initializa group matrices for each pair of regions and each condition
grp_pacx_cond1 = zeros(size(pacz,1), size(pacz,2), length(group));
grp_pacx_cond2 = zeros(size(pacz,1), size(pacz,2), length(group));
grp_pacx_cond3 = zeros(size(pacz,1), size(pacz,2), length(group));
grp_pacx_cond4 = zeros(size(pacz,1), size(pacz,2), length(group));

rng1 = 120
rng2= 3.4
close all
    for reg1 = 1:length(phase_regions_list) %  loop thru pairs of regions
        phase_region = phase_regions_list{reg1};
        
        for reg2 = 1:length(power_regions_list)
            power_region = power_regions_list{reg2};
            
            for ptnt= 1:length(group) % loop thru ptnts
            % go to patient dir
            cd(['/mnt/yassamri/iEEG/sandra/subj_' group{ptnt} '/figures/cfc'])
            
            % save as group data matrix
            load(['pacz_' phase_region '_' power_region '_cond1.mat'])
            grp_pacx_cond1(:,:, ptnt) = nanmean(pacz,3);

            load(['pacz_' phase_region '_' power_region '_cond2.mat'])
            grp_pacx_cond2(:,:, ptnt) = nanmean(pacz,3);
            
            load(['pacz_' phase_region '_' power_region '_cond3.mat'])
            grp_pacx_cond3(:,:, ptnt) = nanmean(pacz,3);
            
            load(['pacz_' phase_region '_' power_region '_cond4.mat'])
            grp_pacx_cond4(:,:, ptnt) = nanmean(pacz,3);
            
            end
            
            
bar_vector_mn1 = [mean(mean(squeeze(nanmean(grp_pacx_cond1(power_freqs<rng1,phas_freqs>rng2,:),1)),1))...
    mean(mean(squeeze(nanmean(grp_pacx_cond2(power_freqs<rng1,phas_freqs>rng2,:),1)),1))...
    mean(mean(squeeze(nanmean(grp_pacx_cond3(power_freqs<rng1,phas_freqs>rng2,:),1)),1))...
    mean(mean(squeeze(nanmean(grp_pacx_cond4(power_freqs<rng1,phas_freqs>rng2,:),1)),1))];


bar_vector_std1 = [nanstd(mean(squeeze(nanmean(grp_pacx_cond1(power_freqs<rng1,phas_freqs>rng2,:),1)),1),0)...
    nanstd(mean(squeeze(nanmean(grp_pacx_cond2(power_freqs<rng1,phas_freqs>rng2,:),1)),1),0)...
    nanstd(mean(squeeze(nanmean(grp_pacx_cond3(power_freqs<rng1,phas_freqs>rng2,:),1)),1),0)...
    nanstd(mean(squeeze(nanmean(grp_pacx_cond4(power_freqs<rng1,phas_freqs>rng2,:),1)),1),0)];


figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std1/sqrt(ptnt), 'rx')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
title([['CFC: phase of ' phase_region '(' num2str(rng2) ':8hz)' ' vs power ' '(42:' num2str(rng1) ')' 'of ' power_region  ]])
ylim([-.4 .57])
ylim([-.3 .4])
ylabel(['mean PACz'])
set(gca, 'FontSize', 16, 'FontWeight', 'bold')           
            
            
   clear bar_vector_mn1  bar_vector_std1        
            
            
        end
        
        
    end
    
    %%
saveas(gcf, 'temp.bmp')
