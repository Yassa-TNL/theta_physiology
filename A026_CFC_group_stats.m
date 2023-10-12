
cond1a      = cell(1,length(group));
cond2a      = cell(1,length(group));
cond3a      = cell(1,length(group));
cond4a      = cell(1,length(group));

cond1a_null = cell(1,length(group));
cond2a_null = cell(1,length(group));
cond3a_null = cell(1,length(group));
cond4a_null = cell(1,length(group));

cond2a_pacz = cell(1,length(group));
cond3a_pacz = cell(1,length(group));
cntr = 0
for reg1 = 1:length(phase_regions_list) %  loop thru pairs of regions
    phase_region = phase_regions_list{reg1};
    
    for reg2 = 1:length(power_regions_list)
        power_region = power_regions_list{reg2};
        
        for ptnt= 1:length(group) % loop thru ptnts
            cntr = cntr+1
            % go to patient dir
            cd(['/mnt/yassamri/iEEG/sandra/subj_' group{ptnt} '/figures/cfc'])
            
            % cond1
            load(['pacz_' phase_region '_' power_region '_cond1.mat'])
            mean_mtx = squeeze(nanmean(Perm_PAC,3));
            std_mtx = squeeze(nanstd(Perm_PAC,0,3));
            raw_mtx = (pacz.*std_mtx)+mean_mtx;
            cond1a{ptnt} = raw_mtx;
            cond1a_null{ptnt} = Perm_PAC;

            % cond2
            load(['pacz_' phase_region '_' power_region '_cond2.mat'])
            cond2a_pacz{ptnt} = pacz;
            mean_mtx = squeeze(nanmean(Perm_PAC,3));
            std_mtx = squeeze(nanstd(Perm_PAC,0,3));
            raw_mtx = (pacz.*std_mtx)+mean_mtx;
            cond2a{ptnt} = raw_mtx;
            cond2a_null{ptnt} = Perm_PAC;

            % cond3
            load(['pacz_' phase_region '_' power_region '_cond3.mat'])
            cond3a_pacz{ptnt} = pacz;
            mean_mtx = squeeze(nanmean(Perm_PAC,3));
            std_mtx = squeeze(nanstd(Perm_PAC,0,3));
            raw_mtx = (pacz.*std_mtx)+mean_mtx;
            cond3a{ptnt} = raw_mtx;
            cond3a_null{ptnt} = Perm_PAC;
            
            % cond4
            load(['pacz_' phase_region '_' power_region '_cond4.mat'])
            mean_mtx = squeeze(nanmean(Perm_PAC,3));
            std_mtx = squeeze(nanstd(Perm_PAC,0,3));
            raw_mtx = (pacz.*std_mtx)+mean_mtx;
            cond4a{ptnt} = raw_mtx;
            cond4a_null{ptnt} = Perm_PAC;
            
        end
    end
end

cond2_pacz= nanmean(cat(3,cond2a_pacz{:}),3);
cond3_pacz= nanmean(cat(3,cond3a_pacz{:}),3);


cond1 = nanmean(cat(3,cond1a{:}),3);
cond2 = nanmean(cat(3,cond2a{:}),3);
cond3 = nanmean(cat(3,cond3a{:}),3);
cond4 = nanmean(cat(3,cond4a{:}),3);

cond1_null_mean = squeeze(nanmean(cat(4,cond1a_null{:}),4));
cond1_null_std  = squeeze(nanstd(cat(4,cond1a_null{:}),0,4));

cond2_null_mean = squeeze(nanmean(cat(4,cond2a_null{:}),4));
cond2_null_std  = squeeze(nanstd(cat(4,cond2a_null{:}),0,4));

cond3_null_mean = squeeze(nanmean(cat(4,cond3a_null{:}),4));
cond3_null_std  = squeeze(nanstd(cat(4,cond3a_null{:}),0,4));

cond4_null_mean = squeeze(nanmean(cat(4,cond4a_null{:}),4));
cond4_null_std  = squeeze(nanstd(cat(4,cond4a_null{:}),0,4));
%
cd('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')

for ph = 1:size(cond1_null_mean,1)
    for po = 1:size(cond1_null_mean,2)
        p_mat_cond1(ph,po) = invprctile([squeeze(cond1_null_mean(ph,po,:))],cond1(ph,po));
        p_mat_cond2(ph,po) = invprctile([squeeze(cond2_null_mean(ph,po,:))],cond2(ph,po));
        p_mat_cond3(ph,po) = invprctile([squeeze(cond3_null_mean(ph,po,:))],cond3(ph,po));
        p_mat_cond4(ph,po) = invprctile([squeeze(cond4_null_mean(ph,po,:))],cond4(ph,po));

    end
end

p_mat_cond1 = (100-p_mat_cond1)/100;
p_mat_cond1(p_mat_cond1>.05)=nan;

p_mat_cond2 = (100-p_mat_cond2)/100;
p_mat_cond2(p_mat_cond2>.05)=nan;

p_mat_cond3 = (100-p_mat_cond3)/100;
p_mat_cond3(p_mat_cond3>.05)=nan;

p_mat_cond4 = (100-p_mat_cond4)/100;
p_mat_cond4(p_mat_cond4>.05)=nan;

%% plot traces as function of low freq phase
figure
plot(phas_freqs,nanmean(cond2_pacz,1), '-k','LineWidth', 3);

hold on
plot(phas_freqs,nanmean(cond3_pacz,1),'-r', 'LineWidth', 3);
xlim([phas_freqs(1) phas_freqs(end)])
legend({'lure-', 'lure+'})
xlabel(['Freq for Phase ']);
ylabel(['PACz value'])
title(['PAC: phase of ' phase_region ' vs. power of ' power_region]) 
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

%% plot p-maps

mx = 0
mn = -.05
figure;
subplot(4,1,1)
imagesc(-1*p_mat_cond1);
set(gca, 'xtick', 1:2:length(phas_freqs));
x=floor(phas_freqs*10)/10;
y=floor(power_freqs*10)/10
set(gca,'xticklabel', x(1:2:length(phas_freqs)),'YDir','normal', 'ytick', 1:2:length(power_freqs),'yticklabel', y(1:2:length(y)),'FontSize', 14, 'FontWeight', 'bold');
xlabel(['Freq for Phase ']);
ylabel(['Freq for Power ' ]);
title('repeat +')
colorbar
caxis([mn mx])
colormap jet

subplot(4,1,2)
imagesc(-1*p_mat_cond2);
set(gca, 'xtick', 1:2:length(phas_freqs));
x=floor(phas_freqs*10)/10;
y=floor(power_freqs*10)/10
set(gca,'xticklabel', x(1:2:length(phas_freqs)),'YDir','normal', 'ytick', 1:2:length(power_freqs),'yticklabel', y(1:2:length(y)),'FontSize', 14, 'FontWeight', 'bold');
xlabel(['Freq for Phase ']);
ylabel(['Freq for Power ' ]);
title('lure -')
colorbar
caxis([mn mx])
colormap jet

subplot(4,1,3)
imagesc(-1*p_mat_cond3);
set(gca, 'xtick', 1:2:length(phas_freqs));
x=floor(phas_freqs*10)/10;
y=floor(power_freqs*10)/10
set(gca,'xticklabel', x(1:2:length(phas_freqs)),'YDir','normal', 'ytick', 1:2:length(power_freqs),'yticklabel', y(1:2:length(y)),'FontSize', 14, 'FontWeight', 'bold');
xlabel(['Freq for Phase ']);
ylabel(['Freq for Power ' ]);
title('lure +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
colorbar
caxis([mn mx])
colormap jet

subplot(4,1,4)
imagesc(-1*p_mat_cond4);
set(gca, 'xtick', 1:2:length(phas_freqs));
x=floor(phas_freqs*10)/10;
y=floor(power_freqs*10)/10
set(gca,'xticklabel', x(1:2:length(phas_freqs)),'YDir','normal', 'ytick', 1:2:length(power_freqs),'yticklabel', y(1:2:length(y)),'FontSize', 14, 'FontWeight', 'bold');
xlabel(['Freq for Phase ']);
ylabel(['Freq for Power ' ]);
title('new +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
colorbar
caxis([mn mx])
colormap jet

