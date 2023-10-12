close all
mn  = -.4
mx  = .99
reg = 'CA1' % HC CA3 CA1 OFC FRO TEMP CING INS EC
clear chan_idx
if strcmp('OFC',reg)
    chan_idx = OFC_chan_idx;
elseif strcmp('FRO',reg)
    chan_idx = fro_chan_idx;
elseif strcmp('CING',reg)
    chan_idx = cingulate_chan_idx;
elseif strcmp('TEMP',reg)
    chan_idx = temp_chan_idx;
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
%%
cond1 = nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_1(:, edge_points+1:end-edge_points,:, chan_idx),3)),3);
cond2 = nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_2(:, edge_points+1:end-edge_points,:, chan_idx),3)),3);
cond3 = nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_3(:, edge_points+1:end-edge_points,:, chan_idx),3)),3);
cond4 = nanmean(squeeze(nanmean(norm_freq_acrs_chan_cond_4(:, edge_points+1:end-edge_points,:, chan_idx),3)),3);



figure
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), cond1)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), cond2)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), cond3)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure +')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), cond4)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')

suptitle(['subj:' num2str(subj) ' - region:' reg])
