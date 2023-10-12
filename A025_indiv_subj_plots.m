%close all
subj = '44'
reg1 = 'OFC'
reg2 = 'HC'
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/figures/cfc'])
phase_region = reg1
power_region = reg2


mn = -.8
mx =.85


h=figure; 
suptitle(['subj:' subj ' ' reg1 ' vs ' reg2])
subplot(4,1,1)
load(['pacz_' reg1 '_' reg2 '_cond1.mat'])
imagesc(squeeze(nanmean(pacz,3))); 
phase_tickmarks = 1:4:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:4:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'XDir', 'normal')
xlabel(['Freq for Phase in ' reg1]);
ylabel(['Freq for Power in ' reg2]);
title('repeat +')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar

subplot(4,1,2)
load(['pacz_' reg1 '_' reg2 '_cond2.mat'])
imagesc(squeeze(nanmean(pacz,3)));
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'XDir', 'normal')
xlabel(['Freq for Phase in ' phase_region]);
ylabel(['Freq for Power in ' power_region]);
title('lure -')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar

subplot(4,1,3)
load(['pacz_' reg1 '_' reg2 '_cond3.mat'])
imagesc(squeeze(nanmean(pacz,3)));
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'XDir', 'normal')
xlabel(['Freq for Phase in ' phase_region]);
ylabel(['Freq for Power in ' power_region]);
title('lure +')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar

subplot(4,1,4)
load(['pacz_' reg1 '_' reg2 '_cond4.mat'])
imagesc(squeeze(nanmean(pacz,3)));
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'XDir', 'normal')
xlabel(['Freq for Phase in ' phase_region]);
ylabel(['Freq for Power in ' power_region]);title('new')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar
%%
saveas(gcf, 'temp.jpg')