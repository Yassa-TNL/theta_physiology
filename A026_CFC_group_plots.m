
phase_regions_list = {'TEMP'} 
power_regions_list = {'CA3'}  

%% initializa group matrices for each pair of regions and each condition

    % cortex w/ MTL
if strcmp('OFC', phase_regions_list) | strcmp('FRO', phase_regions_list) | strcmp('TEMP', phase_regions_list) && strcmp('MTL', power_regions_list)
    group = {'39' '57'  '66'  '44' '63' '84'}

elseif strcmp('CING', phase_regions_list) && strcmp('MTL', power_regions_list)
    group = {'39' '57'  '66'  '44' '63'}
    
    % cortex w/ HC and cortex w/ CA3
elseif strcmp('OFC', phase_regions_list) | strcmp('FRO', phase_regions_list) | strcmp('TEMP', phase_regions_list) | strcmp('CING', phase_regions_list) && strcmp('HC', power_regions_list) | strcmp('CA3', power_regions_list)
    group = {'39' '57'  '66'  '44' }
end


cd(['/mnt/yassamri/iEEG/sandra/subj_39/figures/cfc'])
load('pacz_FRO_MTL_cond1.mat')
grp_pacx_cond1 = zeros(size(pacz,1), size(pacz,2), length(group));
grp_pacx_cond2 = zeros(size(pacz,1), size(pacz,2), length(group));
grp_pacx_cond3 = zeros(size(pacz,1), size(pacz,2), length(group));
grp_pacx_cond4 = zeros(size(pacz,1), size(pacz,2), length(group));
%close all
mn = -.3
mx = 2

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
        
        % save grp data
        cd('/mnt/yassamri/iEEG/sandra/groupdata_cfc')
        save(['group_pacz_' phase_region '_' power_region '_cond1'], 'grp_pacx_cond1');
        save(['group_pacz_' phase_region '_' power_region '_cond2'], 'grp_pacx_cond2');
        save(['group_pacz_' phase_region '_' power_region '_cond3'], 'grp_pacx_cond3');
        save(['group_pacz_' phase_region '_' power_region '_cond4'], 'grp_pacx_cond4');
        
        %% plot
        figure;
        suptitle(['CFC: phase of ' phase_region ' vs power of ' power_region ])

        subplot(4,1,1)
        imagesc(nanmean(grp_pacx_cond1,3));
        % imagesc(nanmean(grp_pacx_cond1,3)./nanstd(grp_pacx_cond1,0,3));
        phase_tickmarks = 1:2:length(phas_freqs);
        set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
        power_tickmarks = 1:2:length(power_freqs);
        set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
        set(gca,'YDir', 'normal')
        xlabel(['Freq for Phase in ' phase_region]);
        ylabel(['Freq for Power in ' power_region]);
        title('repeat +')
        set(gca, 'FontSize', 12, 'FontWeight', 'bold')
        colorbar
        colormap jet
        caxis([mn mx])
        colorbar
        
        subplot(4,1,2)
        imagesc(nanmean(grp_pacx_cond2,3));
        %  imagesc(nanmean(grp_pacx_cond2,3)./nanstd(grp_pacx_cond2,0,3));
        set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
        set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
        set(gca,'YDir', 'normal')
        xlabel(['Freq for Phase in ' phase_region]);
        ylabel(['Freq for Power in ' power_region]);
        title('lure -')
        set(gca, 'FontSize', 12, 'FontWeight', 'bold')
        colorbar
        colormap jet
        caxis([mn mx])
        colorbar
        
        subplot(4,1,3)
        imagesc(nanmean(grp_pacx_cond3,3));
        % imagesc(nanmean(grp_pacx_cond3,3)./nanstd(grp_pacx_cond3,0,3));
        set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
        set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
        set(gca,'YDir', 'normal')
        xlabel(['Freq for Phase in ' phase_region]);
        ylabel(['Freq for Power in ' power_region]);
        title('lure +')
        set(gca, 'FontSize', 12, 'FontWeight', 'bold')
        colorbar
        colormap jet
        caxis([mn mx])
        colorbar
        
        subplot(4,1,4)
        imagesc(nanmean(grp_pacx_cond4,3));
        % imagesc(nanmean(grp_pacx_cond4,3)./nanstd(grp_pacx_cond4,0,3));
        set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
        set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
        set(gca,'YDir', 'normal')
        xlabel(['Freq for Phase in ' phase_region]);
        ylabel(['Freq for Power in ' power_region]);title('new')
        set(gca, 'FontSize', 12, 'FontWeight', 'bold')
        colorbar
        colormap jet
        caxis([mn mx])
        colorbar
        % suptitle(['CFC: phase of ' phase_region ' vs power of ' power_region ; '           N = ' num2str(ptnt) '                 '])
        %  saveas(gcf,[ 'group cfc' phase_region '_' power_region '.bmp'])
        
            
        end
           
    end

%%
cd('/mnt/yassamri/iEEG/sandra/neuroblitz_figures')


saveas(gcf, ['CFC_plots_' phase_regions_list{:} '_' power_regions_list{:} '.bmp'])




 