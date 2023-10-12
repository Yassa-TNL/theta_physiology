%clear all

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

% run stats-pool all matrices from all elecs from all patients, AND get raw non z-scored values

cond1a      = cell(1,length(group)*length(phase_regions_list));
cond2a      = cell(1,length(group)*length(phase_regions_list));
cond3a      = cell(1,length(group)*length(phase_regions_list));
cond4a      = cell(1,length(group)*length(phase_regions_list));

cond1a_null = cell(1,length(group)*length(phase_regions_list));
cond2a_null = cell(1,length(group)*length(phase_regions_list));
cond3a_null = cell(1,length(group)*length(phase_regions_list));
cond4a_null = cell(1,length(group)*length(phase_regions_list));

cond1a_pacz = cell(1,length(group)*length(phase_regions_list));
cond2a_pacz = cell(1,length(group)*length(phase_regions_list));
cond3a_pacz = cell(1,length(group)*length(phase_regions_list));
cond4a_pacz = cell(1,length(group)*length(phase_regions_list));

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
            cond1a_pacz{cntr} = pacz;
            mean_mtx = squeeze(nanmean(Perm_PAC,3));
            std_mtx = squeeze(nanstd(Perm_PAC,0,3));
            raw_mtx = (pacz.*std_mtx)+mean_mtx;
            cond1a{cntr} = raw_mtx;
            cond1a_null{cntr} = Perm_PAC;

            % cond2
            load(['pacz_' phase_region '_' power_region '_cond2.mat'])
            cond2a_pacz{cntr} = pacz;
            mean_mtx = squeeze(nanmean(Perm_PAC,3));
            std_mtx = squeeze(nanstd(Perm_PAC,0,3));
            raw_mtx = (pacz.*std_mtx)+mean_mtx;
            cond2a{cntr} = raw_mtx;
            cond2a_null{cntr} = Perm_PAC;

            % cond3
            load(['pacz_' phase_region '_' power_region '_cond3.mat'])
            cond3a_pacz{cntr} = pacz;
            mean_mtx = squeeze(nanmean(Perm_PAC,3));
            std_mtx = squeeze(nanstd(Perm_PAC,0,3));
            raw_mtx = (pacz.*std_mtx)+mean_mtx;
            cond3a{cntr} = raw_mtx;
            cond3a_null{cntr} = Perm_PAC;
            
            % cond4
            load(['pacz_' phase_region '_' power_region '_cond4.mat'])
            cond4a_pacz{cntr} = pacz;
            mean_mtx = squeeze(nanmean(Perm_PAC,3));
            std_mtx = squeeze(nanstd(Perm_PAC,0,3));
            raw_mtx = (pacz.*std_mtx)+mean_mtx;
            cond4a{cntr} = raw_mtx;
            cond4a_null{cntr} = Perm_PAC;
            
        end
    end
end

% z-scored mtx pooled, avg
cond1_pacz= nanmean(cat(3,cond1a_pacz{:}),3);
cond2_pacz= nanmean(cat(3,cond2a_pacz{:}),3);
cond3_pacz= nanmean(cat(3,cond3a_pacz{:}),3);
cond4_pacz= nanmean(cat(3,cond4a_pacz{:}),3);

%% plot traces as function of low freq phase
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
figure
hold on
temp1 = cat(3,cond2a_pacz{:});
temp  = squeeze(nanmean(temp1,1))';
stdshade(temp,.1,'k',phas_freqs, [],[])
h1=plot(phas_freqs,nanmean(nanmean(cond2_pacz,3),1),'-k', 'LineWidth', 3);

temp1 = cat(3,cond3a_pacz{:});
temp  = squeeze(nanmean(temp1,1))';stdshade(temp,.1,'r',phas_freqs, [],[])
h2    = plot(phas_freqs,nanmean(nanmean(cond3_pacz,3),1),'-r', 'LineWidth', 3);


xlim([phas_freqs(1) phas_freqs(end)])
%ylim([0 .55])
xlabel(['Freq for Phase ']);
ylabel(['PACz value'])
title([' Functional Coupling between ' phase_regions_list{:} ' and ' power_regions_list{:} ]) 
set(gca, 'FontSize', 14, 'FontWeight', 'bold')


% stats from cluster
stat_vec = logical(nanmean(zmapthresh,1))
find(stat_vec)

% TEMP HC
% line([phas_freqs(2) phas_freqs(5)],[.6 .6], 'color','b','LineWidth', 3)
% line([phas_freqs(7) phas_freqs(9)],[.6 .6], 'color','b','LineWidth', 3)

% CING HC
% line([phas_freqs(3) phas_freqs(5)],[.6 .6], 'color','b','LineWidth', 3)
% line([phas_freqs(9) phas_freqs(13)],[.6 .6], 'color','b','LineWidth', 3)


% FRO HC
% line([phas_freqs(1) phas_freqs(6)],[.6 .6], 'color','b','LineWidth', 3)
% line([phas_freqs(8) phas_freqs(11)],[.6 .6], 'color','b','LineWidth', 3)
% line([phas_freqs(14) phas_freqs(18)],[.6 .6], 'color','b','LineWidth', 3)


% OFC HC
%  line([phas_freqs(1) phas_freqs(18)],[1 1], 'color','b','LineWidth', 3)

% fro ca3
% line([phas_freqs(1) phas_freqs(6)],[.5 .5], 'color','b','LineWidth', 3)
% line([phas_freqs(8) phas_freqs(11)],[.5 .5], 'color','b','LineWidth', 3)
% line([phas_freqs(13) phas_freqs(18)],[.5 .5], 'color','b','LineWidth', 3)

% temp ca3
% line([phas_freqs(2) phas_freqs(5)],[.5 .5], 'color','b','LineWidth', 3)
% line([phas_freqs(16) phas_freqs(16)],[.5 .5], 'color','b','LineWidth', 3)

% cing ca3
% line([phas_freqs(1) phas_freqs(2)],[.9 .9], 'color','b','LineWidth', 3)
% line([phas_freqs(9) phas_freqs(11)],[.9 .9], 'color','b','LineWidth', 3)
% line([phas_freqs(14) phas_freqs(18)],[.9 .9], 'color','b','LineWidth', 3)


cd('/mnt/yassamri/iEEG/sandra/neuroblitz_figures')
saveas(gcf, ['CFC_trace_' phase_regions_list{:} '_' power_regions_list{:}  '.bmp'])

%% plot average with all chans pooled together

mx = 1.2
mn = -.8

figure;
subplot(4,1,1)
temp1 = cat(3,cond1a_pacz{:});
imagesc(nanmean(temp1,3));
phase_tickmarks = 1:2:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:2:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'YDir', 'normal')
ylabel(['Freq for Power in ' power_regions_list{:}]);
title('repeat +')
set(gca, 'FontSize', 13, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar

subplot(4,1,2)
temp2 = cat(3,cond2a_pacz{:});
imagesc(nanmean(temp2,3));
phase_tickmarks = 1:2:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:2:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'YDir', 'normal')
ylabel(['Freq for Power in ' power_regions_list{:}]);
title('lure -')
set(gca, 'FontSize', 13, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar

subplot(4,1,3)
temp3 = cat(3,cond3a_pacz{:});
imagesc(nanmean(temp3,3));
phase_tickmarks = 1:2:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:2:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'YDir', 'normal')
ylabel(['Freq for Power in ' power_regions_list{:}]);
title('lure +')
set(gca, 'FontSize', 13, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar

subplot(4,1,4)
temp4 = cat(3,cond4a_pacz{:});
imagesc(nanmean(temp4,3));
phase_tickmarks = 1:2:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:2:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'YDir', 'normal')
xlabel(['Freq for Phase in ' phase_regions_list{:}]);
ylabel(['Freq for Power in ' power_regions_list{:}]);
title('new +')
caxis([mn mx])
set(gca, 'FontSize', 13, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar

suptitle([' Functional Coupling between ' phase_regions_list{:} ' and ' power_regions_list{:} ])
cd('/mnt/yassamri/iEEG/sandra/neuroblitz_figures')
saveas(gcf, ['CFC_heatmap_' phase_regions_list{:} '_' power_regions_list{:}  '.bmp'])




%% % get p-mtx

% raw mtx pooled, avg
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
p_mat_cond1(p_mat_cond1>(.05/length(phas_freqs)))=nan;

p_mat_cond2 = (100-p_mat_cond2)/100;
p_mat_cond2(p_mat_cond2>(.05/length(phas_freqs)))=nan;

p_mat_cond3 = (100-p_mat_cond3)/100;
p_mat_cond3(p_mat_cond3>(.05/length(phas_freqs)))=nan;

p_mat_cond4 = (100-p_mat_cond4)/100;
p_mat_cond4(p_mat_cond4>(.05/length(phas_freqs)))=nan;



%% plot p-maps

mx = 0
mn = -.02
figure;
subplot(4,1,1)
imagesc(-1*p_mat_cond1);
phase_tickmarks = 1:2:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:2:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'YDir', 'normal')
xlabel(['Freq for Phase in NC']);
ylabel(['Freq for Power in HC']);
title('repeat +')
set(gca, 'FontSize', 13, 'FontWeight', 'bold')
colorbar
%colormap jet
caxis([mn mx])
colorbar




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
%colormap jet

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
%colormap jet

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
%colormap jet




         %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% one val from each subj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get one avg mtx frm each subj
cond1a      = cell(1,length(phase_regions_list));
cond2a      = cell(1,length(phase_regions_list));
cond3a      = cell(1,length(phase_regions_list));
cond4a      = cell(1,length(phase_regions_list));

for reg1 = 1:length(phase_regions_list) %  loop thru pairs of regions
    phase_region = phase_regions_list{reg1};
    
    for reg2 = 1:length(power_regions_list)
        power_region = power_regions_list{reg2};
        
        cd('/mnt/yassamri/iEEG/sandra/groupdata_cfc')

        
        load(['group_pacz_' phase_region '_HC_cond1.mat'])
        cond1a{reg1} = grp_pacx_cond1;
        
        load(['group_pacz_' phase_region '_HC_cond2.mat'])
        cond2a{reg1} = grp_pacx_cond2;
        
        load(['group_pacz_' phase_region '_HC_cond3.mat'])
        cond3a{reg1} = grp_pacx_cond3;
        
        load(['group_pacz_' phase_region '_HC_cond4.mat'])
        cond4a{reg1} = grp_pacx_cond4;
    end
end

cond1 = cat(3,cond1a{:});
cond2 = cat(3,cond2a{:});
cond3 = cat(3,cond3a{:});
cond4 = cat(3,cond4a{:});


%% plot average, with one matrix per subj per region

mx = 1
mn = -.3
figure;

subplot(4,1,1)
imagesc(nanmean(cond1,3));
phase_tickmarks = 1:2:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:2:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'YDir', 'normal')
xlabel(['Freq for Phase in NC']);
ylabel(['Freq for Power in HC']);
title('repeat +')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar

subplot(4,1,2)
imagesc(nanmean(cond2,3));
phase_tickmarks = 1:2:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:2:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'YDir', 'normal')
xlabel(['Freq for Phase in NC']);
ylabel(['Freq for Power in HC']);
title('lure -')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar


subplot(4,1,3)
imagesc(nanmean(cond3,3));
phase_tickmarks = 1:2:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:2:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'YDir', 'normal')
xlabel(['Freq for Phase in NC']);
ylabel(['Freq for Power in HC']);
title('lure +')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar



subplot(4,1,4)
imagesc(nanmean(cond4,3));
phase_tickmarks = 1:2:length(phas_freqs);
set(gca, 'XTick',phase_tickmarks,'XTickLabel',round(phas_freqs(phase_tickmarks)))
power_tickmarks = 1:2:length(power_freqs);
set(gca, 'YTick',power_tickmarks,'YTickLabel',round(power_freqs(power_tickmarks)))
set(gca,'YDir', 'normal')
xlabel(['Freq for Phase in NC']);
ylabel(['Freq for Power in HC']);
title('New +')
caxis([mn mx])
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
colorbar
colormap jet
caxis([mn mx])
colorbar

suptitle([' '])


%% plot traces as function of low freq phase - 1 datapoint per subj

figure
hold on

temp = squeeze(nanmean(cond2,1))';
stdshade(temp,.1,'k',phas_freqs, [],[])
h1=plot(phas_freqs,nanmean(nanmean(cond2,3),1),'-k', 'LineWidth', 3);

temp = squeeze(nanmean(cond3,1))';
stdshade(temp,.1,'r',phas_freqs, [],[])
h2=plot(phas_freqs,nanmean(nanmean(cond3,3),1),'-r', 'LineWidth', 3);


xlim([phas_freqs(1) phas_freqs(end)])
xlabel(['Freq for Phase ']);
ylabel(['PACz value'])
title(['PAC: phase of neocotex vs. power of hippocampus']) 
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
saveas(gcf, 'temp.jpg')
