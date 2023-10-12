
x_reg = 'INS'
y_reg = 'CA3'

% get regional elecs
[fro_chan_idx,MTL_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,OFC_chan_idx,CA3_chan_idx,...
    CA1_chan_idx,HC_chan_idx, ACC_chan_idx] = get_elecs(subj);
if strcmp('OFC',x_reg)
    x_elecs = OFC_chan_idx;
elseif strcmp('FRO',x_reg)
    x_elecs = fro_chan_idx;
elseif strcmp('TEMP',x_reg)
    x_elecs = temp_chan_idx;
elseif strcmp('CING',x_reg)
    x_elecs = cingulate_chan_idx;
elseif strcmp('HC',x_reg)
    x_elecs = HC_chan_idx;
elseif strcmp('INS',x_reg)
    x_elecs = insula_chan_idx;
end

if strcmp('MTL',y_reg)
    y_elecs = MTL_chan_idx;
elseif  strcmp('HC',y_reg)
    y_elecs = HC_chan_idx;
elseif strcmp('CA3',y_reg)
    y_elecs = CA3_chan_idx;
end

% data
if strcmp('tuning_correct', exp_type)
    cond_2_x = norm_freq_acrs_chan_cond_2(:,501:1501,:,x_elecs);
    cond_3_x = norm_freq_acrs_chan_cond_3(:,501:1501,:,x_elecs);
    
    cond_2_y = norm_freq_acrs_chan_cond_2(:,501:1501,:,y_elecs);
    cond_3_y = norm_freq_acrs_chan_cond_3(:,501:1501,:,y_elecs);
end

% get y_axis cluster
cd(['/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_' lock '_' exp_type '/cluster_matrices'])
load('CA3_onset_group_cluster_delta.mat')
y_axis_clust = zmapthresh;
load('OFC_onset_group_cluster_deltatheta.mat')
%%
x_axis_clust = zmapthresh;

% 1 value oper trials

cond_2_x = nanmean(cond_2_x,4);
cond_2_y = nanmean(cond_2_y,4);
cond_3_x = nanmean(cond_3_x,4);
cond_3_y = nanmean(cond_3_y,4);

% cond2
trls_2 = nan(size(cond_2_x,3),2);
for trl = 1:size(cond_2_x,3)
    tempx = squeeze(cond_2_x(:,:, trl));
    tempy = squeeze(cond_2_y(:,:, trl));
    trls_2(trl,:) = [nanmean(nanmean(tempx(logical(x_axis_clust)),2)) nanmean(nanmean(tempy(logical(y_axis_clust)),2))];
end


% cond3
trls_3 = nan(size(cond_3_x,3),2);
for trl = 1:size(cond_3_x,3)
    tempx = squeeze(cond_3_x(:,:, trl));
    tempy = squeeze(cond_3_y(:,:, trl));
    trls_3(trl,:) = [nanmean(nanmean(tempx(logical(x_axis_clust)),2)) nanmean(nanmean(tempy(logical(y_axis_clust)),2))];
end


all_trial = [trls_2 ; trls_3]
cntr = 0
for a = 1:size(all_trial,1)
    if ~isnan(all_trial(a,:))
        cntr=cntr+1;
        temp(cntr,:) = all_trial(a,:);
    end
end
[R, P ]= corrcoef(temp, 'Alpha', .05);

figure
subplot(1,2,1)
plot(trls_2(:,1), trls_2(:,2), 'b*')
hold on
plot(trls_3(:,1), trls_3(:,2), 'r*')
xlabel([x_reg ' cluster'])
ylabel([y_reg ' cluster'])
title(['R = ' num2str(R(2)) ', P = ' num2str(P(2)) ])
legend({'lure -' 'lure +'})
set(gca, 'FontSize',15, 'FontWeight', 'bold')
saveas(gcf,'temp.jpg')

subplot(1,2,2)
imagesc(x_axis_clust)
hold on
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
title('Cluster-corrected Z map')
xlabel('Time (s)'), ylabel('Frequency (Hz)')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
%%
saveas(gcf,'temp.jpg')

%% keep indiv trials
% cond 2 
trls_2 = nan(size(cond_2_x,4)*size(cond_2_y,4)*size(cond_2_y,3),2);
cntr=0
for ex = 1:size(cond_2_x,4)
    for ey = 1:size(cond_2_y,4)
        for t = 1:size(cond_2_x,3)
            cntr = cntr+1;
            tempx = squeeze(cond_2_x(:,:, t, ex));
            tempy = squeeze(cond_2_y(:,:, t, ey));
            trls_2(cntr,:) =[nanmean(nanmean(tempx(logical(x_axis_clust)),2)) nanmean(nanmean(tempy(logical(y_axis_clust)),2))];
            clear tempx tempy
        end
    end
end
% cond 3
trls_3 = nan(size(cond_3_x,4)*size(cond_3_y,4)*size(cond_3_x,3),2);
cntr=0
for ex = 1:size(cond_3_x,4)
    for ey = 1:size(cond_3_y,4)
        for t = 1:size(cond_3_x,3)
            cntr = cntr+1
            tempx = squeeze(cond_3_x(:,:, t, ex));
            tempy = squeeze(cond_3_y(:,:, t, ey));
            trls_3(cntr,:) =[nanmean(nanmean(tempx(logical(x_axis_clust)),2)) nanmean(nanmean(tempy(logical(y_axis_clust)),2))];
            clear tempx tempy
        end
    end
end

all_trial = [trls_2 ; trls_3]
cntr = 0
for a = 1:size(all_trial,1)
    if ~isnan(all_trial(a,:))
        cntr=cntr+1;
        temp(cntr,:) = all_trial(a,:);
    end
end
[R, P ]= corrcoef(temp, 'Alpha', .05);

figure
plot(trls_2(:,1), trls_2(:,2), 'b*')
hold on
plot(trls_3(:,1), trls_3(:,2), 'r*')
xlabel([x_reg ' cluster'])
ylabel('CA3 cluster')
title(['R = ' num2str(R(2)) ', P = ' num2str(P(2)) ])
legend({'lure -' 'lure +'})
set(gca, 'FontSize',15, 'FontWeight', 'bold')
saveas(gcf,'temp.jpg')