CA3_chan_idx
elec = 1
mn   = -.3
mx   = .6
titles = {'repeat'  'lure -' 'lure +' 'new'}

figure
for cond_cnt = 1:cond_num
    if cond_cnt ==1
        mtx = norm_freq_acrs_chan_cond_1(:, edge_points+1:end-edge_points, :, CA3_chan_idx(elec));
    elseif cond_cnt==2
        mtx = norm_freq_acrs_chan_cond_2(:, edge_points+1:end-edge_points, :, CA3_chan_idx(elec));
    elseif cond_cnt==3
        mtx = norm_freq_acrs_chan_cond_3(:, edge_points+1:end-edge_points, :, CA3_chan_idx(elec));
    elseif cond_cnt ==4
        mtx = norm_freq_acrs_chan_cond_4(:, edge_points+1:end-edge_points, :, CA3_chan_idx(elec));
    end
    subplot (cond_num,1,cond_cnt)
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq),nanmean(mtx,3))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    colormap jet
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title(titles{cond_cnt})
end

tickmarks = 1:2:length(freq);

figure
cond_cnt = 1
mtx = norm_freq_acrs_chan_cond_1(:, edge_points+1:end-edge_points, :, CA3_chan_idx(elec));
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq),nanmean(mtx,3))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(titles{cond_cnt})

figure
cond_cnt = 2
mtx = norm_freq_acrs_chan_cond_2(:, edge_points+1:end-edge_points, :, CA3_chan_idx(elec));
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq),nanmean(mtx,3))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(titles{cond_cnt})

figure
cond_cnt = 3
mtx = norm_freq_acrs_chan_cond_3(:, edge_points+1:end-edge_points, :, CA3_chan_idx(elec));
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq),nanmean(mtx,3))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(titles{cond_cnt})

figure
cond_cnt = 4
mtx = norm_freq_acrs_chan_cond_4(:, edge_points+1:end-edge_points, :, CA3_chan_idx(elec));
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq),nanmean(mtx,3))
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title(titles{cond_cnt})
%%
figure
hold on
for trl = 1:5
plot(linspace(-pre_stim, post_stim, size(cond2,2)), cond3(1,:,CA3_chan_idx(elec)), 'm')
end

plot(-.5:1/fs:2, nanmean(cond1(:,:,chan),1),'b', 'LineWidth', 2)
plot(-.5:1/fs:2, nanmean(cond2(:,:,chan),1),'g', 'LineWidth', 2)
plot(-.5:1/fs:2, nanmean(cond3(:,:,chan),1),'r', 'LineWidth', 2)
plot(-.5:1/fs:2, nanmean(cond4(:,:,chan),1),'m', 'LineWidth', 2)
legend({'repeat', 'incorr lure', 'corr lure', 'new'})