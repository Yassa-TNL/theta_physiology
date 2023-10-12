% draw time courses for onset timing for OFC, FRO, TEMP, CA3

addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
cd('/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_onset')

grp = {'39' '57' '44' '63' '66'}
%theta = [2.9 5.5; 2.9 5; 3 4;3 5 ; 3.5 4.9] %OFC, FRO, temp, MTL, CA3
theta = [3 5; 3 5;3 5;3 5 ;3 5]

freq_name = []
win = []
stim_onset=301
% OFC
traces_OFC_subjs = nan(5, data_length-stim_onset+1);
for ptnt = 1:length(grp)
    load(['subj' grp{ptnt} 'spectrograms.mat'])
    temp = nanmean(OFC_cond3,3);
    traces_OFC_subjs(ptnt,:) = nanmean(temp(freq>theta(1,1) & freq<theta(1,2),stim_onset:end),1);
    clear temp
    clear OFC_cond3
end

% FRO
traces_fro_subjs =nan(5, data_length-stim_onset+1);
for ptnt = 1:length(grp)
    load(['subj' grp{ptnt} 'spectrograms.mat'])
    temp = nanmean(fro_cond3,3);
    traces_fro_subjs(ptnt,:) = nanmean(temp(freq>theta(2,1) & freq<theta(2,2),stim_onset:end),1);
    clear temp
    clear fro_cond3
end

% temp
traces_temp_subjs = nan(5, data_length-stim_onset+1);
for ptnt = 1:length(grp)
    load(['subj' grp{ptnt} 'spectrograms.mat'])
    temp = nanmean(temp_cond3,3);
    traces_temp_subjs(ptnt,:) = nanmean(temp(freq>theta(3,1) & freq<theta(3,2),stim_onset:end),1);
    clear temp
    clear temp_cond3
end

% MTL
traces_MTL_subjs = nan(5, data_length-stim_onset+1);
for ptnt = 1:length(grp)
    load(['subj' grp{ptnt} 'spectrograms.mat'])
    temp = nanmean(MTL_cond3,3);
    traces_MTL_subjs(ptnt,:) = nanmean(temp(freq>theta(4,1) & freq<theta(4,2),stim_onset:end),1);
    clear temp
    clear MTL_cond3
end

% CA3
traces_CA3_subjs =nan(5, data_length-stim_onset+1);
for ptnt = 1:length(grp)
    load(['subj' grp{ptnt} 'spectrograms.mat'])
    temp = nanmean(CA3_cond3,3);
    traces_CA3_subjs(ptnt,:) = nanmean(temp(freq>theta(5,1) & freq<theta(5,2),stim_onset:end),1);
    clear temp
    clear CA3_cond3
end

%
figure
hold on
stdshade(traces_OFC_subjs,.1,'g',[],[] ,freq_name,win)
h1 = plot(nanmean(traces_OFC_subjs,1), 'g', 'LineWidth', 3)

stdshade(traces_fro_subjs,.1,'b',[],[] ,freq_name,win)
h2 = plot(nanmean(traces_fro_subjs,1), 'b', 'LineWidth', 3)

stdshade(traces_temp_subjs,.1,'m',[],[] ,freq_name,win)
h3 = plot(nanmean(traces_temp_subjs,1), 'm', 'LineWidth', 3)

stdshade(traces_MTL_subjs,.1,'r',[],[] ,freq_name,win)
h4 = plot(nanmean(traces_MTL_subjs,1), 'r', 'LineWidth', 2)

stdshade(traces_CA3_subjs,.1,'k',[],[] ,freq_name,win)
h5 = plot(nanmean(traces_CA3_subjs,1), 'k', 'LineWidth', 2)

%legend([h1 h2 h3 h4 h5],'OFC', 'FRO', 'TEMP', 'MTL', 'CA3')
legend([h1 h2 h3],'OFC', 'Frontal Cortex', 'Temporal Cortex')

set(gca, 'FontSize', 13, 'FontWeight', 'bold', 'FontName', 'Arial')
title(['         Power in 3-5 Hz       ' ; 'Correct Identification of Lures'])
ylabel('z-score power')
xlabel('time')
xlim([0 1500])


saveas(gcf, 'group_time_traces_onset.bmp') 


