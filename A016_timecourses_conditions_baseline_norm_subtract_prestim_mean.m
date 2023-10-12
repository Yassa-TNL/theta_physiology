%freq_ranges   = [2.9 6]  % low theta
%freq_name =[];
freq_ranges   = [40 200];freq_name = 'g'
strt_time = 500
end_time  = 800
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')

% elecs
if strcmp('57', subj)
    fro_chan_idx       = [1 4:8 13 16:17 25:27 35 37 44 46:47 70:83 89 91 98 99 101 102 108 109 110 112];
    MTL_chan_idx       = [48:50 58:61 67:68 113:115 122:124 132:134];
    temp_chan_idx      = [55 56 65 66 72 73 118:121 128:131 136:140];
    insula_chan_idx    = [29:32 38:40 93:96 103 104]
    cingulate_chan_idx = [2 9:11 74:77 ]
    OFC_chan_idx       = [18:19 84:86]
elseif strcmp('39', subj)
    fro_chan_idx       = [8 16:19 25:26 28 63:64 66 72:74];
    MTL_chan_idx       = [29:32 39:42 51:52 77:79 85:88 94:96];
    temp_chan_idx      = [34:38 45:48 54:57 83:84 91:93 100:103 ];
    insula_chan_idx    = [20 21]
    cingulate_chan_idx = [11 12]
    OFC_chan_idx       = [1 2 58];
    
elseif strcmp('66', subj)
    fro_chan_idx       = [37:40 49:50 95:97 101 108 109 110];
    MTL_chan_idx       = [1:5 11:15 21:23 61:65 71:75 81:83 ];
    temp_chan_idx      = [7:9 16:20 24 25 29 30 48 68:70 77 78 84:90 ];
    insula_chan_idx    = [];
    cingulate_chan_idx = [31:33 35 51:54 91:93 111:114];
    OFC_chan_idx       = [41 42 102];
    
elseif strcmp('63', subj)
    fro_chan_idx       = [45 52:56 116 ];
    MTL_chan_idx       = [5:8 17:19 25:28 67:70 77:80 91:94];
    temp_chan_idx      = [14 22 23 24 31 34 73:76 84:86];
    insula_chan_idx    = [];
    cingulate_chan_idx = [47 49 97:101];
    OFC_chan_idx       = [37:39];
    
elseif strcmp('44', subj)
    fro_chan_idx       = [36:38 44:47 84 :86];
    MTL_chan_idx       = [1:4 11:14 21 22 23 49 58:62 68 70:72];
    temp_chan_idx      = [7:10 16:19 25:28 54:57 64:67 74:76];
    insula_chan_idx    = [];
    cingulate_chan_idx = [29:32];
    OFC_chan_idx       = [39 40 77:80];
    
elseif strcmp('83', subj)
    fro_chan_idx       = [7:10 12 15:20 25:28 32:35 42 44 50 52 29 ];
    MTL_chan_idx       = [];
    temp_chan_idx      = [5 6 21 47];
    insula_chan_idx    = [];
    cingulate_chan_idx = [];
    OFC_chan_idx       = [];
    
elseif strcmp('84', subj)
    fro_chan_idx       = [ 27 29 31 32  28 30 40];
    MTL_chan_idx       = [17:20 49:52];
    temp_chan_idx      = [11:13 15:16 41:44 47 48 55 56];
    insula_chan_idx    = [];
    cingulate_chan_idx = [];
    OFC_chan_idx       = [1 2 3 33 34 36];
end

power_series_mn_cond1 = [];
power_series_mn_cond2 = [];
power_series_mn_cond3 = [];
power_series_mn_cond4 = [];


stim_onset  = 301;
for chan = 1:chan_counter % loop thru chan
    for cntr = 1:size(freq_ranges,1) % go thru 1 freq
        counter_freq_range             = freq>freq_ranges(cntr,1) & freq<freq_ranges(cntr,2);
        
        power_series_mn_cond1 (:,chan) = nanmean(mn_acrs_trials_cond1(counter_freq_range,stim_onset:end,chan),1)-...
                                           nanmean(nanmean(mn_acrs_trials_cond1(counter_freq_range,1:stim_onset,chan),1)); 
        power_series_std_cond1(:,chan) = nanstd(mn_acrs_trials_cond1 (counter_freq_range,stim_onset:end,chan),0,1);
                                          
        
        power_series_mn_cond2 (:,chan) = nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset:end,chan),1)-...
                                    nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,1:stim_onset,chan),1));
        power_series_std_cond2(:,chan) = nanstd(mn_acrs_trials_cond2 (counter_freq_range,stim_onset:end,chan),0,1);

        
        power_series_mn_cond3 (:,chan) = nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset:end,chan),1)-...
                                         nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,1:stim_onset,chan),1));
        power_series_std_cond3(:,chan) = nanstd(mn_acrs_trials_cond3 (counter_freq_range,stim_onset:end,chan),0,1);

        
        power_series_mn_cond4 (:,chan) = nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset:end,chan),1)-...
                                        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,1:stim_onset,chan),1));
        power_series_std_cond4(:,chan) = nanstd(mn_acrs_trials_cond4 (counter_freq_range,stim_onset:end,chan),0,1);
        
    end
end



% find selective theta chans: average over time range after subtracting
% prestim mean


MTL_chan_idx_new=[];
for chan =MTL_chan_idx
    temp = [nanmean(power_series_mn_cond1(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond2(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond3(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond4(strt_time:end_time,chan))];

    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        MTL_chan_idx_new = [MTL_chan_idx_new chan]
    end
end


fro_chan_idx_new=[];
for chan =fro_chan_idx
    temp = [nanmean(power_series_mn_cond1(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond2(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond3(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond4(strt_time:end_time,chan))];

    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        fro_chan_idx_new = [fro_chan_idx_new chan]
    end
end


temp_chan_idx_new=[];
for chan =temp_chan_idx
    temp = [nanmean(power_series_mn_cond1(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond2(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond3(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond4(strt_time:end_time,chan))];

    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        temp_chan_idx_new = [temp_chan_idx_new chan]
    end
end


insula_chan_idx_new=[];
for chan =insula_chan_idx
   temp = [nanmean(power_series_mn_cond1(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond2(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond3(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond4(strt_time:end_time,chan))];
    
    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        insula_chan_idx_new = [insula_chan_idx_new chan]
    end
    clear temp
end

cingulate_chan_idx_new=[];
for chan =cingulate_chan_idx
   temp = [nanmean(power_series_mn_cond1(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond2(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond3(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond4(strt_time:end_time,chan))];
    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
       cingulate_chan_idx_new = [cingulate_chan_idx_new chan]
    end
end

OFC_chan_idx_new=[];
for chan =OFC_chan_idx
   temp = [nanmean(power_series_mn_cond1(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond2(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond3(strt_time:end_time,chan))...
        nanmean(power_series_mn_cond4(strt_time:end_time,chan))];
    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        OFC_chan_idx_new = [OFC_chan_idx_new chan]
    end
end

clc
MTL_chan_idx_new
temp_chan_idx_new
fro_chan_idx_new

% delta/theta
%frontal cortex
figure
suptitle([num2str(freq_ranges(1)) '-' num2str(freq_ranges(2)) ' hz power'])
a = .45
b = -.48

if strcmp('g',freq_name)
    a = .4
    b = -.12
end

% mtl
subplot(3,2,1)
hold on
chan_idx = MTL_chan_idx_new;


temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,[], [])
h1 = plot(nanmean(temp,1), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,[], [])
h2 = plot(nanmean(temp,1), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,[], [])
h3 = plot(nanmean(temp,1), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,[], [])
h4 = plot(nanmean(temp,1), 'm', 'LineWidth', 2)
clear temp

xlim([0 1700])
ylim([b a])
title ('mtl')
xlabel('time from stim onset')
ylabel('z-score power')
%legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')

% temporal lobe
subplot(3,2,2)
hold on

chan_idx = temp_chan_idx_new;


temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,[], [])
h1 = plot(nanmean(temp,1), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,[], [])
h2 = plot(nanmean(temp,1), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,[], [])
h3 = plot(nanmean(temp,1), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,[], [])
h4 = plot(nanmean(temp,1), 'm', 'LineWidth', 2)
clear temp


xlim([0 1700])
ylim([b a])
title ('temp lobe')
xlabel('time from stim onset')
ylabel('z-score power')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')

% fro
subplot(3,2,3)
hold on
chan_idx = fro_chan_idx_new;


temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,[], [])
h1 = plot(nanmean(temp,1), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,[], [])
h2 = plot(nanmean(temp,1), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,[], [])
h3 = plot(nanmean(temp,1), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,[], [])
h4 = plot(nanmean(temp,1), 'm', 'LineWidth', 2)
clear temp

xlim([0 1700])
ylim([b a])
title ('fro')
xlabel('time from stim onset')
ylabel('z-score power')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')


% OFC_chan_idx_new
subplot(3,2,4)
hold on
chan_idx = OFC_chan_idx_new;


temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,[], [])
h1 = plot(nanmean(temp,1), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,[], [])
h2 = plot(nanmean(temp,1), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,[], [])
h3 = plot(nanmean(temp,1), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,[], [])
h4 = plot(nanmean(temp,1), 'm', 'LineWidth', 2)
clear temp

xlim([0 1700])
ylim([b a])
title ('ofc')
xlabel('time from stim onset')
ylabel('z-score power')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')

% cingulate
subplot(3,2,5)
hold on
chan_idx = cingulate_chan_idx_new;



temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,[], [])
h1 = plot(nanmean(temp,1), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,[], [])
h2 = plot(nanmean(temp,1), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,[], [])
h3 = plot(nanmean(temp,1), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,[], [])
h4 = plot(nanmean(temp,1), 'm', 'LineWidth', 2)
clear temp



xlim([0 1700])
ylim([b a])
title ('cing')
xlabel('time from stim onset')
ylabel('z-score power')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')

% insula
subplot(3,2,6)
hold on
chan_idx = insula_chan_idx_new;


temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,[], [])
h1 = plot(nanmean(temp,1), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,[], [])
h2 = plot(nanmean(temp,1), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,[], [])
h3 = plot(nanmean(temp,1), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,[], [])
h4 = plot(nanmean(temp,1), 'm', 'LineWidth', 2)
clear temp

xlim([0 1700])
ylim([b a])
title ('insula')
xlabel('time from stim onset')
ylabel('z-score power')
legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')




% gamma
% mtl

figure
hold on
win= .1*fs


chan_idx = MTL_chan_idx_new;
temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,freq_name,win)
h1 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,freq_name,win)
h2 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,freq_name,win)
h3 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,freq_name,win)
h4 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'm', 'LineWidth', 2)
clear temp

xlim([0 1700])
% ylim([b a])
title ('mtl')
xlabel('time from stim onset')
ylabel('z-score power')
legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')
set(gca, 'FontSize', 12, 'FontWeight', 'bold')


% temporal lobe
figure
hold on

chan_idx = temp_chan_idx_new;

temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,freq_name,win)
h1 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,freq_name,win)
h2 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,freq_name,win)
h3 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,freq_name,win)
h4 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'm', 'LineWidth', 2)
clear temp


xlim([0 1700])
%ylim([b a])
title ('temp lobe')
xlabel('time from stim onset')
ylabel('z-score power')
legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')

%
figure
hold on
chan_idx = fro_chan_idx_new;

temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,freq_name,win)
h1 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,freq_name,win)
h2 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,freq_name,win)
h3 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,freq_name,win)
h4 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'm', 'LineWidth', 2)
clear temp

xlim([0 1700])
%ylim([b a])
legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')
title ('fro')
xlabel('time from stim onset')
ylabel('z-score power')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')


% OFC_chan_idx_new
figure
hold on

chan_idx = OFC_chan_idx_new;

temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,freq_name,win)
h1 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,freq_name,win)
h2 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,freq_name,win)
h3 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,freq_name,win)
h4 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'm', 'LineWidth', 2)
clear temp

xlim([0 1700])
%ylim([b a])
title ('ofc')
legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')
xlabel('time from stim onset')
ylabel('z-score power')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')

% cingulate
figure
hold on
chan_idx = cingulate_chan_idx_new;


temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,freq_name,win)
h1 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,freq_name,win)
h2 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,freq_name,win)
h3 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,freq_name,win)
h4 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'm', 'LineWidth', 2)
clear temp


xlim([0 1700])
%ylim([b a])
title ('cing')
xlabel('time from stim onset')
ylabel('z-score power')
legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')

set(gca, 'FontSize', 9, 'FontWeight', 'bold')


% insula
figure
hold on
chan_idx = insula_chan_idx_new;

temp = power_series_mn_cond1(:,chan_idx)';
stdshade(temp,.1,'b',[],[] ,freq_name,win)
h1 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'b', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond2(:,chan_idx)';
stdshade(temp,.1,'g',[],[] ,freq_name,win)
h2 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'g', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond3(:,chan_idx)';
stdshade(temp,.1,'r',[],[] ,freq_name,win)
h3 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'r', 'LineWidth', 2)
clear temp

temp = power_series_mn_cond4(:,chan_idx)';
stdshade(temp,.1,'m',[],[] ,freq_name,win)
h4 = plot(conv(nanmean(temp,1),ones(1,win)/win,'same'), 'm', 'LineWidth', 2)
clear temp

xlim([0 1700])
%ylim([b a])
title ('insula')
xlabel('time from stim onset')
legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')
ylabel('z-score power')
legend([h1 h2 h3 h4 ], 'rep','incorr lure', 'corr lure', 'new')
set(gca, 'FontSize', 9, 'FontWeight', 'bold')





%% save timecourses
chan_idx = MTL_chan_idx_new;
%chan_idx  = MTL_chan_idx;


clear cond1 cond3 cond2 cond4
cond1    = nanmean(power_series_mn_cond1(:,chan_idx),2)';
cond2    = nanmean(power_series_mn_cond2(:,chan_idx),2)';
cond3    = nanmean(power_series_mn_cond3(:,chan_idx),2)';
cond4    = nanmean(power_series_mn_cond4(:,chan_idx),2)';

all_conds = [cond1; cond2; cond3;cond4]
cd('/mnt/yassamri/iEEG/sandra/group_data')

save(['MTL_gamma_' num2str(strt_time) '_' num2str(end_time) '_subj_' subj], 'all_conds')
%save(['MTL_gamma_allchans_subj_' subj], 'all_conds')





%%

chan_idx = temp_chan_idx_new;
%chan_idx = temp_chan_idx;
clear cond1 cond3 cond2 cond4
cond1    = nanmean(power_series_mn_cond1(:,chan_idx),2)';
cond2    = nanmean(power_series_mn_cond2(:,chan_idx),2)';
cond3    = nanmean(power_series_mn_cond3(:,chan_idx),2)';
cond4    = nanmean(power_series_mn_cond4(:,chan_idx),2)';

all_conds = [cond1; cond2; cond3;cond4]
cd('/mnt/yassamri/iEEG/sandra/group_data')



save(['temporal_lobe_gamma_' num2str(strt_time) '_' num2str(end_time) '_subj_' subj], 'all_conds')
%save(['temporal_lobe_gamma_allchans_subj_' subj], 'all_conds')



%%
chan_idx = fro_chan_idx_new;
%chan_idx = fro_chan_idx;

clear cond1 cond3 cond2 cond4
cond1    = nanmean(power_series_mn_cond1(:,chan_idx),2)';
cond2    = nanmean(power_series_mn_cond2(:,chan_idx),2)';
cond3    = nanmean(power_series_mn_cond3(:,chan_idx),2)';
cond4    = nanmean(power_series_mn_cond4(:,chan_idx),2)';

all_conds = [cond1; cond2; cond3;cond4]
cd('/mnt/yassamri/iEEG/sandra/group_data')

%  save(['frontal_lobe_gamma_allchans_subj_' subj], 'all_conds')
save(['frontal_lobe_gamma_' num2str(strt_time) '_' num2str(end_time) '_subj_' subj], 'all_conds')


%%
%save(['MTL_theta_' num2str(strt_time) '_' num2str(end_time) '_subj_' subj], 'all_conds')
%save(['MTL_theta_allchans_subj_' subj], 'all_conds')


%save(['temporal_lobe_theta_' num2str(strt_time) '_' num2str(end_time) '_subj_' subj], 'all_conds')
%save(['temporal_lobe_theta_allchans_subj_' subj], 'all_conds')


%save(['frontal_lobe_theta_allchans_subj_' subj], 'all_conds')
% save(['frontal_lobe_theta_' num2str(strt_time) '_' num2str(end_time) '_subj_' subj], 'all_conds')
