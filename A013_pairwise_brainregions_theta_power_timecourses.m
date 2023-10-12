% organize chan labels to create excel column
for chan = 1:length(chan_label)
    chan_label_new{chan} = chan_label{chan}(4:end-3)
end
m = chan_label_new'
open m

%%
freq_ranges   = [2.9 6]  % low theta
% freq_ranges = [40 200] % gamma

% elecs
if strcmp('57', subj)
    fro_chan_idx  = [1 4:8 13 16:17 25:27 35 37 44 46:47 70:83 89 91 98 99 101 102 108 109 110 112];
    MTL_chan_idx  = [48:50 58:61 67:68 113:115 122:124 132:134];
    temp_chan_idx = [55 56 65 66 72 73 118:121 128:131 136:140];
    insula_chan_idx = [29:32 38:40 93:96 103 104]
    cingulate_chan_idx = [2 9:11 74:77 ]
    OFC_chan_idx = [18:19 84:86]
elseif strcmp('39', subj)
    fro_chan_idx  = [8 16:19 25:26 28 63:64 66 72:74];
    MTL_chan_idx  = [29:32 39:42 51:52 77:79 85:88 94:96];
    temp_chan_idx = [34:38 45:48 54:57 83:84 91:93 100:103 ];
    insula_chan_idx = [20 21]
    cingulate_chan_idx = [11 12]
    OFC_chan_idx = [1 2 58]
    
elseif strcmp('66', subj)
    fro_chan_idx  = [37:40 49:50 95:97 101 108 109 110];
    MTL_chan_idx  = [1:5 11:15 21:23 61:65 71:75 81:83 ];
    temp_chan_idx = [7:9 16:20 24 25 29 30 48 68:70 77 78 84:90 ];
    insula_chan_idx = [];
    cingulate_chan_idx = [31:33 35 51:54 91:93 111:114];
    OFC_chan_idx = [41 42 102];
    
elseif strcmp('63', subj)
    fro_chan_idx  = [45 52:56 116 ];
    MTL_chan_idx  = [5:8 17:19 25:28 67:70 77:80 91:94];
    temp_chan_idx = [14 22 23 24 31 34 73:76 84:86];
    insula_chan_idx = [];
    cingulate_chan_idx = [47 49 97:101];
    OFC_chan_idx = [37:39];
    
elseif strcmp('44', subj)
    fro_chan_idx  = [36:38 44:47 84 :86];
    MTL_chan_idx  = [1:4 11:14 21 22 23 49 58:62 68 70:72];
    temp_chan_idx = [7:10 16:19 25:28 54:57 64:67 74:76];
    insula_chan_idx = [];
    cingulate_chan_idx = [29:32];
    OFC_chan_idx = [39 40 77:80];
    
elseif strcmp('83', subj)
    fro_chan_idx  = [7:10 12 15:20 25:28 32:35 42 44 50 52 29 ];
    MTL_chan_idx  = [];
    temp_chan_idx = [5 6 21 47];
    insula_chan_idx    = [];
    cingulate_chan_idx = [];
    OFC_chan_idx       = [];
    
elseif strcmp('84', subj)
    fro_chan_idx  = [ 27 29 31 32  28 30 40];
    MTL_chan_idx  = [17:20 49:52];
    temp_chan_idx = [11:13 15:16 41:44 47 48 55 56];
    insula_chan_idx    = [];
    cingulate_chan_idx = [];
    OFC_chan_idx       = [1 2 3 33 34 36];
end

total_chans = length(fro_chan_idx)+length(MTL_chan_idx)+length(temp_chan_idx);

stim_onset  = 301;

% theta
for chan = 1:chan_counter
    for cntr = 1:size(freq_ranges,1)
        counter_freq_range = freq>freq_ranges(cntr,1) & freq<freq_ranges(cntr,2);
        power_series_mn_cond3 (:,chan) = nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset:end,chan),1);
        power_series_std_cond3(:,chan) = nanstd(mn_acrs_trials_cond3(counter_freq_range,stim_onset:end,chan),0,1);
    end
end

% find selective theta chans
strt_time = 200
end_time  = 400 
fro_chan_idx_new=[]
for chan =fro_chan_idx
    temp = [nanmean(nanmean(mn_acrs_trials_cond1(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))]
    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        fro_chan_idx_new = [fro_chan_idx_new chan]
    end
end

MTL_chan_idx_new=[]
for chan =MTL_chan_idx
    temp = [nanmean(nanmean(mn_acrs_trials_cond1(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))]
    sum(temp(3) > [temp(1) temp(2) temp(4)])
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        MTL_chan_idx_new = [MTL_chan_idx_new chan]
    end
end

temp_chan_idx_new=[]
for chan =temp_chan_idx
    temp = [nanmean(nanmean(mn_acrs_trials_cond1(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))]
    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        temp_chan_idx_new = [temp_chan_idx_new chan]
    end
end


insula_chan_idx_new=[]
for chan =insula_chan_idx
    temp = [nanmean(nanmean(mn_acrs_trials_cond1(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))]
    
    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        insula_chan_idx_new = [insula_chan_idx_new chan]
    end
    clear temp
end

cingulate_chan_idx_new=[]
for chan =cingulate_chan_idx
    temp = [nanmean(nanmean(mn_acrs_trials_cond1(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))]
    sum(temp(3) > [temp(1) temp(2) temp(4)])
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
       cingulate_chan_idx_new = [cingulate_chan_idx_new chan]
    end
end

OFC_chan_idx_new=[]
for chan =OFC_chan_idx
    temp = [nanmean(nanmean(mn_acrs_trials_cond1(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))]
    sum(temp(3) > [temp(1) temp(2) temp(4)]);
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        OFC_chan_idx_new = [OFC_chan_idx_new chan]
    end
end
   
% plot delta - selective chans

figure
hold on
temp = nanmean(power_series_mn_cond3(:,fro_chan_idx_new),2)
zerod_sgn = temp - temp(1)
h1 = plot(zerod_sgn, 'k', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,temp_chan_idx_new),2)
zerod_sgn = temp - temp(1)
h2= plot(zerod_sgn, 'b', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,MTL_chan_idx_new),2)
zerod_sgn = temp - temp(1)
h3=plot(zerod_sgn, 'r', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,insula_chan_idx_new),2)
zerod_sgn = temp - temp(1)
h4=plot(zerod_sgn, 'm', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,cingulate_chan_idx_new),2)
zerod_sgn = temp - temp(1)
h5=plot(zerod_sgn, 'g', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,OFC_chan_idx_new),2)
zerod_sgn = temp - temp(1)
h6=plot(zerod_sgn, 'y', 'LineWidth', 2)


xlim([0 800])
%ylim([-.3 .6])
title ([num2str(freq_ranges(1)) '-' num2str(freq_ranges(2)) ' hz power'])
xlabel('time from stim onset')
ylabel('z-score power')
legend([h2 h1 h4 h5 h3 h6], 'TEMP','FRO', 'INS',  'CING', 'MTL', 'OFC')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

%% plot gamma selective chans

win= .2*fs
figure
hold on
temp = conv(nanmean(power_series_mn_cond3(:,fro_chan_idx_new),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h1 = plot(conv(zerod_sgn, ones(1,win)/win,'same'), 'k', 'LineWidth', 2)

temp = conv(nanmean(power_series_mn_cond3(:,temp_chan_idx_new),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h2= plot(zerod_sgn, 'b', 'LineWidth', 2)

temp = conv(nanmean(power_series_mn_cond3(:,MTL_chan_idx_new),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h3= plot(zerod_sgn, 'r', 'LineWidth', 2)


temp = conv(nanmean(power_series_mn_cond3(:,insula_chan_idx_new),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h4= plot(zerod_sgn, 'm', 'LineWidth', 2)

temp = conv(nanmean(power_series_mn_cond3(:,cingulate_chan_idx_new),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h5= plot(zerod_sgn, 'g', 'LineWidth', 2)

temp = conv(nanmean(power_series_mn_cond3(:,OFC_chan_idx_new),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h6= plot(zerod_sgn, 'y', 'LineWidth', 2)


xlim([0 800])
%ylim([-.3 .6])
title ([num2str(freq_ranges(1)) '-' num2str(freq_ranges(2)) ' hz power'])
xlabel('time from stim onset')
ylabel('z-score power')
%legend([h2 h1 h3 ], 'TEMP','FRO', 'MTL')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

legend([h2 h1 h4 h5 h3 h6], 'TEMP','FRO', 'INS',  'CING', 'MTL', 'OFC')

%% plot delta - average across channels

figure
hold on
temp = nanmean(power_series_mn_cond3(:,fro_chan_idx),2)
zerod_sgn = temp - temp(1)
h1 = plot(zerod_sgn, 'k', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,temp_chan_idx),2)
zerod_sgn = temp - temp(1)
h2= plot(zerod_sgn, 'b', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,MTL_chan_idx),2)
zerod_sgn = temp - temp(1)
h3=plot(zerod_sgn, 'r', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,insula_chan_idx),2)
zerod_sgn = temp - temp(1)
h4=plot(zerod_sgn, 'm', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,cingulate_chan_idx),2)
zerod_sgn = temp - temp(1)
h5=plot(zerod_sgn, 'g', 'LineWidth', 2)

temp = nanmean(power_series_mn_cond3(:,OFC_chan_idx),2)
zerod_sgn = temp - temp(1)
h6=plot(zerod_sgn, 'y', 'LineWidth', 2)


xlim([0 2000])
%ylim([-.3 .6])
title ([num2str(freq_ranges(1)) '-' num2str(freq_ranges(2)) ' hz power'])
xlabel('time from stim onset')
ylabel('z-score power')
legend([h2 h1 h4 h5 h3 h6], 'TEMP','FRO', 'INS',  'CING', 'MTL', 'OFC')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

%%
win= .1*fs
figure
hold on
temp = conv(nanmean(power_series_mn_cond3(:,fro_chan_idx),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h1 = plot(conv(zerod_sgn, ones(1,win)/win,'same'), 'k', 'LineWidth', 2)

temp = conv(nanmean(power_series_mn_cond3(:,temp_chan_idx),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h2= plot(zerod_sgn, 'b', 'LineWidth', 2)

temp = conv(nanmean(power_series_mn_cond3(:,MTL_chan_idx),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h3= plot(zerod_sgn, 'r', 'LineWidth', 2)


temp = conv(nanmean(power_series_mn_cond3(:,insula_chan_idx),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h4= plot(zerod_sgn, 'm', 'LineWidth', 2)

temp = conv(nanmean(power_series_mn_cond3(:,cingulate_chan_idx),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h5= plot(zerod_sgn, 'g', 'LineWidth', 2)

temp = conv(nanmean(power_series_mn_cond3(:,OFC_chan_idx),2), ones(1,win)/win,'same')
zerod_sgn = temp - temp(1)
h6= plot(zerod_sgn, 'y', 'LineWidth', 2)


xlim([0 2000])
%ylim([-.3 .6])
title ([num2str(freq_ranges(1)) '-' num2str(freq_ranges(2)) ' hz power'])
xlabel('time from stim onset')
ylabel('z-score power')
%legend([h2 h1 h3 ], 'TEMP','FRO', 'MTL')
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

legend([h2 h1 h4 h5 h3 h6], 'TEMP','FRO', 'INS',  'CING', 'MTL', 'OFC')

%% plot delta - for each channel
a = .7
b = 1.4
figure
hold on
h1=plot(power_series_mn_cond3(:,fro_chan_idx(1)), 'k')
plot(power_series_mn_cond3(:,fro_chan_idx), 'k')
h2=plot(power_series_mn_cond3(:,temp_chan_idx(1))+a, 'b')
plot(power_series_mn_cond3(:,temp_chan_idx)+a, 'b')
h3=plot(power_series_mn_cond3(:,MTL_chan_idx(1))+b, 'r')
plot(power_series_mn_cond3(:,MTL_chan_idx)+b, 'r')
 
xlim([0 800])
title ([num2str(freq_ranges(1)) '-' num2str(freq_ranges(2)) ' hz power'])
xlabel('time from stim onset')
ylabel('z-score power')
legend([h1 h2 h3],'FC', 'TC', 'MTL')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')


%%
figure
hold on
a = .1
for chan = fro_chan_idx
    a = a+.2
   plot(power_series_mn_cond3(:,chan)+a, 'k') 
end
h1 = plot(power_series_mn_cond3(:,chan)+a, 'k') 
h2 = plot(nanmean(power_series_mn_cond3(:,fro_chan_idx),2), 'r')
xlim([0 800])
%ylim([-.3 8])

title ([num2str(freq_ranges(1)) '-' num2str(freq_ranges(2)) ' hz power'])
xlabel('time from stim onset')
ylabel('z-score power')
legend([h1 h2 ],'indiv chans', 'chan avrg' )
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
