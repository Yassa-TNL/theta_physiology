close all;clc
%freq_ranges   = [2.9 6]; freq_name = 'low theta'% low theta
freq_ranges   = [40 200]; freq_name = 'gamma'

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
    OFC_chan_idx       = [1 2 58]
    
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

% remove prestim mean
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
        std_cond1(:,chan) = nanstd(mn_acrs_trials_cond1 (counter_freq_range,stim_onset:end,chan),0,1);
                                          
        
        power_series_mn_cond2 (:,chan) = nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset:end,chan),1)-...
                                    nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,1:stim_onset,chan),1));
       std_cond2(:,chan) = nanstd(mn_acrs_trials_cond2 (counter_freq_range,stim_onset:end,chan),0,1);

        
        power_series_mn_cond3 (:,chan) = nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset:end,chan),1)-...
                                         nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,1:stim_onset,chan),1));
        std_cond3(:,chan) = nanstd(mn_acrs_trials_cond3 (counter_freq_range,stim_onset:end,chan),0,1);

        
        power_series_mn_cond4 (:,chan) = nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset:end,chan),1)-...
                                        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,1:stim_onset,chan),1));
        std_cond4(:,chan) = nanstd(mn_acrs_trials_cond4 (counter_freq_range,stim_onset:end,chan),0,1);
        
    end
end

% get t-ratio for all chans
t_ratio_cond1=[]
t_ratio_cond2=[]
t_ratio_cond3=[]
t_ratio_cond4=[]
for chan = 1:chan_counter
t_ratio_cond1(chan,:)= power_series_mn_cond1 (:,chan)./ std_cond1(:,chan);
t_ratio_cond2(chan,:)= power_series_mn_cond2 (:,chan)./ std_cond1(:,chan);
t_ratio_cond3(chan,:)= power_series_mn_cond3 (:,chan)./ std_cond1(:,chan);
t_ratio_cond4(chan,:)= power_series_mn_cond4 (:,chan)./ std_cond1(:,chan);
end

% MTL
wind_counter = 50:25:975
end_idx      = 200; % moving win of 200ms w/ 75% overlap
chan_mgntd_win = zeros(length(MTL_chan_idx), 2);
counter=0;
for chan = MTL_chan_idx
    counter =counter+1
    specif_mgntd = [];
    win_spec     = [];
    
    % loop thru all win
   for win_idx = wind_counter 
    temp = [nanmean(t_ratio_cond1(chan,win_idx:win_idx+end_idx))...
            nanmean(t_ratio_cond2(chan,win_idx:win_idx+end_idx))...
            nanmean(t_ratio_cond3(chan,win_idx:win_idx+end_idx))...
            nanmean(t_ratio_cond4(chan,win_idx:win_idx+end_idx))];
        
        % if chan is cond specific, measure magnitude of specificity
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        specif_mgntd = [ specif_mgntd  mean([temp(3)-temp(1) temp(3)-temp(2) temp(3)-temp(4)])];
        win_spec              = [win_spec win_idx];
    end
   end
   
   % find max win
   if ~isempty([max(specif_mgntd) find(specif_mgntd==max(specif_mgntd))])
       chan_mgntd_win(counter,:) = [max(specif_mgntd) find(specif_mgntd==max(specif_mgntd))];
   end
end
figure
chan_mgntd_win =chan_mgntd_win(find(chan_mgntd_win(:,2)),2)
hist(wind_counter(chan_mgntd_win))
title(['subj ' subj  ' - '  'MTL - ' freq_name])
% xlim([0 1000])
 ylim([0 14])
set(gca, 'FontSize', 14, 'FontWeight', 'bold')



% Temporal Lobe
wind_counter = 50:25:975
end_idx      = 200; % moving win of 200ms w/ 75% overlap
chan_mgntd_win = zeros(length(temp_chan_idx), 2);
counter=0;
for chan = temp_chan_idx
    counter =counter+1
    specif_mgntd = [];
    win_spec     = [];
    
    % loop thru all win
   for win_idx = wind_counter 
    temp = [nanmean(t_ratio_cond1(chan,win_idx:win_idx+end_idx))...
            nanmean(t_ratio_cond2(chan,win_idx:win_idx+end_idx))...
            nanmean(t_ratio_cond3(chan,win_idx:win_idx+end_idx))...
            nanmean(t_ratio_cond4(chan,win_idx:win_idx+end_idx))];
        
        % if chan is cond specific, measure magnitude of specificity
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        specif_mgntd = [ specif_mgntd  mean([temp(3)-temp(1) temp(3)-temp(2) temp(3)-temp(4)])];
        win_spec              = [win_spec win_idx];
    end
   end
   
   % find max win
   if ~isempty([max(specif_mgntd) find(specif_mgntd==max(specif_mgntd))])
       chan_mgntd_win(counter,:) = [max(specif_mgntd) find(specif_mgntd==max(specif_mgntd))];
   end
end
chan_mgntd_win =chan_mgntd_win(find(chan_mgntd_win(:,2)),2)
figure
hist(wind_counter(chan_mgntd_win))
% xlim([0 1000])
 ylim([0 14])
title(['subj ' subj  ' - '  'temporal lobe - ' freq_name])
set(gca, 'FontSize', 14, 'FontWeight', 'bold')


% frontal
wind_counter = 50:25:975
end_idx      = 200; % moving win of 200ms w/ 75% overlap
chan_mgntd_win = zeros(length(fro_chan_idx), 2);
counter=0;
for chan = fro_chan_idx
    counter =counter+1
    specif_mgntd = [];
    win_spec     = [];
    
    % loop thru all win
   for win_idx = wind_counter 
    temp = [nanmean(t_ratio_cond1(chan,win_idx:win_idx+end_idx))...
            nanmean(t_ratio_cond2(chan,win_idx:win_idx+end_idx))...
            nanmean(t_ratio_cond3(chan,win_idx:win_idx+end_idx))...
            nanmean(t_ratio_cond4(chan,win_idx:win_idx+end_idx))];
        
        % if chan is cond specific, measure magnitude of specificity
    if sum(temp(3) > [temp(1) temp(2) temp(4)])> 2 % if condition spec
        specif_mgntd = [ specif_mgntd  mean([temp(3)-temp(1) temp(3)-temp(2) temp(3)-temp(4)])];
        win_spec              = [win_spec win_idx];
    end
   end
   
   % find max win
   if ~isempty([max(specif_mgntd) find(specif_mgntd==max(specif_mgntd))])
       chan_mgntd_win(counter,:) = [max(specif_mgntd) find(specif_mgntd==max(specif_mgntd))];
   end
end
chan_mgntd_win =chan_mgntd_win(find(chan_mgntd_win(:,2)),2)
figure
hist(wind_counter(chan_mgntd_win))
% xlim([0 1000])
 ylim([0 14])
title(['subj ' subj  ' - '  'frontal lobe - ' freq_name])
set(gca, 'FontSize', 14, 'FontWeight', 'bold')














