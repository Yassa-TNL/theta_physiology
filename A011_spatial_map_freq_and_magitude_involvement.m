% first run A003_condition_wise_spectral_anaylsis
save('elec_loc_wrksp', '-v7.3')
%%
%load('elec_loc_wrksp.mat')
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
load('freq_ranges.mat')

% calculate median freq
all_poss_unique_freq = [];
for f = 1:size(freq_ranges,1)
    all_poss_unique_freq (f,1) = mean(freq(freq>freq_ranges(f,1) & freq<freq_ranges(f,2)));
end

stim_onset  = 301;
strt_time   = 100;
end_time    = 1000;

freq_discrim_all_chan = zeros(size(freq_ranges,1), cond_num, chan_counter);

% gather cond specfici freq averages for each chan
for chan = 1:chan_counter
    for cntr = 1:size(freq_ranges,1)
        counter_freq_range = freq>freq_ranges(cntr,1) & freq<freq_ranges(cntr,2);
        freq_discrim_all_chan(cntr,:,chan) = [nanmean(nanmean(mn_acrs_trials_cond1(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))];
    end
end



% find chan peak freq
chan_peak_freq_disc_magnitude = zeros(chan_counter,2);
for chan = 1:chan_counter
    diff_freq_power = zeros(size(freq_ranges,1),1);
    for cntr = 1:size(freq_ranges,1)
        pair_1_diff = freq_discrim_all_chan(cntr, 3,chan)-freq_discrim_all_chan(cntr, 1,chan);
        pair_2_diff = freq_discrim_all_chan(cntr, 3,chan)-freq_discrim_all_chan(cntr, 2,chan);
        pair_3_diff = freq_discrim_all_chan(cntr, 3,chan)-freq_discrim_all_chan(cntr, 4,chan);
        
        if [pair_1_diff pair_2_diff pair_3_diff]>0
            diff_freq_power(cntr,1) = mean([pair_1_diff pair_2_diff pair_3_diff]);
        end

        clear pair_1_diff pair_2_diff pair_3_diff
    end
    
    if any(diff_freq_power) % if there is a freq that only discriminates PS cond
    freq_idx = find(diff_freq_power==max(diff_freq_power ));
    peak_mean_freq = freq_ranges(freq_idx,:);   
    counter_freq_range = mean(freq(freq>peak_mean_freq(1) & freq<peak_mean_freq(2)));
   
    % get peak freq and disc value
    chan_peak_freq_disc_magnitude(chan,:) = [ counter_freq_range  max(diff_freq_power) ]; 
   
    end 
end

% normal disc magnitude by max diameter from all time win
cd(['/mnt/yassamri/iEEG/sandra/subj_' num2str(subj)])
load('max_diameter.mat')

temp = zeros(size(chan_peak_freq_disc_magnitude,1),1)
for a = 1:size(chan_peak_freq_disc_magnitude,1)
    temp(a) = chan_peak_freq_disc_magnitude(a,2)/max_diameter
end

chan_peak_freq_disc_magnitude(:,2) = temp;



%% store diameters across win
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
load('freq_ranges.mat')

% calculate median freq
all_poss_unique_freq = [];
for f = 1:size(freq_ranges,1)
    all_poss_unique_freq (f,1) = mean(freq(freq>freq_ranges(f,1) & freq<freq_ranges(f,2)));
end

stim_onset  = 301;
wind_counter = 0:50:1700
chan_peak_freq_disc_magnitude = zeros(chan_counter,2, length(wind_counter));

for time_idx = 1:length(wind_counter)

freq_discrim_all_chan = zeros(size(freq_ranges,1), cond_num, chan_counter);
strt_time             = wind_counter(time_idx);
end_time              = strt_time+100;

% gather cond specfici freq averages for each chan
for chan = 1:chan_counter
    for cntr = 1:size(freq_ranges,1)
        counter_freq_range = freq>freq_ranges(cntr,1) & freq<freq_ranges(cntr,2);
        freq_discrim_all_chan(cntr,:,chan) = [nanmean(nanmean(mn_acrs_trials_cond1(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond2(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond3(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))...
        nanmean(nanmean(mn_acrs_trials_cond4(counter_freq_range,stim_onset+strt_time:stim_onset+end_time,chan),1))];
    end
end



% find chan peak freq
for chan = 1:chan_counter
    diff_freq_power = zeros(size(freq_ranges,1),1);
    for cntr = 1:size(freq_ranges,1)
        pair_1_diff = freq_discrim_all_chan(cntr, 3,chan)-freq_discrim_all_chan(cntr, 1,chan);
        pair_2_diff = freq_discrim_all_chan(cntr, 3,chan)-freq_discrim_all_chan(cntr, 2,chan);
        pair_3_diff = freq_discrim_all_chan(cntr, 3,chan)-freq_discrim_all_chan(cntr, 4,chan);
        
        if [pair_1_diff pair_2_diff pair_3_diff]>0
            diff_freq_power(cntr,1) = mean([pair_1_diff pair_2_diff pair_3_diff]);
        end

        clear pair_1_diff pair_2_diff pair_3_diff
    end
    
    if any(diff_freq_power) % if there is a freq that only discriminates PS cond
    
    freq_idx           = find(diff_freq_power==max(diff_freq_power));
    peak_mean_freq     = freq_ranges(freq_idx,:);   
    counter_freq_range = mean(freq(freq>peak_mean_freq(1) & freq<peak_mean_freq(2)));
   
    % get peak freq and disc value
    chan_peak_freq_disc_magnitude(chan,:,time_idx) = [ counter_freq_range  max(diff_freq_power) ]; 
   
    end 
end

end


% get max
max_diameter = max(max(squeeze(chan_peak_freq_disc_magnitude(:,2,:))));
save('max_diameter','max_diameter')