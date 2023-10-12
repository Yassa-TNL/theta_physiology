% output:
%1. data :notch filtered LFP data for elecs of interest
%2. trial_data :epoched data based on photodiode 

clear all; close all; clc
subj = '84'
addpath(genpath('/media/SSanDra/Pattern_Separation/analysis_pipeline_final/')) % FT folder
ft_defaults
script_path = '/media/SSanDra/Pattern_Separation/analysis_pipeline_final/fieldtrip-20181130/fileio/private' %path for functions to open files
file_path   = '/media/SSanDra/Pattern_Separation/subj_84/dataset' % data folder
cd(script_path)

%% get photodiode
photodiode   = read_neuralynx_ncs([file_path '/photo1.ncs']);
fs           = unique(photodiode.SampFreq);
photodiode_signal = -1*(photodiode.dat(:))';
clear photodiode

% downsample photodiode
D = 8;
photodiode_signal = decimate(photodiode_signal, D);
fs = fs/D % new fs
figure; hold on
plot(photodiode_signal)

%% manually indicate start point and end point using data cursor and previous fig
if strcmp ('84', subj)
    start.DataIndex = 3154321;
end
photodiode_signal = photodiode_signal(start.DataIndex : end);
%%
thresh_mltpl = 7
a            = diff(photodiode_signal);
meana        = nanmean(a);
stda         = nanstd(a);
onset        = (a>(meana+(thresh_mltpl*stda)));
offset       = (a<(meana-(thresh_mltpl*stda)));

% measure offset
off_idx      = (diff(offset)==1);
offsets      = find(off_idx);
offsets_dur  = diff(offsets);

real_offsets = offsets(offsets_dur>.5*mean(offsets_dur))
offset_idx   = nan(1,length(photodiode_signal));
offset_idx(real_offsets) = 1;

% measure onset
on_idx = (diff(onset)==1);
onsets = find(on_idx);
onsets_dur  = diff(onsets);
real_onsets = onsets(onsets_dur>.5*mean(onsets_dur));
onset_idx   = nan(1,length(photodiode_signal));
onset_idx(real_onsets) = 1;

% change these variables to indicate the correct start and end onsets and
% offsets
if subj == '84'
    % offest
    a = 3
    b = 0
    
    % onset
    c = 2
    d = 0
end

real_offsets = real_offsets(a:end-b);
real_onsets  = real_onsets(c:end-d);

onset_idx  = nan(1,length(photodiode_signal));
offset_idx = nan(1,length(photodiode_signal));
offset_idx(real_offsets) = 1;
onset_idx(real_onsets)   = 1;
if onsets_dur(end)> .5*mean(onsets_dur)
   onset_idx(onsets(end))   = 1;
   offset_idx(offsets(end)) = 1;
end

figure; hold on
plot(photodiode_signal)
plot(onset_idx-3500, 'r*')
plot(offset_idx-3500, 'g*')

% variable to  get trials
on_idx  = find(onset_idx==1);
off_idx = find(offset_idx==1);


on_idx_check  = on_idx+(start.DataIndex-1);
off_idx_check = off_idx+(start.DataIndex-1);

% check identified onsets out of total
photodiode   = read_neuralynx_ncs([file_path '/photo1.ncs']);
photodiode_signal = -1*(photodiode.dat(:))';
D = 8;
photodiode_signal = decimate(photodiode_signal, D);
onset_idx  = nan(1,length(photodiode_signal));
offset_idx = nan(1,length(photodiode_signal));
onset_idx(on_idx_check)  = 1;
offset_idx(off_idx_check) = 1;

figure; hold on
plot(photodiode_signal)
plot(onset_idx-3500, 'r*')
plot(offset_idx-3500, 'g*')

% timestamps are on_idx off_idx

%% get LFP data 
chans_per_strip = 8;
LFP_path      =  [file_path '/LFP'];
desired_chans = {'LHH', 'RHH'};
chan_counter  = length(desired_chans)*chans_per_strip;
lfp      = read_neuralynx_ncs([LFP_path '/' desired_chans{chan_label} num2str(a) '.ncs']);
temp_data = lfp.dat(:);
temp_LFP  = decimate(temp_data,D);
temp_LFP  = temp_LFP(start.DataIndex:end);
LFP_data  = zeros(chan_counter,length(temp_LFP));
fs           = unique(lfp.SampFreq)/D;

for chan_label = 1:length(desired_chans)
    for a = 1:chans_per_strip
        lfp= read_neuralynx_ncs([LFP_path '/LHH' num2str(a) '.ncs']);
        temp_data = lfp.dat(:);
        temp_LFP  = decimate(temp_data,D);
        LFP_data(((chan_label-1)*chans_per_strip) + a,:) = temp_LFP(start.DataIndex:end);
    end
end
rmpath(genpath('/media/SSanDra/Pattern_Separation/analysis_pipeline_final/fieldtrip-20181130'))

%% generate notch filter
flt1 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);
flt2 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
    'DesignMethod','butter','SampleRate',fs);
flt3 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',179,'HalfPowerFrequency2',181, ...
    'DesignMethod','butter','SampleRate',fs);

%% notch filter and preview plot 
figure;
hold on
yscale = 800
data = nan(size(LFP_data));

for chan = 1:chan_counter
    data (chan,:) = filtfilt(flt3, filtfilt(flt2, filtfilt(flt1,LFP_data(chan,:))));
    plot(0:1/fs:(size(LFP_data,2)-1)/fs, data(chan,:)+(yscale*chan))
end

%% get trial data
pre_stim  = .5;
post_stim = 2;
strt_time = on_idx' - pre_stim*fs;
end_time  = on_idx' + post_stim*fs;

trial_data     = single(nan(length(on_idx), length(strt_time(1):end_time(1)), size(data,1))); % trail x time x chan
for chan = 1:chan_counter
  for trial = 1:length(on_idx)
  trial_data(trial,:,chan) = data(chan, strt_time(trial):end_time(trial));
  end
end



