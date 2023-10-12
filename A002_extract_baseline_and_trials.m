% Output:
% 1. chanel mean and std for later normalization
% 2. trial_data
tic
clear all; close all;clc
patient_num   = 1
subjects      = {'39', '57', '66', '63', '44', '83'};
lock          = 'onset'%'response'
ref           = 'LM'%[] %:'' 'WM' 'LM' 'ESR [for ir83]' : la placian montage' 'CAR'
par_num       = 3;
if patient_num==6
    ref           = 'ESR'
end
downsample    = 'yes';
select_chan   = 3; % 1: get specified chans indicated by chan_name 2: uses logical vector of chans of interest (HC) 3: gets all chans
chan_name_1   = 'HH'; %'AM' %'HH'
chan_name_2   = 'TH'; %'AC' %'TH'
normalization = 1; % 1 = entire task 2 = prestim
freq_analysis = 'wavelet';
norm = {'entire_recording' 'prestim' 'cond_spec_prestim'}
path_name = '/mnt/yassamri/iEEG/sandra'
addpath([path_name '/analysis_pipeline_final'])

% convert besa2edf
filenames    = {'2016042511_0003.edf' '2017032319_0024.edf' '2017121913_0027.edf' '2017092119_0028.edf' '2016073015_0001.edf' '2018101311_0002.edf'  };%'2016073015_0002.edf'
filename     = filenames{patient_num}
subj         = subjects {patient_num}

event_chan   = 1;
if ismember(patient_num,[2 3 4 6])
    event_chan  = 2;
end

% load edf file data
cd([path_name '/subj_' subj])
[hdr, record]  = edfread(filename);
chan_label     = hdr.label;
fs             = hdr.frequency(1);
events         = -1*record(event_chan,:);

% load trial info
load(['behavior_subj' subj])
nTrials_study = size(training_behav_matrix,1);
nTrials_test  = size(testing_behav_matrix,1);
nBlocks       = 2;

% get & plot event channel
events = -1*record(event_chan,:);

if strcmp('yes',downsample)
    D = 5;
else
    D = 1;
end
events = decimate(events, D);


% figure; hold on
% plot(events)
a = diff(events);

meana  = mean(a);
stda   = std(a);
onset  = (a>(meana+(10*stda)));
offset = (a<(meana-(10*stda)));
plot(onset*50000)
plot(offset*5000)

% measure offset
off_idx = (diff(offset)==1);
offsets = find(off_idx);
offsets_dur = diff(offsets);
real_offsets = offsets(offsets_dur > .5*mean(offsets_dur));
offset_idx = nan(1,length(events));
offset_idx(real_offsets) = 1;

if strcmp('39', subj)
    if offsets_dur(end)> .5*mean(offsets_dur)
        offset_idx(offsets(end)) = 1;
    end
elseif  strcmp('66', subj) | strcmp('57', subj)| strcmp('44', subj)
    offset_idx(offsets(end)) = 1;
    real_offsets = [real_offsets offsets(end)]
end

% measure onset
on_idx = (diff(onset)==1);
onsets = find(on_idx);
onsets_dur = diff(onsets);
real_onsets = onsets(onsets_dur > .5*mean(onsets_dur));
onset_idx = nan(1,length(events));
max_error = 0 %0.018; %.024sec, 24msec

onset_idx(real_onsets+ceil(max_error*(fs/D))) = 1;


% plot to check
figure; hold on
plot(events)
plot(onset_idx, 'r*')
plot(offset_idx, 'g*')

% remove extra offsets [a b] - onset [c d]
if strcmp('39', subj)
    % offest
    a = 2
    b = 0
    
    % onset
    c = 2
    d = 0
elseif strcmp('57', subj)
    %offset
    a = 1
    b = 1
    
    % onset
    c = 1
    d = 0
elseif strcmp('66', subj)
    %offset
    a = 1
    b = 0
    
    % onset
    c = 1
    d = 0
elseif  strcmp('44', subj)
    %offset
    a = 2
    b = 0
    
    % onset
    c = 1
    d = 1
    
elseif  strcmp('63', subj)
    %offset
    a = 2
    b = 0
    
    % onset
    c = 1
    d = 1
elseif strcmp('83', subj)
    %offset
    a = 4
    b = 0
    
    % onset
    c = 5
    d = 0
end

real_offsets = real_offsets(a:end-b);
real_onsets  = real_onsets (c:end-d);

% redefine indices
onset_idx  = nan(1,length(events));
offset_idx = nan(1,length(events));
offset_idx(real_offsets) = 1;
if strcmp('39', subj)
    if offsets_dur(end)> .5*mean(offsets_dur)
        offset_idx(offsets(end)) = 1;
    end
elseif  strcmp('66', subj) || strcmp('44', subj)||  strcmp('63', subj) %strcmp('83', subj) ||
    offset_idx(offsets(end)) = 1;
    real_offsets = [real_offsets offsets(end)]
end
onset_idx(real_onsets+(max_error*fs/D)) = 1;

% plot final detected timestamps
figure; hold on
plot(events)
plot(onset_idx-3000/D, 'r*')
plot(offset_idx-3000/D, 'g*')

% redefine w/ log vec
onset_idx  = zeros(1,length(events));
offset_idx = zeros(1,length(events));
offset_idx(real_offsets) = 1;
if strcmp('39', subj)
    if offsets_dur(end)> .5*mean(offsets_dur)
        offset_idx(offsets(end)) = 1;
    end
elseif  strcmp('66', subj) || strcmp('44', subj) %|| strcmp('83', subj)
    offset_idx(offsets(end)) = 1;
    real_offsets = [real_offsets offsets(end)]
end
onset_idx(real_onsets+(max_error*fs/D)) = 1;

% variable to  get trials
on_idx  = find(onset_idx);
off_idx = find(offset_idx);
difference = off_idx(:) - on_idx(:)

% selec sites of interest
% chans from strings
if select_chan == 1
    if strcmp('83', subj)
        chan_set_1 = regexp(chan_label, 'LRH');
        chan_set_2 = regexp(chan_label, 'LLH');
        
    else
        chan_set_1 = regexp(chan_label, chan_name_1);
        chan_set_2 = regexp(chan_label, chan_name_2);
        chan_idx   = zeros(length(chan_set_1),1);
    end
    for a =1:length(chan_set_1)
        if ~isempty(chan_set_1{a}) ||  ~isempty(chan_set_2{a})
            chan_idx(a) = 1;
        end
    end
    
    desired_chans    = logical(chan_idx);
    
    % verified HC chans
elseif select_chan == 2
    
    hippocampus_chans = zeros(1,length(chan_label));
    if strcmp('44',subj)
        hippocampus_chans([6:8  27:28 16:18 76:77 88]) = 1; % IR44
    elseif strcmp('39',subj)
        hippocampus_chans([45:48 57:58 159:160 169:170]) = 1; % IR39
    elseif strcmp('57',subj)
        hippocampus_chans([68:71 79:80 150 158:159 ]) = 1; %57
    end
    desired_chans= (logical(hippocampus_chans));
    
    % all chans in brain
elseif select_chan == 3
    if strcmp('39',subj)
        out_of_brain = [1:4 24 34 64:128 147 157 167 178];
    elseif strcmp('57',subj)
        out_of_brain = [1:4 13:14 24 65:67 77 85:97 117 147 168:184 ];
    elseif strcmp('66',subj)
        out_of_brain = [1:4 65 66 127:137];
    elseif strcmp('83',subj)
        out_of_brain = [];
    elseif strcmp('63',subj)
        out_of_brain = [];
    elseif strcmp('44',subj)
        out_of_brain = [1:4 24 34:36 57:63 73 84 104:135];
    end
    
    desired_chans = ones(1,length(chan_label));
    desired_chans(out_of_brain) = 0;
end


data = record(logical(desired_chans),:);
chan_label   = chan_label(logical(desired_chans));
chan_counter = length(chan_label);

% REREF: WM, LM, CAR re ref
if strcmp('WM',ref)
    [WM_ref_vector] = get_WM_reref(subj, chan_label);
    data_reref = data - data(WM_ref_vector,:);
    clear data
    data = data_reref;
    clear data_reref
elseif strcmp('LM', ref)
    [LM_ref_vector] = get_LM_reref(subj);
    data_reref = nan(size(data));
    
    for iProbe = 1:size(LM_ref_vector,1) % loop thru probe
        for iElec = LM_ref_vector(iProbe,1):LM_ref_vector(iProbe,2)
            if iElec == LM_ref_vector(iProbe,1) % if on first elec
                data_reref(iElec,:) = data(iElec,:) - data(iElec+1,:);
            elseif iElec == LM_ref_vector(iProbe,2) % if on last elec
                data_reref(iElec,:) = data(iElec,:) - data(iElec-1,:);
            else
                data_reref(iElec,:) = data(iElec,:) - (data(iElec-1,:)+data(iElec+1,:))/2;
            end
            disp(iElec)
        end
    end
    
    clear data
    data = data_reref;
elseif strcmp('CAR',ref)
    CAR_ref = nanmean(data,1) ;
    CAR_data = zeros(size(data));
    for chan = 1:size(data,1)
        CAR_data(chan,:) = data(chan,:)-CAR_ref;
    end
    clear data
    data = CAR_data;
    clear CAR_data
elseif strcmp('ESR',ref)
    data_reref = nan(size(data));
    if strcmp('83',subj)
        ref_vector = [5 65; 67 83; 84 99];
    end
    
    for iProbe =  1:size(ref_vector,1)
        for iElec = ref_vector(iProbe,1):ref_vector(iProbe,2)
            data_reref(iElec,:) =  data(iElec,:) - mean(data(ref_vector(iProbe,1):ref_vector(iProbe,2),:) ,1);
        end
    end
    clear data
    data = data_reref;
end

% generate notch filter
flt1 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',59,'HalfPowerFrequency2',61, ...
    'DesignMethod','butter','SampleRate',fs);
flt2 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',119,'HalfPowerFrequency2',121, ...
    'DesignMethod','butter','SampleRate',fs);
flt3 = designfilt('bandstopiir','FilterOrder',2, ...
    'HalfPowerFrequency1',179,'HalfPowerFrequency2',181, ...
    'DesignMethod','butter','SampleRate',fs);

% notch filter
data_new = nan(size(data));
for chan =1:chan_counter
    data_new (chan,:) = filtfilt(flt3, filtfilt(flt2, filtfilt(flt1,data(chan,:))));
end
clear data
data = data_new;


ds_data = nan(size(data,1), ceil(size(data,2)/D));
for chan = 1:size(data,1)
    ds_data(chan,:) = decimate(data(chan,:), D);
end
fs      = fs/D;
clear data
data = ds_data;

% down sample data and timestamps to fs = 500
if strcmp('yes', downsample)
    D2 = 2;
    LFP_data = zeros(size(data,1), length(decimate(data(1,:), D2)));
    for iElec = 1:size(data,1)
        LFP_data(iElec,:) = decimate(data(iElec,:), D2);
    end
    fs = fs/D2;
    %timestamp data ***** THESE ARE FINAL DS TS
    on_idx  = round(on_idx*1/D2);
    off_idx = round(off_idx*1/D2);
end
clear data
data = LFP_data;

% get artifact indices - chan
artifact_points      = zeros(size(data,1), size(data,2));
chan_artifact_thresh = zeros(chan_counter,2);
sec   = 1; % duration to remove around artifacts
if strcmp('44', subj)
    mltpl = 3;  %mutliple of std above mean for artifact thresholding 2.8 --> 3 for IR44
else
    mltpl = 4.5;
end
for chan =1:chan_counter
    x = data(chan,:);
    
    %get mean and std
    mn = mean(x);
    sd = std(x);
    pos_thresh = mn+(mltpl*sd);
    neg_thresh = mn-(mltpl*sd);
    idx_above_thresh = (x>pos_thresh);
    idx_below_thresh = (x<neg_thresh);
    
    %     figure;hold on; plot(0:1/fs:(size(data,2)-1)/fs, x)
    %     line([0 (size(data,2)-1)/fs] ,[pos_thresh pos_thresh], 'Color','k','LineWidth',4)
    %     line([0 (size(data,2)-1)/fs] ,[neg_thresh neg_thresh], 'Color','k','LineWidth',4)
    
    % get artifact indices/chan
    temp = idx_above_thresh+idx_below_thresh;
    artifact_idx = find(temp);
    for counter = 1:length(artifact_idx)
        if artifact_idx(counter)<=sec*fs
            temp(1:(artifact_idx(counter)+sec*fs)) = 1; %set sec bef and sec after = artifact
        elseif artifact_idx(counter)+sec*fs>length(temp)
            temp(artifact_idx(counter):end) = 1;
        else
            temp((artifact_idx(counter)-sec*fs):(artifact_idx(counter)+sec*fs)) = 1; %set sec bef and sec after = artifact
        end
    end
    artifact_points(chan,:) = temp';
    chan_artifact_thresh(chan,:) = [neg_thresh pos_thresh];
    time_domain_chan_mn_std(chan,:) = [mn sd];
    clear temp mn sd artifact_idx x pos_thresh neg_thresh idx_above_thresh idx_below_thresh
    
end

% wavelet params
dt = 1/fs;
NumVoices = 32;
a0 = 2^(1/NumVoices);
wavCenterFreq = 6/(2*pi);
minfreq  = 3;
maxfreq  = 200;
minscale = wavCenterFreq/(maxfreq*dt);
maxscale = wavCenterFreq/(minfreq*dt);
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale));
scales   = a0.^(minscale:maxscale).*dt;
freq     = wavCenterFreq./(fs*scales.*dt);

% get trial data
% get start and end trial indices
ISI                = .5;
if strcmp ('onset',lock)
    pre_stim  = .9;
    post_stim = 2.4;
    strt_time = on_idx' - pre_stim*fs;
    end_time  = on_idx' + post_stim*fs;
elseif strcmp ('onset2',lock)
    pre_stim  = .9;
    post_stim = 3;
    strt_time = on_idx' - pre_stim*fs;
    end_time  = on_idx' + post_stim*fs;
elseif strcmp ('offset',lock)
    pre_stim  = 1;
    post_stim = .5;
    strt_time = off_idx' - pre_stim*fs
    end_time  = off_idx' + post_stim*fs;
elseif strcmp('response',lock)
    % find response indices
    responses          = [training_behav_matrix(:,4); testing_behav_matrix(:,4)];
    resp_idx_tmp       = fs*([training_behav_matrix(:,4); testing_behav_matrix(:,4)]+ISI);
    test_trial_idx     = 1:nTrials_test+nTrials_study;
    resp_idx           = off_idx(test_trial_idx)+resp_idx_tmp';
    pre_stim           = 1.5;
    post_stim          = .5;
    strt_time          = resp_idx' - pre_stim*fs;
    end_time           = resp_idx' + post_stim*fs;
end
%%
trial_data     = single(nan(length(on_idx), length(strt_time(1):end_time(1)), size(data,1))); % trail x time x chan
for chan = 1:chan_counter
    for trial = 1:length(on_idx)
        trial_data(trial,:,chan) = data(chan, strt_time(trial):end_time(trial));
    end
end

cd ([path_name '/subj_' subj ])

if strcmp('83',subj)
    ref = 'LM'
end

if strcmp('onset',lock) || strcmp('onset2',lock)
    if strcmp('',downsample)
        save(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_' num2str(select_chan) '_NotDownsampled'],...
            'trial_data', 'fs', 'chan_counter', 'nTrials_study', 'nTrials_test', 'nBlocks', 'chan_label',...
            'pre_stim', 'post_stim',  'lock', 'on_idx' , 'off_idx', 'artifact_points', 'chan_artifact_thresh',...
            'chan_artifact_thresh', 'freq',  'minfreq', 'maxfreq', '-v7.3');
    elseif strcmp('yes',downsample)
        save(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_' num2str(select_chan) '_fs_' num2str(fs)],...
            'trial_data', 'fs', 'chan_counter', 'nTrials_study', 'nTrials_test', 'nBlocks', 'chan_label',...
            'pre_stim', 'post_stim', 'lock', 'on_idx' , 'off_idx', 'artifact_points', 'chan_artifact_thresh',...
            'chan_artifact_thresh', 'freq',  'minfreq', 'maxfreq');
    end
elseif strcmp('response',lock)
    
    if strcmp('',downsample)
        save(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_' num2str(select_chan) '_NotDownsampled'],...
            'trial_data', 'fs', 'chan_counter', 'nTrials_study', 'nTrials_test', 'nBlocks', 'chan_label',...
            'pre_stim', 'post_stim', 'lock', 'on_idx' , 'off_idx', 'artifact_points', 'chan_artifact_thresh',...
            'chan_artifact_thresh', 'freq',  'minfreq', 'maxfreq','responses', 'resp_idx', '-v7.3');
    elseif  strcmp('yes',downsample)
        save(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_' num2str(select_chan) '_fs_' num2str(fs)],...
            'trial_data', 'fs', 'chan_counter', 'nTrials_study', 'nTrials_test', 'nBlocks', 'chan_label',...
            'pre_stim', 'post_stim', 'lock', 'on_idx' , 'off_idx', 'artifact_points', 'chan_artifact_thresh',...
            'chan_artifact_thresh', 'freq',  'minfreq', 'maxfreq','responses', 'resp_idx');
    end
end

%% do spectrogram analysis instead
if strcmp('spectrogram',freq_analysis)
    data(logical(artifact_points)) = nan;
    win   = 1*fs %win = 1 sec
    F     = freq;
    ovrlp = .99*win
    [~,f,t,ps] = spectrogram(data(1, :),win,ovrlp,F,fs);
    baseline_spectrogram = nan(size(ps,1), size(ps,2), chan_counter);
    chan_powr_mn_spectrogram  = nan(size(baseline_spectrogram,1),chan_counter);
    chan_powr_std_spectrogram = nan(size(baseline_spectrogram,1),chan_counter);
    
    % calc baseline for each chan
    parfor chan =1:chan_counter
        [~,f,t,ps] = spectrogram(data(chan, :),win,ovrlp,F,fs);
        baseline_spectrogram(:,:,chan) = 20*log10(ps);
    end
    
    % save mean & std for each freq
    for chan = 1:chan_counter
        for f = 1:size(baseline_spectrogram,1)
            chan_powr_mn_spectrogram (f,chan) =  nanmean(baseline_spectrogram(f, :, chan),2);
            chan_powr_std_spectrogram(f,chan) =  nanstd(baseline_spectrogram(f, :, chan),0,2);
        end
    end
    
    cd ([path_name '/subj_' subj ])
end

% calc chan power, remove artifacts, get mean and std
clear data_new a a0 b  chan_name_1 chan_name_2 data_new ds_data desired_chans...
    flt1 flt2 flt3 hdr left_chans lock trial_data  training_behav_matrix_label train_images time_domain_chan_mn_std...
    testing_behav_matrix_label  test_images strt_time real_onsets real_offsets onsets_dur onsets onset_idx onset record...
    offset offset_idx offsets offsets_dur  end_time right_chans

%%  mean and std
chan_powr_mn  = zeros(length(freq),chan_counter);
chan_powr_std = zeros(length(freq),chan_counter);

chan_powr_mn3a  = zeros(length(freq),chan_counter);
chan_powr_mn2a  = zeros(length(freq),chan_counter);
chan_powr_mn1  = zeros(length(freq),chan_counter);
chan_powr_mn2  = zeros(length(freq),chan_counter);
chan_powr_mn3  = zeros(length(freq),chan_counter);
chan_powr_mn4  = zeros(length(freq),chan_counter);
chan_powr_std3a  = zeros(length(freq),chan_counter);
chan_powr_std2a  = zeros(length(freq),chan_counter);
chan_powr_std1  = zeros(length(freq),chan_counter);
chan_powr_std2  = zeros(length(freq),chan_counter);
chan_powr_std3  = zeros(length(freq),chan_counter);
chan_powr_std4  = zeros(length(freq),chan_counter);


cond2a_prestims_trial_artif =[];
cond3a_prestims_trial_artif = [];
cond1_prestims_trial_artif = [];
cond2_prestims_trial_artif = [];
cond3_prestims_trial_artif = [];
cond4_prestims_trial_artif = [];
        
        
if normalization ==1

    parfor chan =  1:chan_counter
        cwt = cwtft({data(chan, :),dt},...
            'scales',scales,'wavelet','morl');
        clean_points =  artifact_points(chan,:) == 0;
        temp =  single(10*log10(abs(cwt.cfs(:,clean_points)).^2));
        chan_powr_mn  (:, chan)  = nanmean(temp, 2);
        chan_powr_std (:, chan)  = nanstd (temp, 1, 2);
    end
    
elseif normalization ==2 ||  normalization ==3
    % find response indices
    responses          = [training_behav_matrix(:,4); testing_behav_matrix(:,4)];
    resp_idx_tmp       = fs*([training_behav_matrix(:,4); testing_behav_matrix(:,4)]+ISI);
    test_trial_idx     = 1:nTrials_test+nTrials_study;
    resp_idx           = off_idx(test_trial_idx)+resp_idx_tmp';
    pre_stim_dur       = [resp_idx(1:end-1)' on_idx(2:end)'];
    
    if normalization ==2
          pre_stim_dur_idx = [];
        for iRow = 1:size(pre_stim_dur,1)
            pre_stim_dur_idx =[pre_stim_dur_idx pre_stim_dur(iRow,1):pre_stim_dur(iRow,2)];
        end
        
        parfor chan =  1:chan_counter
            cwt = cwtft({data(chan, :),dt},...
                'scales',scales,'wavelet','morl');
            temp =  single(10*log10(abs(cwt.cfs).^2));
            bad_points =  artifact_points(chan,:) == 1;
            temp(:,logical(bad_points)) = nan;
            chan_powr_mn  (:, chan)  = nanmean(temp(:,round(pre_stim_dur_idx)'), 2);
            chan_powr_std (:, chan)  = nanstd (temp(:,round(pre_stim_dur_idx)'), 1, 2);
        end
    elseif normalization ==3

        % all prestim durs
        pre_stim_dur = [[0 0]; pre_stim_dur];
        
        % log vec for all conds
        [cond3a_log_vec_temp,cond2a_log_vec,cond1_log_vec,cond2_log_vec,cond3_log_vec,cond4_log_vec] =...
            get_trial_condition_label(subj);
        
        % corr log vec size to include training and test
        cond3a_log_vec = [cond3a_log_vec_temp' zeros(1,length(cond1_log_vec))];
        cond2a_log_vec = [cond2a_log_vec' zeros(1,length(cond1_log_vec))];
        cond1_log_vec  = [zeros(1,length(cond3a_log_vec_temp)) cond1_log_vec'];
        cond2_log_vec  = [zeros(1,length(cond3a_log_vec_temp)) cond2_log_vec'];
        cond3_log_vec  = [zeros(1,length(cond3a_log_vec_temp)) cond3_log_vec'];
        cond4_log_vec  = [zeros(1,length(cond3a_log_vec_temp)) cond4_log_vec'];
         
        cond2a_prestims = pre_stim_dur(logical(cond2a_log_vec),:);
        cond3a_prestims = pre_stim_dur(logical(cond3a_log_vec),:);
        cond1_prestims = pre_stim_dur(logical(cond1_log_vec),:);
        cond2_prestims = pre_stim_dur(logical(cond2_log_vec),:);
        cond3_prestims = pre_stim_dur(logical(cond3_log_vec),:);
        cond4_prestims = pre_stim_dur(logical(cond4_log_vec),:);
        
        cond2a_prestims_trial_artif = ones(chan_counter, size(cond2a_prestims,1));
        cond3a_prestims_trial_artif = ones(chan_counter,size(cond3a_prestims,1));
        cond1_prestims_trial_artif = ones(chan_counter,size(cond1_prestims,1));
        cond2_prestims_trial_artif = ones(chan_counter,size(cond2_prestims,1));
        cond3_prestims_trial_artif = ones(chan_counter,size(cond3_prestims,1));
        cond4_prestims_trial_artif = ones(chan_counter,size(cond4_prestims,1));

        % find which prestims have artifact
        for chan =  1:chan_counter
            
            % condition 1
            trial_dur = 2.4;
            for iTrl = 1:size(cond2a_prestims,1)
                if cond2a_prestims(iTrl,1)==0
                    continue
                else
                    if sum(artifact_points(chan,round(cond2a_prestims(iTrl,1)):round(cond2a_prestims(iTrl,2))+trial_dur*fs))==0;
                        cond2a_prestims_trial_artif(chan,iTrl) = 0;
                    end
                end
            end
            
            % condition 2
            for iTrl = 1:size(cond3a_prestims,1)
                if cond3a_prestims(iTrl,1)==0
                    cond3a_prestims_trial_artif(chan,iTrl) = 1;
                    continue
                else
                    if sum(artifact_points(chan,round(cond3a_prestims(iTrl,1)):round(cond3a_prestims(iTrl,2))+trial_dur*fs))==0;
                        cond3a_prestims_trial_artif(chan,iTrl) = 0;
                    end
                end
            end
            
            % condition 3
            for iTrl = 1:size(cond1_prestims,1)
                if cond1_prestims(1,1)==0
                    continue
                else
                    if sum(artifact_points(chan,round(cond1_prestims(iTrl,1)):round(cond1_prestims(iTrl,2))+trial_dur*fs))==0;
                        cond1_prestims_trial_artif(chan,iTrl) = 0;
                    end
                end
            end
            
            % condition 4
            for iTrl = 1:size(cond2_prestims,1)
                if cond2_prestims(1,1)==0
                    continue
                else
                    if sum(artifact_points(chan,round(cond2_prestims(iTrl,1)):round(cond2_prestims(iTrl,2))+trial_dur*fs))==0;
                        cond2_prestims_trial_artif(chan,iTrl) = 0;
                    end
                end
            end
            
            % condition 5
            for iTrl = 1:size(cond3_prestims,1)
                if cond3_prestims(1,1)==0
                    continue
                else
                    if sum(artifact_points(chan,round(cond3_prestims(iTrl,1)):round(cond3_prestims(iTrl,2))+trial_dur*fs))==0;
                        cond3_prestims_trial_artif(chan,iTrl) = 0;
                    end
                end
            end
            
            % condition 6
            for iTrl = 1:size(cond4_prestims,1)
                if cond4_prestims(1,1)==0
                    continue
                else
                    if sum(artifact_points(chan,round(cond4_prestims(iTrl,1)):round(cond4_prestims(iTrl,2))+trial_dur*fs))==0;
                        cond4_prestims_trial_artif(chan,iTrl) = 0;
                    end
                end
            end
            
        end
        
        
        
        parfor chan =  1:chan_counter
            cwt = cwtft({data(chan, :),dt},...
                'scales',scales,'wavelet','morl');
            temp =  single(10*log10(abs(cwt.cfs).^2));
            
            % cond2a
            cond2a_prestims_temp = cond2a_prestims(cond2a_prestims_trial_artif(chan,:)==0,:);
            cond2a_prestim_all_indices = [];
            for iPrestims = 1:size(cond2a_prestims_temp,1)
                cond2a_prestim_all_indices = [cond2a_prestim_all_indices cond2a_prestims_temp(iPrestims,1):cond2a_prestims_temp(iPrestims,2)];
            end
            chan_powr_mn2a (:, chan)= nanmean(temp(:,round(cond2a_prestim_all_indices)), 2);
            chan_powr_std2a(:, chan)= nanstd(temp(:,round(cond2a_prestim_all_indices)),1,2);
            
            
            % cond3a
            cond3a_prestims_temp = cond3a_prestims(cond3a_prestims_trial_artif(chan,:)==0,:);
            cond3a_prestim_all_indices = [];
            for iPrestims = 1:size(cond3a_prestims_temp,1)
                cond3a_prestim_all_indices = [cond3a_prestim_all_indices cond3a_prestims_temp(iPrestims,1):cond3a_prestims_temp(iPrestims,2)];
            end
            
            chan_powr_mn3a (:, chan)= nanmean(temp(:,round(cond3a_prestim_all_indices)), 2);
            chan_powr_std3a(:, chan)= nanstd(temp(:,round(cond3a_prestim_all_indices)),1,2);
            
            % cond1
            cond1_prestims_temp = cond1_prestims(cond1_prestims_trial_artif(chan,:)==0,:);
            cond1_prestim_all_indices = [];
            for iPrestims = 1:size(cond1_prestims_temp,1)
                cond1_prestim_all_indices = [cond1_prestim_all_indices cond1_prestims_temp(iPrestims,1):cond1_prestims_temp(iPrestims,2)];
            end
            
            chan_powr_mn1 (:, chan)= nanmean(temp(:,round(cond1_prestim_all_indices)), 2);
            chan_powr_std1(:, chan)= nanstd(temp(:,round(cond1_prestim_all_indices)),1,2);
            
            % cond2
            cond2_prestims_temp = cond2_prestims(cond2_prestims_trial_artif(chan,:)==0,:);
            cond2_prestim_all_indices = [];
            for iPrestims = 1:size(cond2_prestims_temp,1)
                cond2_prestim_all_indices = [cond2_prestim_all_indices cond2_prestims_temp(iPrestims,1):cond2_prestims_temp(iPrestims,2)];
            end
            
            chan_powr_mn2 (:, chan)= nanmean(temp(:,round(cond2_prestim_all_indices)), 2);
            chan_powr_std2(:, chan)= nanstd(temp(:,round(cond2_prestim_all_indices)),1,2);
            
            % cond3
            cond3_prestims_temp = cond3_prestims(cond3_prestims_trial_artif(chan,:)==0,:);
            cond3_prestim_all_indices = [];
            for iPrestims = 1:size(cond3_prestims_temp,1)
                cond3_prestim_all_indices = [cond3_prestim_all_indices cond3_prestims_temp(iPrestims,1):cond3_prestims_temp(iPrestims,2)];
            end
            chan_powr_mn3 (:, chan)= nanmean(temp(:,round(cond3_prestim_all_indices)), 2);
            chan_powr_std3(:, chan)= nanstd(temp(:,round(cond3_prestim_all_indices)),1,2);
            
            chan_powr_mn3 (:, chan)= nanmean(temp(:,round(cond3_prestim_all_indices)), 2);
            chan_powr_std3(:, chan)= nanstd(temp(:,round(cond3_prestim_all_indices)),1,2);
            
            % cond4
            cond4_prestims_temp = cond4_prestims(cond4_prestims_trial_artif(chan,:)==0,:);
            cond4_prestim_all_indices = [];
            for iPrestims = 1:size(cond4_prestims_temp,1)
                cond4_prestim_all_indices = [cond4_prestim_all_indices cond4_prestims_temp(iPrestims,1):cond4_prestims_temp(iPrestims,2)];
            end
            
            chan_powr_mn4 (:, chan)= nanmean(temp(:,round(cond4_prestim_all_indices)), 2);
            chan_powr_std4(:, chan)= nanstd(temp(:,round(cond4_prestim_all_indices)),1,2);
        end
        
        
%         %% concatinate pre_stims assoc with each cond
%        
%         cond2a_prestim_all_indices=[];
%         for iPrestims = 1:size(cond2a_prestims,1)
%             if  cond2a_prestims(iPrestims,1)==0 % remove first prestim 0
%                 disp('cond2a_prestims')
%                 continue
%             else
%                 cond2a_prestim_all_indices = [cond2a_prestim_all_indices cond2a_prestims(iPrestims,1):cond2a_prestims(iPrestims,2)];
%             end
%         end
%         
%         cond3a_prestim_all_indices=[];
%         for iPrestims = 1:size(cond3a_prestims,1)
%             if  cond3a_prestims(iPrestims,1)==0  % remove first prestim 0
%                  disp('cond3a_prestims')
%                 continue
%             else
%                 cond3a_prestim_all_indices = [cond3a_prestim_all_indices cond3a_prestims(iPrestims,1):cond3a_prestims(iPrestims,2)];
%             end
%         end
%         
%         cond1_prestim_all_indices=[];
%         for iPrestims = 1:size(cond1_prestims,1)
%             cond1_prestim_all_indices = [cond1_prestim_all_indices cond1_prestims(iPrestims,1):cond1_prestims(iPrestims,2)];
%         end
% 
%         cond2_prestim_all_indices=[];
%         for iPrestims = 1:size(cond2_prestims,1)
%             cond2_prestim_all_indices = [cond2_prestim_all_indices cond2_prestims(iPrestims,1):cond2_prestims(iPrestims,2)];
%         end
% 
%         cond3_prestim_all_indices=[];
%         for iPrestims = 1:size(cond3_prestims,1)
%             cond3_prestim_all_indices = [cond3_prestim_all_indices cond3_prestims(iPrestims,1):cond3_prestims(iPrestims,2)];
%         end
% 
%         cond4_prestim_all_indices=[];
%         for iPrestims = 1:size(cond4_prestims,1)
%             cond4_prestim_all_indices = [cond4_prestim_all_indices cond4_prestims(iPrestims,1):cond4_prestims(iPrestims,2)];
%         end
%         
        
        
                
       
    end
end
save(['normalization_' norm{normalization} '_ref_' ref '_baseline_info_wavelet_' num2str(NumVoices) 'num_20logdb_' num2str(minfreq) 'hz_' num2str(maxfreq) 'hz'  '_notched_artifact_reject_subj_'...
    subj '_select_chan_' num2str(select_chan)  '_fs_' num2str(fs)], 'chan_powr_mn', 'chan_powr_std','chan_artifact_thresh', 'freq', 'on_idx', 'off_idx', 'chan_label',...
    'chan_powr_mn3a','chan_powr_mn2a','chan_powr_mn1','chan_powr_mn2','chan_powr_mn3','chan_powr_mn4',...
    'chan_powr_std3a','chan_powr_std2a','chan_powr_std1','chan_powr_std2','chan_powr_std3','chan_powr_std4', ...
    'cond1_prestims_trial_artif', 'cond2_prestims_trial_artif', 'cond3_prestims_trial_artif', 'cond4_prestims_trial_artif',...
    'cond2a_prestims_trial_artif', 'cond3a_prestims_trial_artif');
disp('done')
toc