clear all; close all; clc
tic
subj = '87' % 84, 85, 87
lock = 'onset'
downsample = 'yes' 
ref  = 'LM'%'WM' 'LM'[];
clinical = 'yes'
select_chan   = 3 % 1: get specified chans indicated by chan_name 3: gets all chans
normalization = 3% 1 = entire task 2 = prestim
norm = {'entire_recording' 'prestim' 'cond_spec_prestim'}
addpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final'))
ft_defaults
script_path = '/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-20181130/fileio/private'
file_path   = (['/mnt/yassamri/iEEG/sandra/subj_' subj '/research_dataset'])
if strcmp('',downsample)
    cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
    load('timestampsDS8_FS1000.mat') % load downsampled timestamps which we will upsample
end
cd(script_path)

% get photodiode
if strcmp('84',subj)
    photodiode_fn= '/photo1.ncs';
    LFP_fn       = '/LHH1.ncs';
elseif strcmp('85',subj)
    photodiode_fn=   '/photo1_0001.ncs';
    LFP_fn       = '/LHH1_0001.ncs';
elseif strcmp('87',subj)
    photodiode_fn=   '/photo1_0004.ncs';
    LFP_fn       = '/LHH1_0004.ncs';
end
photodiode   = read_neuralynx_ncs([file_path photodiode_fn]);
fs           = unique(photodiode.SampFreq);

if strcmp ('84', subj)
    photodiode_signal = -1*(photodiode.dat(:))';
else
    photodiode_signal = (photodiode.dat(:))';  
end


clear photodiode

if strcmp('yes',downsample)
   D = 8;
   fctr=1;
else
   D = 1; 
   fctr=8;
end

photodiode_signal = decimate(photodiode_signal, D);
fs = fs/D
figure; hold on
plot(photodiode_signal)

% manually indicate start point and end point
if strcmp ('84', subj)
    end_DataIndex     = length(photodiode_signal);
    start.DataIndex = 3154321*fctr;
    thresh_mltpl =  7;
    photodiode_signal = photodiode_signal(start.DataIndex : end_DataIndex);
    
elseif strcmp ('85', subj)
    start.DataIndex   = 388726*fctr;
    end_DataIndex     = 1589996*fctr;
    photodiode_signal = photodiode_signal(start.DataIndex : end_DataIndex);
    thresh_mltpl = 7.9;
    
elseif strcmp ('87', subj)
    start.DataIndex   = 11111346*fctr;
    end_DataIndex     = 12308003*fctr;
    photodiode_signal = photodiode_signal(start.DataIndex : end_DataIndex);
    thresh_mltpl = 7.9;
end


% dont re-calc timestamps for non-downsampled data
if strcmp('yes',downsample)
    a      = diff(photodiode_signal);
    meana  = nanmean(a);
    stda   = nanstd(a);
    
    %clear real_onsets real_offsets on_idx
    onset  = (a>(meana+(thresh_mltpl*stda)));
    offset = (a<(meana-(thresh_mltpl*stda)));
    
    % measure offset
    off_idx = (diff(offset)==1);
    offsets = find(off_idx);
    offsets_dur  = diff(offsets);
    
    real_offsets = offsets(offsets_dur>.5*mean(offsets_dur));
    offset_idx   = nan(1,length(photodiode_signal));
    offset_idx(real_offsets) = 1;
    
    % measure onset
    on_idx = (diff(onset)==1);
    onsets = find(on_idx);
    onsets_dur  = diff(onsets);
    real_onsets = onsets(onsets_dur>.5*mean(onsets_dur));
    onset_idx   = nan(1,length(photodiode_signal));
    onset_idx(real_onsets) = 1;
    
    sum(on_idx)
    sum(off_idx)
    
    
    % if there are a few extra pnsets/offsets --> adjust
    if strcmp('84',subj)
        % offest
        a = 3
        b = 0
        
        % onset
        c = 2
        d = 0
    elseif strcmp('85',subj)
        % offest
        a = 2
        b = 0
        
        % onset
        c = 1
        d = 0
      elseif strcmp('87',subj)
        % offest
        a = 1
        b = 0
        
        % onset
        c = 1
        d = 0      
    end
    real_offsets = real_offsets(a:end-b);
    real_onsets  = real_onsets(c:end-d);
    
    % correct first real onset for IR95
    if strcmp('85', subj)
        real_onsets(1) =  1874;
    end
    
    onset_idx  = nan(1,length(photodiode_signal));
    offset_idx = nan(1,length(photodiode_signal));
    offset_idx(real_offsets) = 1;
    onset_idx(real_onsets) = 1;
    if onsets_dur(end)> .5*mean(onsets_dur)
        onset_idx(onsets(end)) = 1;
        offset_idx(offsets(end)) = 1;
    end
    
    % plot truncated photodiote signal uncluding your task
    
    if strcmp('84',subj)
        shift = -3500;
    elseif strcmp('85',subj)
        shift = -4800;
    elseif strcmp('87',subj)
        shift = 2500;
    end
    
    figure; hold on
    plot(photodiode_signal)
    plot(onset_idx+shift, 'r*')
    plot(offset_idx+shift, 'g*')
    
    % variable to  get trials
    on_idx  = find(onset_idx==1);
    off_idx = find(offset_idx==1);

end
on_idx  = on_idx*fctr;
off_idx = off_idx*fctr;
on_idx_check  = (on_idx+((start.DataIndex-1)));
off_idx_check = (off_idx+((start.DataIndex-1)));

% check identified onsets out of total
photodiode   = read_neuralynx_ncs([file_path photodiode_fn]);
photodiode_signal = -1*(photodiode.dat(:))';
photodiode_signal = decimate(photodiode_signal, D);
onset_idx  = nan(1,length(photodiode_signal));
offset_idx = nan(1,length(photodiode_signal));

onset_idx(on_idx_check)   = 1;
offset_idx(off_idx_check) = 1;

figure; hold on
plot(photodiode_signal)
if strcmp('84', subj)
    plot(onset_idx-3500, 'r*')
    plot(offset_idx-3500, 'g*')
elseif strcmp('84', subj)
    plot(onset_idx+4700, 'r*')
    plot(offset_idx+4700, 'g*')
elseif strcmp('87', subj)
    plot(onset_idx-1400, 'r*')
    plot(offset_idx-1400, 'g*')    
    
end

% timestamps are on_idx off_idx
% behav info
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
load(['behavior_subj' subj])
nTrials_study = size(training_behav_matrix,1);
nTrials_test  = size(testing_behav_matrix,1);
nBlocks       = 2;

% get LFP data: 
cd(script_path)
LFP_path      =  [file_path '/LFP'];

if strcmp('84', subj)
    desired_chans = {'ROF','RHH','RAM', 'RAC', 'LOF','LHH', 'LAM', 'LAC'};
elseif strcmp('85', subj)
    desired_chans = {'LAM','LHH','RAC', 'RAM', 'RHH','ROF', 'RPC'};
    ext = '_0001';
elseif strcmp('87', subj)
    desired_chans = {'LAC','LAM','LHH', 'LOF', 'RAC','RAM', 'RHH', 'ROF'};
    ext = '_0004';
end

% get LFP data either from 1. research recording 2. saved (parsed) clinical recordings

if strcmp('', clinical) % if using research recording
    save_nm = '_research'
    
    chan_counter  = length(desired_chans)*8;
    lfp       = read_neuralynx_ncs([LFP_path LFP_fn]);
    temp_data = lfp.dat(:);
    temp_LFP  = decimate(temp_data,D);
    temp_LFP  = temp_LFP(start.DataIndex:end_DataIndex);
    LFP_data_temp  = zeros(chan_counter,length(temp_LFP));
    fs        = unique(lfp.SampFreq)/D;
    
    i = 0;
    chan_label = [];
    for chan = 1:length(desired_chans)
        for a = 1:8
            % data
            i = i+1;
            chan_label{i} = [desired_chans{chan} num2str(a)];
            lfp= read_neuralynx_ncs([LFP_path '/' desired_chans{chan} num2str(a) ext  '.ncs']);
            temp_data = lfp.dat(:);
            temp_LFP  = decimate(temp_data,D);
            LFP_data_temp(((chan-1)*8) + a,:) = temp_LFP(start.DataIndex:end_DataIndex);
            
            % get right and left index
            right_chans(i) = chan_label{i}(1)== 'R';
            left_chans(i)  = chan_label{i}(1)== 'L';
        end
    end
    cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
elseif strcmp('yes', clinical)
    load(['/mnt/yassamri/iEEG/sandra/subj_' subj '/clinical_dataset/parsed_data.mat'])
    LFP_data_temp    = clinical_full_data ;
    chan_label  = chan_label_clinical;
    chan_counter = length(chan_label);
    save_nm = '_clinical'
end

if strcmp('yes', downsample)
    % down sample data and timestamps to fs = 500
    D2 = 2;
    if strcmp('84', subj) % clin data is already at 500
        LFP_data = clinical_full_data;
    else
        LFP_data = zeros(size(LFP_data_temp,1), length(decimate(LFP_data_temp(1,:), D2)));
        for iElec = 1:size(LFP_data_temp,1)
            LFP_data(iElec,:) = decimate(LFP_data_temp(iElec,:), D2);
        end
    end
    
    fs = fs/D2;
    %timestamp data ***** THESE ARE FINAL DS TS
    on_idx  = round(on_idx*1/D2);
    off_idx = round(off_idx*1/D2);
end
rmpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-20181130'))

% re ref
if strcmp('WM',ref)
    [WM_ref_vector] = get_WM_reref(subj, chan_label);
    data_reref = LFP_data - LFP_data(WM_ref_vector,:);
elseif strcmp('LM', ref)
    [LM_ref_vector] = get_LM_reref(subj);
    data_reref = nan(size(LFP_data));
    for iProbe = 1:size(LM_ref_vector,1) % loop thru probe
        for iElec = LM_ref_vector(iProbe,1):LM_ref_vector(iProbe,2)
            if iElec == LM_ref_vector(iProbe,1) % if on first elec
                if (strcmp('87', subj) || strcmp('85', subj) || strcmp('84', subj)) && iElec==25
                data_reref(iElec,:) = LFP_data(iElec,:) - LFP_data(22,:); % LAC8-LAC7
                else
                data_reref(iElec,:) = LFP_data(iElec,:) - LFP_data(iElec+1,:);
                end
                
            elseif ((iElec == LM_ref_vector(iProbe,2)) && (iElec~=LM_ref_vector(iProbe,1)))% if on last elec
                data_reref(iElec,:) = LFP_data(iElec,:) - LFP_data(iElec-1,:);
            else
                data_reref(iElec,:) = LFP_data(iElec,:) - (LFP_data(iElec-1,:)+LFP_data(iElec+1,:))/2;
            end
            disp(iElec)
        end
    end
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

% notch filter and plot 
figure;
hold on
yscale = 800;
data = nan(size(data_reref));

for chan = 1:chan_counter
    data (chan,:) = filtfilt(flt3, filtfilt(flt2, filtfilt(flt1,data_reref(chan,:))));
    plot(0:1/fs:(size(data_reref,2)-1)/fs, data_reref(chan,:)+(yscale*chan))
end

artifact_points      = zeros(size(data,1), size(data,2));
chan_artifact_thresh = zeros(chan_counter,2);
sec   = 1; % duration to remove around artifacts

% artif thresh
if strcmp('44', subj); mltpl = 3;
elseif strcmp('84', subj) || strcmp('84', subj); mltpl = 3.5;
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
    
    close all
    %  preview last chan
    figure; hold on; plot(0:1/fs:(size(data,2)-1)/fs, x)
    line([0 (size(data,2)-1)/fs] ,[pos_thresh pos_thresh], 'Color','k','LineWidth',4)
    line([0 (size(data,2)-1)/fs] ,[neg_thresh neg_thresh], 'Color','k','LineWidth',4)
    
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
    % close all
end

% wavelet params
dt = 1/fs;
NumVoices = 32;
a0 = 2^(1/NumVoices);
wavCenterFreq = 6/(2*pi);
minfreq = 3;
maxfreq = 200;
minscale = wavCenterFreq/(maxfreq*dt); 
maxscale = wavCenterFreq/(minfreq*dt); 
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale)); 
scales = a0.^(minscale:maxscale).*dt;
freq = wavCenterFreq./(fs*scales.*dt);

ISI                = .5;
% get trial data
if strcmp ('onset',lock)
    pre_stim  = 0.9; 
    post_stim = 2.4; 
    strt_time = on_idx' - pre_stim*fs;
    end_time  = on_idx' + post_stim*fs;
elseif strcmp ('onset2',lock)
    pre_stim  = 0.9; 
    post_stim = 3; 
    strt_time = on_idx' - pre_stim*fs;
    end_time  = on_idx' + post_stim*fs;
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

trial_data     = single(nan(length(on_idx), length(strt_time(1):end_time(1)), size(data,1))); % trail x time x chan
for chan = 1:chan_counter
    for trial = 1:length(on_idx)
        trial_data(trial,:,chan) = data(chan, strt_time(trial):end_time(trial));
    end
end

% save trial data
cd (['/mnt/yassamri/iEEG/sandra/subj_' subj ])
if strcmp('onset',lock) || strcmp('onset2',lock)
    if strcmp('',downsample)
        save(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_' num2str(select_chan) '_NotDownsampled' save_nm],...
            'trial_data', 'fs', 'chan_counter', 'nTrials_study', 'nTrials_test', 'nBlocks', 'chan_label',...
            'pre_stim', 'post_stim', 'lock', 'on_idx' , 'off_idx', 'artifact_points', 'chan_artifact_thresh',...
            'chan_artifact_thresh', 'freq',  'minfreq', 'maxfreq', '-v7.3');
    elseif strcmp('yes',downsample)
        save(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_' num2str(select_chan)  save_nm '_fs_' num2str(fs)],...
            'trial_data', 'fs', 'chan_counter', 'nTrials_study', 'nTrials_test', 'nBlocks', 'chan_label',...
            'pre_stim', 'post_stim', 'lock', 'on_idx' , 'off_idx', 'artifact_points', 'chan_artifact_thresh',...
            'chan_artifact_thresh', 'freq',  'minfreq', 'maxfreq');
    end
elseif strcmp('response',lock)
    if strcmp('',downsample)
        save(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_' num2str(select_chan) '_NotDownsampled' save_nm],...
            'trial_data', 'fs', 'chan_counter', 'nTrials_study', 'nTrials_test', 'nBlocks', 'chan_label',...
            'pre_stim', 'post_stim',  'lock', 'on_idx' , 'off_idx', 'artifact_points', 'chan_artifact_thresh',...
            'chan_artifact_thresh', 'freq',  'minfreq', 'maxfreq','responses', 'resp_idx','-v7.3');
    elseif strcmp('yes',downsample)
        save(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_' num2str(select_chan) save_nm '_fs_' num2str(fs)],...
            'trial_data', 'fs', 'chan_counter', 'nTrials_study', 'nTrials_test', 'nBlocks', 'chan_label',...
            'pre_stim', 'post_stim', 'lock', 'on_idx' , 'off_idx', 'artifact_points', 'chan_artifact_thresh',...
            'chan_artifact_thresh', 'freq',  'minfreq', 'maxfreq','responses', 'resp_idx');
    end
end


%% get baseline values
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

if normalization == 1
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
    end
end
cd (['/mnt/yassamri/iEEG/sandra/subj_' subj ])
save(['normalization_' norm{normalization} '_ref_' ref '_baseline_info_wavelet_' num2str(NumVoices) 'num_20logdb_' num2str(minfreq) 'hz_' num2str(maxfreq) 'hz'  '_notched_artifact_reject_subj_'...
    subj '_select_chan_' num2str(select_chan) save_nm '_fs_' num2str(fs)], 'chan_powr_mn', 'chan_powr_std','chan_artifact_thresh', 'freq', 'on_idx', 'off_idx', 'chan_label',...
     'chan_powr_mn3a','chan_powr_mn2a','chan_powr_mn1','chan_powr_mn2','chan_powr_mn3','chan_powr_mn4',...
    'chan_powr_std3a','chan_powr_std2a','chan_powr_std1','chan_powr_std2','chan_powr_std3','chan_powr_std4',...
        'cond1_prestims_trial_artif', 'cond2_prestims_trial_artif', 'cond3_prestims_trial_artif', 'cond4_prestims_trial_artif',...
    'cond2a_prestims_trial_artif', 'cond3a_prestims_trial_artif');
disp('done')
toc
toc
disp('done')