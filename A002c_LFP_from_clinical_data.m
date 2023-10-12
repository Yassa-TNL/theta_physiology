% patient specific
clear all; close all;clc
subj = '85'
chans_to_examin = 6;
downsample  = 'yes'

script_path = '/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-20181130/fileio/private'
file_path   = (['/mnt/yassamri/iEEG/sandra/subj_' subj '/research_dataset'])

% use fieldtrip
addpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final'))
ft_defaults
cd(script_path)

% get research data
if strcmp('84',subj)
    photodiode_fn= '/photo1.ncs';
    LFP_fn       = '/LHH1.ncs';
elseif strcmp('85',subj)
    photodiode_fn=   '/photo1_0001.ncs';
    LFP_fn       = '/LHH1_0001.ncs';
elseif strcmp('87',subj)
    photodiode_fn= '/photo1_0004.ncs';
    LFP_fn       = '/LHH1_0004.ncs';
end

photodiode   = read_neuralynx_ncs([file_path photodiode_fn]);
fs = unique(photodiode.SampFreq);
if strcmp ('84', subj)
    photodiode_signal = -1*(photodiode.dat(:))';
else
    photodiode_signal = (photodiode.dat(:))';
end
clear photodiode

% get LFP data:
cd(script_path); LFP_path = [file_path '/LFP'];

if strcmp('84', subj)
    desired_chans = {'ROF','RHH','RAM', 'RAC', 'LOF','LHH', 'LAM', 'LAC'};
    ext = '';
elseif strcmp('85', subj)
    desired_chans = {'LAM','LHH','RAC', 'RAM', 'RHH','ROF', 'RPC'};
    ext = '_0001';
elseif strcmp('87', subj)
    desired_chans = {'LAC','LAM','LHH', 'LOF', 'RAC','RAM', 'RHH', 'ROF'};
    ext = '_0004';
end

% check whether discon is during recording; go to end of script to plot
i = 0;
for chan = 1:length(desired_chans)
    for a = 1:8
        i = i+1;
        lfp = read_neuralynx_ncs([LFP_path '/' desired_chans{chan} num2str(a)  ext '.ncs']);
        idx = find(lfp.NumValidSamp<512);
        nValSamp = lfp.NumValidSamp(idx);
        nTrials  = length(idx);
        trl           = zeros(nTrials, 2);
        trl(1, 1)     = 1;
        trl(:, 2)     = ((idx-1)*512 + nValSamp)';
        trl(2:end, 1) = (idx(1:end-1)*512+1)';
        trls{i} = trl; % timestamps of discontinuties 
        clear trl
    end
end
m = lfp.TimeStamp;
mode(diff(m));

if strcmp('yes',downsample)
    D = 8; fctr = 1;
else
    D = 1; fctr = 8;
end

photodiode_signal = decimate(photodiode_signal, D);
fs = fs/D;
figure; hold on
plot(photodiode_signal)

% manually indicate start point and end point
if strcmp ('84', subj)
    end_DataIndex     = length(photodiode_signal);
    start.DataIndex   = 3154321*fctr;
    thresh_mltpl      =  7;
    photodiode_signal = photodiode_signal(start.DataIndex : end_DataIndex);
elseif strcmp ('85', subj)
    start.DataIndex   = 388726*fctr;
    end_DataIndex     = 1589996*fctr;
    photodiode_signal = photodiode_signal(start.DataIndex : end_DataIndex);
    thresh_mltpl      = 7.9;
elseif strcmp ('87', subj)
    start.DataIndex   = 11111346*fctr;
    end_DataIndex     = 12308003*fctr;
    photodiode_signal = photodiode_signal(start.DataIndex : end_DataIndex);
    thresh_mltpl      = 7.9;
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
    sum(on_idx);
    sum(off_idx);
    
    % if there are a few extra onsets/offsets --> adjust
    if strcmp('84',subj)
        % offest   % onset
        a = 3; b = 0; c = 2; d = 0;
    elseif strcmp('85',subj)
        % offest   % onset
        a = 2; b = 0; c = 1; d = 0;
    elseif strcmp('87',subj)
        % offest    % onset
        a = 1; b = 0;  c = 1; d = 0;
    end
    real_offsets = real_offsets(a:end-b);
    real_onsets  = real_onsets(c:end-d);
    
    % correct first real onset for IR85
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
    
    % plot truncated photodiote signal including your task
        if strcmp('84',subj)
        shift = -3500;
    elseif strcmp('85',subj)
        shift = -4800;
    elseif strcmp('87',subj)
        shift = 2500;
    end
    
    figure; hold on;plot(photodiode_signal);plot(onset_idx+shift, 'r*'); plot(offset_idx+shift, 'g*')
        
    % variable to  get trials
    on_idx  = find(onset_idx==1);
    off_idx = find(offset_idx==1);
    
end
on_idx  = on_idx*fctr;
off_idx = off_idx*fctr;


% check identified onsets out of total
on_idx_check  = (on_idx+((start.DataIndex-1)));
off_idx_check = (off_idx+((start.DataIndex-1)));
photodiode    = read_neuralynx_ncs([file_path photodiode_fn]);
photodiode_signal_all = -1*(photodiode.dat(:))';
photodiode_signal_all = decimate(photodiode_signal_all, D);

onset_idx  = nan(1,length(photodiode_signal_all));
offset_idx = nan(1,length(photodiode_signal_all));
onset_idx(on_idx_check)   = 1;
offset_idx(off_idx_check) = 1;

figure; hold on
plot(photodiode_signal_all)
if strcmp('84', subj)
    plot(onset_idx-3500, 'r*')
    plot(offset_idx-3500, 'g*')
elseif strcmp('85', subj)
    plot(onset_idx+4700, 'r*')
    plot(offset_idx+4700, 'g*')
elseif strcmp('87', subj)
    plot(onset_idx-1400, 'r*')
    plot(offset_idx-1400, 'g*')
    
end
% timestamps are on_idx off_idx

% get LFP data:
LFP_path      =  [file_path '/LFP'];

if strcmp('84', subj)
    desired_chans = {'ROF','RHH','RAM', 'RAC', 'LOF','LHH', 'LAM', 'LAC'};
    ext = '';
elseif strcmp('85', subj)
    desired_chans = {'LAM','LHH','RAC', 'RAM', 'RHH','ROF', 'RPC'};
    ext = '_0001';
elseif strcmp('87', subj)
    desired_chans = {'LAC','LAM','LHH', 'LOF', 'RAC','RAM', 'RHH', 'ROF'};
    ext = '_0004';
end

chan_counter  = length(desired_chans)*8;
lfp       = read_neuralynx_ncs([LFP_path LFP_fn]);
temp_data = lfp.dat(:);
temp_LFP  = decimate(temp_data,D);

% get an approp range of neuralynx data to croscorr w/ clinical data
strt = start.DataIndex;
ed   = end_DataIndex;
lfp_research_temp = temp_LFP(strt:ed); % for clinical data align

% neuralynxchanel label
fs         = unique(lfp.SampFreq)/D;
i          = 0;
chan_label_research = [];
lfp_research  = zeros(chan_counter,length(lfp_research_temp));

for chan = 1:length(desired_chans)
    for a = 1:8
        % chans
        i = i+1;
        chan_label_research{i} = [desired_chans{chan} num2str(a)];
        
        % data
        lfp = read_neuralynx_ncs([LFP_path '/' desired_chans{chan} num2str(a)  ext '.ncs']);
        temp_data = lfp.dat(:);
        temp_LFP  = decimate(temp_data,D);
        lfp_research(((chan-1)*8) + a,:) = temp_LFP(strt:ed);
    end
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% clinical data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
addpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final'));
ft_defaults
if strcmp('87',subj)
    date = 'PH2 Feb12-Feb21';
    fna  = 'BA03176';%finalized
    fnb  = {'V'};
elseif strcmp('84',subj) 
    date = 'PH2 OCT23-NOV1';
    fna = 'BA0316R';%finalized
    fnb = {'M'};
elseif strcmp('85',subj)
    date = 'PH2 jan7-jan18';
    fna  = 'BA0316Z';
    fnb  = {'V'}; % finalized
end

if strcmp('',subj) % use this set up to loop thru all files for new patients
    direc  = ['/mnt/yassamri/iEEG/sandra/subj_' subj '/clinical_dataset/' date '/NKT/EEG2100/'];
    cd(direc)
    files=dir('*.EEG*');
    listing = cell(1,size(files,1));
    for iFile = 1:size(files,1)
        listing{iFile}=files(iFile).name;
    end
else % patient for which file is found
    direc  = ['/mnt/yassamri/iEEG/sandra/subj_' subj '/clinical_dataset/' date '/NKT/EEG2100/'];
    listing = 1;
end
%%
for fnb_cntr=1
    addpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final'));
    ft_defaults
    if strcmp('',subj)
        file = [direc  listing{fnb_cntr} ];
    else
        file = [direc fna fnb{1} '.EEG'];
    end
    % init lfp mtx
    lfp_clinical_all = ft_read_data (file); % all elecs for this 2-hr recording
    hdr_info  = ft_read_header (file);
    fs_clin = hdr_info.Fs;
    new_D = fs/fs_clin;
    
    % downsample orig data if clinical data has dif sam rate
    if ~(new_D==1)
        %orig data
        for iElec = 1:size(lfp_research,1)
            lfp_research_DS(iElec,:) = decimate(lfp_research(iElec,:), new_D);
        end
        %timestamp data ***** THESE ARE FINAL DS TS
        on_idx_DS  = round(on_idx*1/new_D);
        off_idx_DS = round(off_idx*1/new_D);
        
        %check downsamples on/off idx
        photodiode_signal_DS = decimate(photodiode_signal, new_D);
        figure
        plot(photodiode_signal_DS)
        if strcmp('84', subj)
            end_DataIndex     = length(photodiode_signal_DS);
            start.DataIndex   = round(3154321*(1/new_D));
        end
        
        onset_idx_DS  = nan(1,length(photodiode_signal_DS));
        offset_idx_DS = nan(1,length(photodiode_signal_DS));
        onset_idx_DS (on_idx_DS)   = 1;
        offset_idx_DS(off_idx_DS)  = 1;
        hold on
        plot(onset_idx_DS+shift, 'r*')
        plot(offset_idx_DS+shift, 'g*')
    else
        lfp_research_DS = lfp_research;
    end
    chan_label_clinical=hdr_info.label;
        
    % pair research and clin labels
    elec_pairs = [];
    for iElec =1:length(chan_label_research)
        for iElec2 = 1:length(chan_label_clinical)
            if strcmp(chan_label_clinical{iElec2},chan_label_research{iElec})
                elec_pairs = [ elec_pairs ; iElec iElec2];
            end
        end
    end
    lfp_clinical = lfp_clinical_all(elec_pairs(1:7,2),:);
    
    rmpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-20181130'));
    rmpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-lite-20190920'));
    cd('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final');     
   
    figure
    i = 3;
    [covar, lags] = xcov(lfp_clinical(i,:), lfp_research_DS(i,:));
    plot(lags, covar); title(['subj' subj ' '  'raw - chan' num2str(i)]);
    %saveas(gcf,[direc '/iEEG/' num2str(fnb_cntr) '_' num2str(listing{fnb_cntr}) '.png'] )
    %close all
end

%%               ************** run this after u find right file **************
% convolve for 1 chan - bipolar references
i = 2
lfp_clinical_bp_ref = lfp_clinical(i,:)    - lfp_clinical(i+1,:);
lfp_research_bp_ref = lfp_research_DS(i,:) - lfp_research_DS(i+1,:);
conv_temp = conv(flip(lfp_clinical_bp_ref),lfp_research_bp_ref, 'full');
figure;
plot(conv_temp);
title(['subj' subj ' ' fna fnb{fnb_cntr} ' conv ']);

if strcmp('84',subj)
    offset_clin_research = 0;
    yshift  = -16^-4;
    yscale = -10^-6;
    conv_temp = conv_temp*-1;

elseif strcmp('85',subj)
    offset_clin_research = 1;
    yshift  = 10^-5;
    yscale = -10^-6;
    conv_temp = conv_temp*-1;
elseif strcmp('87',subj)
    offset_clin_research = 2; % 5 prev
    yshift  = 10^-5;
    yscale = -10^-6;
    conv_temp = conv_temp*-1;
end

% index data and preview 1 second in each data type
for i = 1:8
lfp_clinical_bp_ref = lfp_clinical(i,:)    - lfp_clinical(i+1,:);
lfp_research_bp_ref = lfp_research_DS(i,:) - lfp_research_DS(i+1,:);
flpClinical = flip(lfp_clinical_bp_ref);
match_lag   = find(conv_temp==max(conv_temp))
clinical_idx_range = match_lag-length(lfp_research_bp_ref)+1-offset_clin_research:match_lag-offset_clin_research;
clinical_region_conv_match = flip(flpClinical(clinical_idx_range));
% overlay re-ref clinical and research
range = 1:2000;
figure;plot(clinical_region_conv_match(range)+yshift)
hold on;plot(lfp_research_bp_ref(range)*yscale)
legend('clinical', 'research')
end


% check if non bipolared should also align
flpClinical_raw = flip(lfp_clinical,2);
clinical_indexed_raw = flip(flpClinical_raw(:,clinical_idx_range),2);
for i = 1:4
figure;plot(clinical_indexed_raw(i,range))
hold on;plot(lfp_research_DS(i,range)*yscale)
legend('clinical', 'research')
end
%% psd of research
[p1,f1] = pspectrum(lfp_research_bp_ref, fs) ;
figure;plot(f1(1:20),p1(1:20))

% psd of clinical
figure
[p2,f2] = pspectrum(clinical_region_conv_match, fs) ;
hold on
plot(f2(1:20),p2(1:20))

%% parse and save clinical data
FLP_lfp_clinical_full_data = flip(lfp_clinical_all,2);
clinical_full_data         = flip(FLP_lfp_clinical_full_data(:,clinical_idx_range),2);

%% overlay entire signal
figure;hold on
if strcmp ('84', subj)
    i1 = 9; i2 = 33;
elseif strcmp ('85', subj)
    i1 = 9; i2 = 9;
elseif strcmp ('87', subj)
    i1 = 15; i2 = 1;
end
subplot(2,1,1);plot(clinical_full_data(i1,:)-clinical_full_data(i1+1,:)+yshift); hold on; plot((lfp_research_DS(i2,:)-lfp_research_DS(i2+1,:))*yscale, 'r')
subplot(2,1,2);plot(clinical_full_data(i1,:)-clinical_full_data(i1+1,:)+yshift); hold on; plot((lfp_research_DS(12,:)-lfp_research_DS(13,:))*yscale,'r')
linkaxes

%%plot end segment
figure;plot(clinical_full_data(i1,size(clinical_full_data,2)-1000:end)-clinical_full_data(i1+1,size(clinical_full_data,2)-1000:end)+yshift); 
hold on; plot((lfp_research_DS(i2,size(clinical_full_data,2)-1000:end)-lfp_research_DS(i2+1,size(clinical_full_data,2)-1000:end))*yscale, 'r')

%% get onset trial data
pre_stim  = 0.5;
post_stim = 2;
strt_time = on_idx' - pre_stim*fs;
end_time  = on_idx' + post_stim*fs;

% clinical trial data onset
clinical_chan_label=chan_label_clinical;
clinical_trial_data = nan(length(on_idx), length(strt_time(1):end_time(1)),length(chan_label_clinical));
for chan = 1:length(chan_label_clinical)
    for trial = 1:length(on_idx)
        clinical_trial_data(trial,:,chan) = clinical_full_data(chan, strt_time(trial):end_time(trial));
    end
end

% research trial data onset
research_chan_label=chan_label_research;
research_trial_data = nan(length(on_idx), length(strt_time(1):end_time(1)),length(research_chan_label));
for chan = 1:length(chan_label_research)
    for trial = 1:length(on_idx)
        research_trial_data(trial,:,chan) = lfp_research_DS(chan, strt_time(trial):end_time(trial));
    end
end

save('research_trial_data', 'research_trial_data_from_clin')

%%
%cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/clinical_dataset']);
%save('parsed_data', 'chan_label_clinical', 'clinical_full_data')

%% solve discontinuity problem - find disc in pd
% ncs = read_neuralynx_ncs([file_path photodiode_fn]);
% idx = find(lfp.NumValidSamp<512); nValSamp = lfp.NumValidSamp(idx);
% nTrials = length(idx);
% trl = zeros(nTrials, 3);
% trl(1, 1) = 1;
% trl(:, 2)     = ((idx-1)*512 + nValSamp)';
% trl(2:end, 1) = (idx(1:end-1)*512+1)';
% photodiode_signal = ncs.dat(:);
% i = 1; photodiode_1 = photodiode_signal(trl(i,1): trl(i,2));
% i = 2; photodiode_2 = photodiode_signal(trl(i,1): trl(i,2));
% figure; plot(photodiode_signal)
% figure; plot(photodiode_1)
% figure; plot(photodiode_2)