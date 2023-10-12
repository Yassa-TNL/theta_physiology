% plots chanel raw data, averaging across trials for encoding, either cond1
% or cond2
close all; clc
subj      = '39'
cond_desired = 2
lock          = 'onset'    % 'response' %'onset'
exp_type      = 'encoding' % {'encoding' 'study_test' 'tuning' 'tuning_correct' 'tuning_incorrect' 'indoor_outdoor'}
ref           = 'LM';
normalization = 3          % 1 = entire recording, 2 = pre_stim 3 = cond_spec_prestim
norm          = {'entire_recording' 'prestim' 'cond_spec_prestim'}
baseline      = norm{normalization}
DS            = 'yes';
clinical      = 'yes';
freq_analysis = 'wavelet';
select_chan   = 3;
minfreq       = 3;
maxfreq       = 200;
fs            = 500;
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')

% indicate whether to pull from clin or research recordings
if strcmp('84', subj) || strcmp('85', subj) || strcmp('87', subj)
    if strcmp('yes', clinical)
        fn_nm = '_clinical';
    else
        fn_nm = '_research';
    end
else
    fn_nm = '';
end

% load baseline
if strcmp('39', subj) || strcmp('44', subj) || strcmp('57', subj) || strcmp('63', subj) || strcmp('66', subj) || strcmp('83', subj)
    load(['normalization_' norm{normalization} '_ref_' ref '_baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_' subj '_select_chan_3_fs_500.mat'])
elseif strcmp('84', subj) || strcmp('85', subj) || strcmp('87', subj)
    load(['normalization_' norm{normalization} '_ref_' ref '_baseline_info_wavelet_32num_20logdb_3hz_200hz_notched_artifact_reject_subj_' subj '_select_chan_3' fn_nm '_fs_500.mat'])
    % load('timestamps.mat')
end
chan_counter = size(chan_powr_mn,2);

% load epoched data and group it into conds
if strcmp('yes',DS)
    load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3' fn_nm '_fs_500.mat'])
elseif strcmp('',DS)
    load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3_NotDownsampled.mat'])
end
[cond1,cond2,cond3,cond4,cond5,cond6] = GetCondData(subj, exp_type, lock, DS, fn_nm, ref);
trial_lengths = [size(cond1,1) size(cond2,1) size(cond3,1) size(cond4,1)]
%
reg = 'HC'
[OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,ACC_chan_idx,EC_chan_idx,...
    CA3_chan_idx,CA1_chan_idx,HC_chan_idx,MTL_chan_idx,NC_chan_idx] = get_elecs_clean(subj);
if strcmp('OFC',reg); chan_idx = OFC_chan_idx;
elseif strcmp('FRO',reg);chan_idx = fro_chan_idx;
elseif strcmp('TEMP',reg);chan_idx = temp_chan_idx;
elseif strcmp('CING',reg); chan_idx = cingulate_chan_idx;
elseif strcmp('ACC',reg);chan_idx = ACC_chan_idx;
elseif strcmp('INS',reg);chan_idx = insula_chan_idx;
elseif strcmp('EC',reg);chan_idx = EC_chan_idx;
elseif strcmp('CA1',reg);chan_idx = CA1_chan_idx;
elseif strcmp('CA3',reg);chan_idx = CA3_chan_idx;
elseif strcmp('HC',reg);chan_idx = HC_chan_idx;
elseif strcmp('MTL',reg);chan_idx = MTL_chan_idx;
end

if cond_desired==1
    artifact_free_cond_matrix = nan(size(cond1));
    artf_mtx = cond3a_prestims_trial_artif;
    cond_mtx =cond1;
elseif cond_desired==2
    artifact_free_cond_matrix = nan(size(cond2));
    artf_mtx = cond2a_prestims_trial_artif;
    cond_mtx =cond2;
end
for iChan = 1:chan_counter
    for iTrial = 1:size(cond_mtx,1)
        if artf_mtx(iChan, iTrial)==0
        artifact_free_cond_matrix(iTrial,:,iChan) = cond_mtx(iTrial,:,iChan);
        end
    end
end
% check ERP
artifact_free_cond_matrix;
figure
hold on
for iChan = chan_idx
    plot(linspace(-pre_stim,post_stim,size(artifact_free_cond_matrix,2)),nanmean(artifact_free_cond_matrix(:,:,iChan),1))
    
end
legend
y=ylim;
line([0 0], [y(1) y(end)])
print('-clipboard','-dbitmap')

