% description:
% 1) plot ERP for both conditions in a desired region, 
% 2) plot region specific spectrograms w/ or w/o ERP removal
% 3) plot ERP for a given chan, and mean(trials-ERP) in the same chan
% 4) plot idividual trials and trials-erp
% 5) individual channel PSD

close all;
reg = 'CA3'
cue_responsive = 'yes'
subj_list     = {'39' '44' '57' '63' '66' '84' '85' '87'}
lock          = 'onset';    % 'response' %'onset'
exp_type      = 'tuning_correct'; % {'encoding' 'study_test' 'tuning' 'tuning_correct' 'tuning_incorrect' 'indoor_outdoor'}
ref           = 'LM';
DS            = 'yes';
clinical      = 'yes';
select_chan   = 3;
addpath('/tmp/yassamri/iEEG/sandra/analysis_pipeline_final')

% initialize group data
cond1AllSubjTrialData = cell(1,length(subj_list));
cond2AllSubjTrialData = cell(1,length(subj_list));
cond3AllSubjTrialData = cell(1,length(subj_list));
cond4AllSubjTrialData = cell(1,length(subj_list));


for iSubj = 1:length(subj_list)
    subj          = subj_list{iSubj}
    cd(['/tmp/yassamri/iEEG/sandra/subj_' subj])
    
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
    
    
    % load trial data
    if strcmp('yes',DS)
        load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3' fn_nm '_fs_500.mat'])
    elseif strcmp('',DS)
        load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3_NotDownsampled.mat'])
    end
    [cond1,cond2,cond3,cond4,cond5,cond6] = GetCondData(subj, exp_type, lock, DS, fn_nm, ref);
    
    %% index region specific and cue responsive elecs
    if strcmp('yes',cue_responsive)
        [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
            ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj);
    else
        [OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
            ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean(subj);
    end
    if strcmp('HC',reg)
        desiredChans = HC_chan_idx;
    elseif strcmp('CA3',reg)
        desiredChans = CA3_chan_idx;
    elseif strcmp('CA1',reg)
        desiredChans = CA1_chan_idx;
    elseif strcmp('OFC',reg)
        desiredChans = OFC_chan_idx;
    elseif strcmp('FRO',reg)
        desiredChans = fro_chan_idx;
    elseif strcmp('TEMP',reg)
        desiredChans = temp_chan_idx;
    elseif strcmp('CING',reg)
        desiredChans = cingulate_chan_idx;
    elseif strcmp('INS',reg)
        desiredChans = insula_chan_idx;
    elseif strcmp('EC',reg)
        desiredChans = EC_chan_idx;
    end

  
    
    %% get cond data across subjects
    if length(desiredChans)==1
        cond1AllSubjTrialData{iSubj} = cond1(:,:,desiredChans)';
        cond2AllSubjTrialData{iSubj} = cond2(:,:,desiredChans)';
        cond3AllSubjTrialData{iSubj} = cond3(:,:,desiredChans)';
        cond4AllSubjTrialData{iSubj} = cond4(:,:,desiredChans)';
    else
        cond1AllSubjTrialData{iSubj} = cond1(:,:,desiredChans);
        cond2AllSubjTrialData{iSubj} = cond2(:,:,desiredChans);
        cond3AllSubjTrialData{iSubj} = cond3(:,:,desiredChans);
        cond4AllSubjTrialData{iSubj} = cond4(:,:,desiredChans);
    end
    
    
end

cond2AllSubjERPCat = cat(2,cond2AllSubjERP{:})';
cond3AllSubjERPCat = cat(2,cond3AllSubjERP{:})';

%% bandpass filter and plot
set(0,'DefaultFigureWindowStyle','normal')

fpass = [4 6];
filteredData3 = (bandpass(temp3', fpass,fs))';
filteredData2 = (bandpass(temp2', fpass,fs))';
 
time = linspace(-pre_stim,post_stim,size(cond3AllSubjERPCat,2));

figure; 
subplot(2,1,1);stdshade(filteredData3,.1,'m',time,[] ,[], []);
xlim([-.2 1])
subplot(2,1,2);stdshade(filteredData2,.1,'b',time,[] ,[], []);
xlim([-.2 1])

print('-clipboard','-dbitmap')


%% get power and smooth it to be able to average across trials
iChan=1
iSubj=5
temp3 = cond3AllSubjTrialData{iSubj}(:,:)';
temp2 = cond2AllSubjTrialData{iSubj}(:,:)';

power3 = temp3.^2;
powerConv3 = nan(size(power3));
for iTrial = 1:size(power3,1)
    powerConv3(iTrial,:)=conv(power3(iTrial,:),[1 1 1 1 1],'same');
end

power2 = temp2.^2;
powerConv2 = nan(size(power2));
for iTrial = 1:size(power2,1)
    powerConv2(iTrial,:)=conv(power2(iTrial,:),[1 1 1 1 1],'same');
end

figure;subplot (3,1,1);hold on
stdshade(powerConv3,.1,'m',time,[] ,[], []);hold on;xlim([-.2 1])
stdshade(powerConv2,.1,'b',time,[] ,[], []);hold on;xlim([-.2 1])
print('-clipboard','-dbitmap')
subplot (3,1,2);
stdshade(powerConv3,.1,'m',time,[] ,[], []);hold on;xlim([-.2 1])
subplot (3,1,3);
stdshade(powerConv2,.1,'b',time,[] ,[], []);hold on;xlim([-.2 1])
print('-clipboard','-dbitmap')
%% plot bandpass filtered individual trials 
%data1 = temp3;
data1 = temp2;

data2 = filteredData3;
%data2 = filteredData2;
figure
for i=1:size(data1,1)

%     
%     subplot(1,2,1);
%     plot(time,data1(i,:),'k')
%     xlim([-.2 1])
%     title('original trial')
%     
%     subplot(1,2,2);
    plot(time,data2(i,:),'k')
    xlim([-.2 1])
    title('original trial Filtered')
    
%     subplot(2,2,3);
%     
%     subplot(2,2,3);
%     plot(time,demeandTrials(i,:),'r')
%     title('trial without ERP')
%     xlim([-.2 1])
    
    %pause
    hold on
end


