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
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')

% initialize group data
cond1AllSubjERP      = cell(1,length(subj_list));
cond2AllSubjERP      = cell(1,length(subj_list));
cond3AllSubjERP      = cell(1,length(subj_list));
cond4AllSubjERP      = cell(1,length(subj_list));

cond1AllSubjTrialData = cell(1,length(subj_list));
cond2AllSubjTrialData = cell(1,length(subj_list));
cond3AllSubjTrialData = cell(1,length(subj_list));
cond4AllSubjTrialData = cell(1,length(subj_list));

cond1_AbsentERP_AllSubj = cell(1,length(subj_list));
cond2_AbsentERP_AllSubj = cell(1,length(subj_list));
cond3_AbsentERP_AllSubj = cell(1,length(subj_list));
cond4_AbsentERP_AllSubj = cell(1,length(subj_list));

for iSubj = 1:length(subj_list)
    subj          = subj_list{iSubj}
    cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
    
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

    %% get ERP across subjects
    if length(desiredChans)==1
        cond1AllSubjERP{iSubj} = squeeze(nanmean(cond1(:,:,desiredChans),1))';
        cond2AllSubjERP{iSubj} = squeeze(nanmean(cond2(:,:,desiredChans),1))';
        cond3AllSubjERP{iSubj} = squeeze(nanmean(cond3(:,:,desiredChans),1))';
        cond4AllSubjERP{iSubj} = squeeze(nanmean(cond4(:,:,desiredChans),1))';
    else
        cond1AllSubjERP{iSubj} = squeeze(nanmean(cond1(:,:,desiredChans),1));
        cond2AllSubjERP{iSubj} = squeeze(nanmean(cond2(:,:,desiredChans),1));
        cond3AllSubjERP{iSubj} = squeeze(nanmean(cond3(:,:,desiredChans),1));
        cond4AllSubjERP{iSubj} = squeeze(nanmean(cond4(:,:,desiredChans),1));
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
    
    %% get erp removed data
    cond1_AbsentERP = nan(size(cond1,1), size(cond1,2),length(desiredChans));
    cond2_AbsentERP = nan(size(cond2,1), size(cond1,2),length(desiredChans));
    cond3_AbsentERP = nan(size(cond3,1), size(cond1,2),length(desiredChans));
    cond4_AbsentERP = nan(size(cond4,1), size(cond1,2),length(desiredChans));
      for iChan = 1:length(desiredChans)
            cond1_AbsentERP(:,:,iChan)=cond1(:,:,desiredChans(iChan))-nanmean(cond1(:,:,desiredChans(iChan)),1);
            cond2_AbsentERP(:,:,iChan)=cond2(:,:,desiredChans(iChan))-nanmean(cond2(:,:,desiredChans(iChan)),1);
            cond3_AbsentERP(:,:,iChan)=cond3(:,:,desiredChans(iChan))-nanmean(cond3(:,:,desiredChans(iChan)),1);
            cond4_AbsentERP(:,:,iChan)=cond4(:,:,desiredChans(iChan))-nanmean(cond4(:,:,desiredChans(iChan)),1);
      end
        
      cond1_AbsentERP_AllSubj {iSubj} = cond1_AbsentERP;
      cond2_AbsentERP_AllSubj {iSubj} = cond1_AbsentERP;
      cond3_AbsentERP_AllSubj {iSubj} = cond1_AbsentERP;
      cond4_AbsentERP_AllSubj {iSubj} = cond1_AbsentERP;
end

cond2AllSubjERPCat = cat(2,cond2AllSubjERP{:})';
cond3AllSubjERPCat = cat(2,cond3AllSubjERP{:})';


 %% plot ERP for 2 cond
figure;
time = linspace(-pre_stim,post_stim,size(cond3AllSubjERPCat,2))
plot(time, nanmean(cond3AllSubjERPCat,1),'Color','m','LineWidth',2.5)
hold on
plot(time, nanmean(cond2AllSubjERPCat,1),'Color','b','LineWidth',2.5)

xlim([-.2 1])
print('-clipboard','-dbitmap')

 %% plot ERP for both cond with std shades
stdshade(cond2AllSubjERPCat,.1,'b',time,[] ,[], []);hold on;
stdshade(cond3AllSubjERPCat,.1,'m',time,[] ,[], []);
print('-clipboard','-dbitmap')
%% get spectrogram data, either w/ or w/o ERP removal
removeERP     = 'yes'
if strcmp('yes' ,removeERP)
    baseline      = [ '_cond_spec_prestim_cue_responseiveyesremoveERP' removeERP] %'_entire_recording'
else
    baseline      = [ '_cond_spec_prestim_cue_responseiveyes' ] %'_entire_recording'
end

Spect_pre_stim = 0.5; Spect_post_stim = 2.0;
[cond1a,cond2a,cond3a,cond4a] = get_power_subj_elecs(subj_list,exp_type,reg,ref,baseline,lock);
cond2Spectrogram = cat(3,cond2a{:});
cond3Spectrogram = cat(3,cond3a{:});
%% plot spectrogram for 1 exemplar channel
close all
ichan=1
mx = .8
mn = -mx
tickmarks = 1:18:length(freq);

figure;
subplot(2,1,1);
imagesc(linspace(-Spect_pre_stim,Spect_post_stim, size(cond2Spectrogram,2)), 1:length(freq),cond2Spectrogram(:,:,ichan));
hold on; colormap jet; title('lure-')
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),'FontSize', 10, 'FontName', 'Arial');colorbar
ylabel('Frequency (Hz)')
caxis([mn mx])
xlim ([-.2 1])

subplot(2,1,2);
imagesc(linspace(-Spect_pre_stim,Spect_post_stim, size(cond3Spectrogram,2)), 1:length(freq),cond3Spectrogram(:,:,ichan));
hold on; colormap jet; title('lure+')
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
set(gca,'YDir','reverse','YTick',tickmarks,'YTickLabel',round(freq(tickmarks)),'FontSize', 10, 'FontName', 'Arial');colorbar
xlabel('Time (s)'), ylabel('Frequency (Hz)')
caxis([mn mx])
xlim ([-.2 1])
 print('-clipboard','-dbitmap')
 %% for the same chan, plot the ERP, which is the mean across trials in that channel and mean across demeaned trials
figure;
time = linspace(-pre_stim,post_stim,size(cond2AllSubjERPCat,2))
ax(1)=subplot(2,1,1)
plot(time, cond3AllSubjERPCat(ichan,:),'Color','k','LineWidth',1.5)
xlim([-.2 1])
title('Channel ERP')

% get desired chan data, demean it, then average across demeaned trials
temp = cond3AllSubjTrialData{1}(:,:,1);
demeandTrials = temp - cond3AllSubjERPCat(ichan,:);
ax(2)=subplot(2,1,2)
plot(time, nanmean(demeandTrials,1),'Color','k','LineWidth',1.5)
xlim([-.2 1])
title('mean (trial-ERP)')

linkaxes(ax,'x')
print('-clipboard','-dbitmap')

%% visualize indiv trial
figure
for i=1:size(temp,1)
    subplot(3,1,1);
    plot(time, cond3AllSubjERPCat(ichan,:),'k')
    xlim([-.2 1])
    title('ERP') 
    
    subplot(3,1,2);
    plot(time,temp(i,:),'k')
    xlim([-.2 1])
    title('original trial')
    
    subplot(3,1,3);
    plot(time,demeandTrials(i,:),'r')
    title('trial without ERP')
    xlim([-.2 1])
    
    pause
end

%% PSD of given channel
timeIdx = (pre_stim+.2)*fs:(pre_stim+1)*fs; % .2 to 1 sec
subj=7
chan=2
temp3 = cond3AllSubjTrialData{subj}(:,timeIdx,chan);
temp2 = cond2AllSubjTrialData{subj}(:,timeIdx,chan);

% use these lines if subject only had 1 channel and dimensions error out
% temp3 = cond3AllSubjTrialData{subj}(timeIdx,:)';
% temp2 = cond2AllSubjTrialData{subj}(timeIdx,:)';

[Pxx3,F] = pwelch(temp3',size(temp3,2),0,flip(freq'),fs,'power');
[Pxx2,F] = pwelch(temp2',size(temp3,2),0,flip(freq'),fs,'power');

figure
%stdshade(Pxx3',.1,'m',F,[] ,[], []);hold on;
plot(F, nanmean(Pxx3,2),'color', 'm','LineWidth', 2)
hold on;
xlim([2 10])
%stdshade(Pxx2',.1,'b',F,[] ,[], []);hold on;
plot(F, nanmean(Pxx2,2),'color', 'b','LineWidth', 2)

xlim([2 10])
xlabel('freq')
ylabel('power')
title('freq power during 0.2 - 1 sec post stim')
print('-clipboard','-dbitmap')