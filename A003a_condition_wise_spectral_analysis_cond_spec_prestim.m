%% indicate variables
tic
for iSubj = 1:4
    clearvars -except iSubj
    close all; 
   
    subj_list = {'39' '44' '57' '63' }
    %subj_list = {'66' '84' '85' '87'}
    subj = subj_list{iSubj}
    
    % adjust theses
    removeERP      = 'yes';
    cue_responsive = 'yes';
    lock          = 'onset';    % 'response' %'onset'
    exp_type      = 'tuning_correct'; % {'encoding' 'study_test' 'tuning' 'tuning_correct' 'tuning_incorrect' 'indoor_outdoor'}
    ref           = 'LM';
    normalization = 3;          % 1 = entire recording, 2 = pre_stim 3 = cond_spec_prestim
    norm          = {'entire_recording' 'prestim' 'cond_spec_prestim'};
    baseline      = norm{normalization}
    theta_test    = '';
    explore_plots = '';
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

        %% remove erp from trials
    if strcmp('yes', removeERP)
        for iChan = 1:size(cond2,3)
            cond1(:,:,iChan)=cond1(:,:,iChan)-nanmean(cond1(:,:,iChan),1);
            cond2(:,:,iChan)=cond2(:,:,iChan)-nanmean(cond2(:,:,iChan),1);
            cond3(:,:,iChan)=cond3(:,:,iChan)-nanmean(cond3(:,:,iChan),1);
            cond4(:,:,iChan)=cond4(:,:,iChan)-nanmean(cond4(:,:,iChan),1);
        end
    end
%%  
    % get trials
    trial_lengths = [size(cond1,1) size(cond2,1) size(cond3,1) size(cond4,1)];
    
    % get wavelet params
    [freq, scales, dt] = get_freq(fs,minfreq, maxfreq);
    
    % exemplar matrix
    cwt = cwtft({cond1(1, :,1),dt},...
        'scales',scales,'wavelet','morl');
    cwt_power_exemp = 10*log10(abs(cwt.cfs).^2);
    
    % initialize variables
    norm_freq_acrs_chan_cond_1 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond1,1), chan_counter);
    norm_freq_acrs_chan_cond_2 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond2,1), chan_counter);
    norm_freq_acrs_chan_cond_3 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond3,1), chan_counter);
    norm_freq_acrs_chan_cond_4 = zeros(size(cwt_power_exemp,1), size(cwt_power_exemp,2), size(cond4,1), chan_counter);

    for chan =1:chan_counter
        % first condition
        if strcmp('encoding',exp_type)
            desired_chan_powr_mn = chan_powr_mn3a;
            desired_chan_powr_std= chan_powr_std3a;
            cond_prestims_trial_artif = cond3a_prestims_trial_artif;
        elseif strcmp('tuning_correct',exp_type) 
            desired_chan_powr_mn = chan_powr_mn1;
            desired_chan_powr_std= chan_powr_std1;
            cond_prestims_trial_artif = cond1_prestims_trial_artif;
        end
        
        %% first condition
        parfor trial = 1:trial_lengths(1)
           
            trial_temp = cond1(trial, :,chan);
            % wavelet spectrogram
            cwt = cwtft({trial_temp,dt},...
                'scales',scales,'wavelet','morl');
            cwt_power = 10*log10(abs(cwt.cfs).^2);
             norm_freq = nan(size(cwt_power));
            % normalize z-score relative to 1-hour session
           if cond_prestims_trial_artif(chan,trial)==0
            for a = 1:size(cwt_power,1) %loop thru freq
                for b = 1:size(cwt_power,2)%loop thru timepoints
                    norm_freq(a, b) = (cwt_power(a,b)-desired_chan_powr_mn(a,chan))/desired_chan_powr_std(a, chan);
                end
            end
           end
            norm_freq_acrs_chan_cond_1(:, :, trial, chan) =  norm_freq;
        end
        clear desired_chan_powr_mn desired_chan_powr_std cond_prestims_trial_artif

        %% second condition
        if strcmp('encoding',exp_type)
            desired_chan_powr_mn = chan_powr_mn2a;
            desired_chan_powr_std= chan_powr_std2a;
            cond_prestims_trial_artif = cond2a_prestims_trial_artif;
        elseif strcmp('tuning_correct',exp_type) 
            desired_chan_powr_mn = chan_powr_mn2;
            desired_chan_powr_std= chan_powr_std2;
            cond_prestims_trial_artif = cond2_prestims_trial_artif;
        end
        parfor trial = 1:trial_lengths(2)
            trial_temp = cond2(trial, :,chan);
            % wavelet spectrogram
            cwt = cwtft({trial_temp,dt},...
                'scales',scales,'wavelet','morl');
            cwt_power = 10*log10(abs(cwt.cfs).^2);
            % normalize z-score relative to 1-hour session
             norm_freq = nan(size(cwt_power));
            % normalize z-score relative to 1-hour session
           if cond_prestims_trial_artif(chan,trial)==0
            for a = 1:size(cwt_power,1) %loop thru freq
                for b = 1:size(cwt_power,2)%loop thru timepoints
                    norm_freq(a, b) = (cwt_power(a,b)-desired_chan_powr_mn(a,chan))/desired_chan_powr_std(a, chan);
                end
            end
           end
            norm_freq_acrs_chan_cond_2(:, :, trial, chan) =  norm_freq;
        end
        clear desired_chan_powr_mn desired_chan_powr_std cond_prestims_trial_artif        
        %% third condition
        if strcmp('tuning_correct',exp_type)
            desired_chan_powr_mn = chan_powr_mn3;
            desired_chan_powr_std= chan_powr_std3;
            cond_prestims_trial_artif = cond3_prestims_trial_artif;
        end
        parfor trial = 1:trial_lengths(3)
            trial_temp = cond3(trial, :,chan);
            % wavelet spectrogram
            cwt = cwtft({trial_temp,dt},...
                'scales',scales,'wavelet','morl');
            cwt_power = 10*log10(abs(cwt.cfs).^2);
            % normalize z-score relative to 1-hour session
             norm_freq = nan(size(cwt_power));
            % normalize z-score relative to 1-hour session
           if cond_prestims_trial_artif(chan,trial)==0
            for a = 1:size(cwt_power,1) %loop thru freq
                for b = 1:size(cwt_power,2)%loop thru timepoints
                    norm_freq(a, b) = (cwt_power(a,b)-desired_chan_powr_mn(a,chan))/desired_chan_powr_std(a, chan);
                end
            end
           end
            norm_freq_acrs_chan_cond_3(:, :, trial, chan) =  norm_freq;
        end
        clear desired_chan_powr_mn desired_chan_powr_std cond_prestims_trial_artif
        
        %% fourth condition
        if strcmp('tuning_correct',exp_type)
            desired_chan_powr_mn = chan_powr_mn4;
            desired_chan_powr_std= chan_powr_std4;
            cond_prestims_trial_artif = cond4_prestims_trial_artif;
        end
          parfor trial = 1:trial_lengths(4)
            trial_temp = cond4(trial, :,chan);
            % wavelet spectrogram
            cwt = cwtft({trial_temp,dt},...
                'scales',scales,'wavelet','morl');
            cwt_power = 10*log10(abs(cwt.cfs).^2);
            % normalize z-score relative to 1-hour session
             norm_freq = nan(size(cwt_power));
            % normalize z-score relative to 1-hour session
           if cond_prestims_trial_artif(chan,trial)==0
            for a = 1:size(cwt_power,1) %loop thru freq
                for b = 1:size(cwt_power,2)%loop thru timepoints
                    norm_freq(a, b) = (cwt_power(a,b)-desired_chan_powr_mn(a,chan))/desired_chan_powr_std(a, chan);
                end
            end
           end
            norm_freq_acrs_chan_cond_4(:, :, trial, chan) =  norm_freq;
        end
        clear desired_chan_powr_mn desired_chan_powr_std cond_prestims_trial_artif      
        
 end
    
    
    % average across trials within chans
    edge_points = 200;
    
    % initialize matrix for all frequencies
    mn_acrs_trials_cond1 = nan(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2)-2*edge_points, chan_counter);
    mn_acrs_trials_cond2 = nan(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2)-2*edge_points, chan_counter);
    mn_acrs_trials_cond3 = nan(size(norm_freq_acrs_chan_cond_1,1), size(norm_freq_acrs_chan_cond_1,2)-2*edge_points, chan_counter);
    mn_acrs_trials_cond4 = nan(size(norm_freq_acrs_chan_cond_2,1), size(norm_freq_acrs_chan_cond_2,2)-2*edge_points, chan_counter);
    
    % remove edges
    norm_freq_acrs_chan_cond_1 = norm_freq_acrs_chan_cond_1(:, edge_points+1:end-edge_points, :, :);
    norm_freq_acrs_chan_cond_2 = norm_freq_acrs_chan_cond_2(:, edge_points+1:end-edge_points, :, :);
    norm_freq_acrs_chan_cond_3 = norm_freq_acrs_chan_cond_3(:, edge_points+1:end-edge_points, :, :);
    norm_freq_acrs_chan_cond_4 = norm_freq_acrs_chan_cond_4(:, edge_points+1:end-edge_points, :, :);
    
    % average across trials within chans
    for elec = 1:size(norm_freq_acrs_chan_cond_1,4) %loop thru chans
        mn_acrs_trials_cond1 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_1(:, :, :, elec), 3);
        mn_acrs_trials_cond2 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_2(:, :, :, elec), 3);
        mn_acrs_trials_cond3 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_3(:, :, :, elec), 3);
        mn_acrs_trials_cond4 (:, :, elec) = nanmean(norm_freq_acrs_chan_cond_4(:, :, :, elec), 3);
    end
    
    % mean acrs chans
    indiv_freq_chans_cond1 = nanmean(mn_acrs_trials_cond1,3);
    indiv_freq_chans_cond2 = nanmean(mn_acrs_trials_cond2,3);
    indiv_freq_chans_cond3 = nanmean(mn_acrs_trials_cond3,3);
    indiv_freq_chans_cond4 = nanmean(mn_acrs_trials_cond4,3);
    
    
    %% figure params
    tickmarks = 1:20:length(freq);
    
    % get start and end trial indices
    if strcmp ('onset',lock)
        if strcmp('yes',theta_test)
            disp('NEED TO CHECK U HAVE DATA FOR THIS TIME RANGE')
            pre_stim  = 2;
        else
            pre_stim  = 0.9;
        end
        post_stim = 2.4; % changed to 1 sec. ISI is 1 sec?
        strt_time = on_idx' - pre_stim*fs;
        end_time  = on_idx' + post_stim*fs;
        
    elseif strcmp ('offset',lock)
        pre_stim  = 1;
        post_stim = .5;
        strt_time = off_idx' - pre_stim*fs
        end_time  = off_idx' + post_stim*fs;
    end
    save_individual_subj_region_data(cue_responsive, ref,baseline, exp_type, lock, subj,removeERP,...
        norm_freq_acrs_chan_cond_1, norm_freq_acrs_chan_cond_2,norm_freq_acrs_chan_cond_3,norm_freq_acrs_chan_cond_4)
    
    save_data_for_average_spectrograms(cue_responsive, ref, baseline,exp_type, lock, subj,removeERP, mn_acrs_trials_cond1,mn_acrs_trials_cond2,...
        mn_acrs_trials_cond3, mn_acrs_trials_cond4)
    
%     save_data_for_average_spectrograms_all_trials_elecs(cue_responsive, ref, baseline, exp_type, lock, subj, norm_freq_acrs_chan_cond_1,norm_freq_acrs_chan_cond_2,...
%         norm_freq_acrs_chan_cond_3, norm_freq_acrs_chan_cond_4)
    toc

end
disp('done')
%%

reg = 'FRO'
cond_num = sum(trial_lengths>0)
[OFC_chan_idx,fro_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
    ACC_chan_idx,EC_chan_idx, HC_chan_idx ,CA3_chan_idx,CA1_chan_idx, MTL_chan_idx, NC_chan_idx]  = get_elecs_clean_cue_resp(subj)
if strcmp('OFC',reg)
    chan_idx = OFC_chan_idx;
elseif strcmp('FRO',reg)
    chan_idx = fro_chan_idx;
elseif strcmp('TEMP',reg)
    chan_idx = temp_chan_idx;
elseif strcmp('CING',reg)
    chan_idx = cingulate_chan_idx;
elseif strcmp('ACC',reg)
    chan_idx = ACC_chan_idx;
elseif strcmp('INS',reg)
    chan_idx = insula_chan_idx;
elseif strcmp('EC',reg)
    chan_idx = EC_chan_idx;
elseif strcmp('CA1',reg)
    chan_idx = CA1_chan_idx;
elseif strcmp('CA3',reg)
    chan_idx = CA3_chan_idx;
elseif strcmp('HC',reg)
    chan_idx = HC_chan_idx;
elseif strcmp('MTL',reg)
    chan_idx = MTL_chan_idx;
end
%%
if strcmp('yes', explore_plots);
    %% ENCODING REGION PLOTS
    cd('/mnt/yassamri/iEEG/sandra/IndivSubjFig/spectrograms')
    reg= 'CA3'
    mn = -0.5
    mx = 0.5
    
    reg_cond1 = mn_acrs_trials_cond1(:,:,chan_idx);
    reg_cond2 = mn_acrs_trials_cond2(:,:,chan_idx);
    reg_cond3 = mn_acrs_trials_cond3(:,:,chan_idx);
    reg_cond4 = mn_acrs_trials_cond4(:,:,chan_idx);
    
    FigH = figure %('Position', get(0, 'Screensize'));
    
    
    % plot first condition
    subplot (cond_num,1,1)
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), nanmean(reg_cond1,3))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    if strcmp('encoding',exp_type)
        title('--> lure +')
    else
        title('repeat')
    end
    colormap jet
    
    subplot (cond_num,1,2)
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), nanmean(reg_cond2,3))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    if strcmp('encoding',exp_type)
        title ('--> lure -')
    else
        title('repeat')
    end
    colormap jet
    suptitle(reg)
    
    
    FigH.PaperUnits    = 'inches';
    %FigH.PaperPosition = [20 20 15 26];
    print([reg '_' subj '.png'],'-dpng','-r0')
    %% RETREIVAL REGION PLOTS
    % close all
    mn  = -0.3
    mx  = 0.3
   
    reg_cond1 = mn_acrs_trials_cond1(:,:,chan_idx);
    reg_cond2 = mn_acrs_trials_cond2(:,:,chan_idx);
    reg_cond3 = mn_acrs_trials_cond3(:,:,chan_idx);
    reg_cond4 = mn_acrs_trials_cond4(:,:,chan_idx);
    
    
    FigH = figure %('Position', get(0, 'Screensize'));
    % plot first condition
    subplot (cond_num,1,1)
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), nanmean(reg_cond1,3))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('repeat')
    colormap jet
    
    subplot (cond_num,1,2)
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), nanmean(reg_cond2,3))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('lure -')
    colormap jet
    
    subplot (cond_num,1,3)
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), nanmean(reg_cond3,3))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('lure +')
    colormap jet
    
    subplot (cond_num,1,4)
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), nanmean(reg_cond4,3))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('new')
    colormap jet
    suptitle(reg)
    
    FigH.PaperUnits = 'inches';
    FigH.PaperPosition = [20 20 15 26];
    print([reg '_' subj '.png'],'-dpng','-r0')
    %%  spectrogram indiv chans
    for chan = 74
        mn   = -.75
        mx   = .75
        
        figure
        % plot first condition
        subplot (cond_num,1,1)
        imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond1(:,:,chan))
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mn mx])
        xlabel('time')
        ylabel('freq')
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title('repeat')
        colormap jet
        
        subplot (cond_num,1,2)
        imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond2(:,:,chan))
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mn mx])
        xlabel('time')
        ylabel('freq')
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title('Pattern Comp: incorr lure')
        colormap jet
        
        subplot (cond_num,1,3)
        imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond3(:,:,chan))
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mn mx])
        xlabel('time')
        ylabel('freq')
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title('Pattern Sep: corr lure')
        colormap jet
        
        subplot (cond_num,1,4)
        imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), mn_acrs_trials_cond4(:,:,chan))
        set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
        colorbar
        caxis([mn mx])
        xlabel('time')
        ylabel('freq')
        hold on
        plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
        title('new')
        colormap jet
        suptitle(chan_label{chan})
    end
    
    %% bar plot test
    figure
    bar([mean(mean(mn_acrs_trials_cond1(freq<6,400:1100,chan),2)) mean(mean(mn_acrs_trials_cond2(freq<6,400:1100,chan),2))...
        mean(mean(mn_acrs_trials_cond3(freq<6,400:1100,chan),2)) mean(mean(mn_acrs_trials_cond4(freq<6,400:1100,chan),2))])
    
    % downsamples
    figure
    bar([mean(mean(mn_acrs_trials_cond1(freq<6,150:850,chan),2)) mean(mean(mn_acrs_trials_cond2(freq<6,150:850,chan),2))...
        mean(mean(mn_acrs_trials_cond3(freq<6,150:850,chan),2)) mean(mean(mn_acrs_trials_cond4(freq<6,150:850,chan),2))])
    
    
    
    %% subtraction
    chan_idx
    chan =124;
    mn = -.8;
    mx = .8;
    subt = mn_acrs_trials_cond1(:,:,chan)-mn_acrs_trials_cond2(:,:,chan);
    figure
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), subt)
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('lure+ - new')
    colormap jet
    %%
    subt = mn_acrs_trials_cond3(:,:,chan)-mn_acrs_trials_cond2(:,:,chan);
    figure
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), subt)
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
    title('lure+ - lure-')
    colormap jet
    
    %% check ERP
    figure
    cntr = 0
    for i = 21:30
        
        cntr = cntr+1
        subplot(10, 1, cntr);plot(linspace(-pre_stim,post_stim,size(cond3,2)),nanmean(trial_data(:,:,i),1))
    end
    figure;plot(linspace(-pre_stim,post_stim,size(cond3,2)),nanmean(nanmean(trial_data,3),1))
    %% plot trials - cond 1-4
    chan  =chan_idx(1)
    mult1 = 1/5
    mult2 = 1/5
    cond_mtx = cond1;
    figure
    for i=1:10
        plot(linspace(-pre_stim,post_stim,size(cond1,2)),cond_mtx(i,:,chan))
        hold on
        %pause
        %title([num2str(i) ' repeat'])
    end
    plot(zeros(1,length(mult1*min(min(cond_mtx(:,:,chan))):max(max(cond_mtx(:,:,chan)))*mult2)), mult1*min(min(cond_mtx(:,:,chan))):max(max(cond_mtx(:,:,chan)))*mult2,'Color','k','LineWidth',1)
    
    %% plot spectrogram of single trial
    mn = -2
    mx = 2
    figure
    imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), norm_freq(:, edge_points+1:end-edge_points))
    set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
    colorbar
    colormap jet
    caxis([mn mx])
    xlabel('time')
    ylabel('freq')
    hold on
    plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
end