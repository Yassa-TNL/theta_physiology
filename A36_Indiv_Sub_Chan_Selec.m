clear all;close all;clc
subj_list     = {'57'}
for iSubj = 1:length(subj_list)
    subj = subj_list{iSubj}
    exp_type      = 'tuning_correct' %encoding, tuning_correct   = 
    reg_list      = {'HC'} %{'OFC' 'FRO' 'TEMP' 'CING' 'INS' 'EC' 'HC' 'CA3'}
    lock          = 'onset'   %onset response
    ref           = 'LM'
    fs            = 500;
    baseline      = '_cond_spec_prestim' % '': entire recording, '_prestim'
    desired_freq  = 'deltatheta'
    desired_freq_lo  = 2
    desired_freq_hi  = 5
    plots_4_conds    =''
    trial_count_limit = ''
    fn_ext        = '_cue_responseive_yes'
    addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
    if strcmp('encoding', exp_type);
        cond_num=2;    titles = {'lure+','lure-'}
        xtick_vals = [-0.5 0 2];
    else cond_num = 4; titles = {'repeat', 'lure-', 'lure+', 'new'}
    end
    
    if strcmp('response', lock)
        pre_stim           = 1.5;
        post_stim          = .5;
        pre_stim_minus_edge  = 1.5-200/fs;
        post_stim_minus_edge = 0.5-200/fs;
        time_range = [0.1*fs 1.2*fs];
        current_pre_stim   = -1;
        current_post_stim  = 0.1;
        xtick_vals         = [-1 0 0.1];        
    elseif strcmp('onset', lock)
        if strcmp('encoding', exp_type)
            pre_stim = 0.5; post_stim = 2;
            xtick_vals = [-pre_stim 0 1 2];
        else
            pre_stim = 0.5; post_stim = 1
            xtick_vals = [-0.2 0 1];
            time_range = [1 (pre_stim+1)*fs]
        end
    elseif strcmp('onset2', lock)
        pre_stim = 0.5; post_stim = 2;
        time_range = [.3*fs (pre_stim+2)*fs]; % was [301 1301]
        xtick_vals = [current_pre_stim 0 1 2];
    end
    minfreq = 3; maxfreq = 200;
    [freq]  = get_freq(fs,minfreq, maxfreq);
    tickmarks = 1:21:length(freq);
    all_chan_sig = zeros(1,length(reg_list));
    cd(['/mnt/yassamri/iEEG/sandra/group_data/IndivSubjData/' ref '_reref/normalization' baseline  '/' exp_type '_' lock '/' subj])
    trial_nums=[];

    for iReg = 1:length(reg_list)
        reg           = reg_list{iReg}
        load(['subj' subj '_' reg '_spectrograms_cond_spec_prestim' fn_ext '.mat'])
        
        chan_sig = zeros(1,size(cond2,4))
        for iElec = 1:size(cond2,4)
            trace2 = squeeze(nanmean(cond2((freq>desired_freq_lo)&(freq<desired_freq_hi), :, :,iElec),1))'; % trial by time
            if strcmp('encoding',exp_type)
                trace3 = squeeze(nanmean(cond1((freq>desired_freq_lo)&(freq<desired_freq_hi), :, :,iElec),1))';
            else
                trace3 = squeeze(nanmean(cond3((freq>desired_freq_lo)&(freq<desired_freq_hi), :, :,iElec),1))';
            end
            
            % remove nans
            trace2 = trace2(~isnan(trace2(:,1)),:);
            trace3 = trace3(~isnan(trace3(:,1)),:);
          %  if strcmp('yes',trial_count_limit)
                %if any([size(trace3,1) size(trace2,1)]<15)
                    %continue
                %end
           % end
           trial_numbers(iElec,:)= [size(trace3,1) size(trace2,1)];
            % subtract to reset
            reset = nanmean(nanmean(trace2(:,0.25*fs:pre_stim*fs),1));
            trace2 = trace2-reset;
            reset = nanmean(nanmean(trace3(:,0.25*fs:pre_stim*fs),1));
            trace3 = trace3-reset;
            
            % get time range of interest
            if strcmp('retrieval', exp_type)
                trace2 = trace2(:,time_range(1):time_range(2));
                trace3 = trace3(:,time_range(1):time_range(2));
            end
            % perm testing
            [zmap,zmapthresh,zmapthresh_for_plot_temp] = permutation_testing_vector(trace3(:,pre_stim*fs+1:end),trace2(:,pre_stim*fs+1:end), 500);
            zmapthresh_for_plot = nan(1,size(trace3,2));
            zmapthresh_for_plot(pre_stim*fs+1:end) = zmapthresh_for_plot_temp;
            chan_sig (iElec) =  nansum(zmapthresh_for_plot_temp)>0;
            
            % plot
            h1 = figure;%set(gcf,'Visible', 'off');
            hold on;
            % 2 conds
            stdshade(trace3,.1,'m',linspace(-pre_stim,post_stim, size(trace3,2)),[] ,[], []);
            h1 = plot(linspace(-pre_stim,post_stim, size(trace3,2)),nanmean(trace3,1), 'm', 'LineWidth', 2);
            stdshade(trace2,.1,'b',linspace(-pre_stim,post_stim, size(trace2,2)),[] ,[], []);
            h2 = plot(linspace(-pre_stim,post_stim, size(trace2,2)),nanmean(trace2,1), 'b', 'LineWidth', 2);
            
            % overlay stas
            scale = max([nanmean(trace3,1) nanmean(trace2,1)])+0.02;
            plot(linspace(-pre_stim,post_stim, size(trace2,2)),scale*zmapthresh_for_plot, 'LineWidth', 2, 'Color', 'k')
            
            % fig prop
            xlim([-pre_stim post_stim]);
            xticks(xtick_vals)
            y=ylim;
            if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
            elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
            end
            xlabel('Time (s)')
            title([reg ' ' exp_type ' subj ' subj])
            line([0 0],y,'color','k')
            %print('-clipboard','-dbitmap')
            saveas(gcf, ['000_cue_resp_method2' reg '_elec' num2str(iElec) '.bmp'])
          %  close all
        end
    end
    
    % plot traces and stats
    if strcmp('tuning_correct',exp_type)&&strcmp('yes',plots_4_conds)
        stdshade(trace1,.1,'r',linspace(current_pre_stim,current_post_stim, size(trace1,2)),[] ,[], []);
        h3 = plot(linspace(current_pre_stim,current_post_stim, size(trace3,2)),nanmean(trace1,1), 'r', 'LineWidth', 2);
        stdshade(trace4,.1,'k',linspace(current_pre_stim,current_post_stim, size(trace4,2)),[] ,[], []);
        h4 = plot(linspace(current_pre_stim,current_post_stim, size(trace2,2)),nanmean(trace4,1), 'k', 'LineWidth', 2);
    end
    all_chan_sig (iReg) = sum(chan_sig);
    
end
