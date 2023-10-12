clear all;close all; clc;
tic
disp('CFC')
subj = '39'
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
[fro_chan_idx,MTL_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,OFC_chan_idx,CA3_chan_idx,CA1_chan_idx,HC_chan_idx] = get_elecs(subj);
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])  
load('cond_data_onset_tuning_correct.mat')

if strcmp('39',subj) || strcmp('57',subj) || strcmp('66',subj) || strcmp('44',subj)
    phase_regions_list = {'OFC', 'FRO', 'TEMP', 'CING'};
    power_regions_list = {'MTL', 'HC'};
elseif strcmp('63',subj) || strcmp('84',subj)
    phase_regions_list = {'OFC', 'FRO', 'TEMP', 'CING'};
    power_regions_list = {'MTL'}; 
end


for reg1 = 1:length(phase_regions_list)
    phase_region = phase_regions_list{reg1}
    
    for reg2 = 1:length(power_regions_list)
        power_region = power_regions_list{reg2}
        clear phase_elecs power_elecs
        % get elecs for phase
            if strcmp('OFC',phase_region)
                phase_elecs = OFC_chan_idx;
                
            elseif strcmp('FRO',phase_region)
                phase_elecs = fro_chan_idx;
                
            elseif strcmp('TEMP',phase_region)
                phase_elecs = temp_chan_idx;
                
            elseif strcmp('CING',phase_region)
                phase_elecs = cingulate_chan_idx;
            end
            
            % get elecs for power
            if strcmp('MTL',power_region)
                power_elecs = MTL_chan_idx;
                
            elseif strcmp('HC',power_region)
                power_elecs = HC_chan_idx;
                
            elseif strcmp('CA3',power_region)
                power_elecs = CA3_chan_idx;
                
            elseif strcmp('CA1',power_region)
                power_elecs = CA1_chan_idx;
                
            elseif strcmp('AMY',power_region)
                power_elecs = AMY_chan_idx;
            end
            
            if isempty (power_elecs) | isempty(phase_elecs)
            
                continue 
            end   
         
        disp([reg1 reg2]) 
        for cond = 1:4
            clear data cond_data 

            % get cond data
            if cond==1
                cond_data = cond1;
            elseif  cond==2
                cond_data = cond2;
            elseif  cond==3
                cond_data = cond3;
            elseif  cond==4
                cond_data = cond4;
            end
            
            time_idx           = [500 1500]; % 1 sec after stim onset
            data = [];
            for elec = 1:size(cond_data,3) % loop thru elecs
                temp = cond_data(:,time_idx(1):time_idx(2),elec);
                temp2 = temp';
                data(elec,:) = temp2(:);
            end
            
            data_for_phase = data(phase_elecs,:);
            data_for_power = data(power_elecs,:);
            
            % Wavelet params
            fs    = 1000; % sampling rate in Hz
            time  = -1:1/fs:1; % time, from -1 to 1 second in steps of 1/sampling-rate
            wavelet_length = length(time);
            half_wavelet   = (wavelet_length-1)/2;
            nconv          = (wavelet_length+size(data, 2)-1);
            
            % Create wavelets to extract low freqs
            phas_freqs  = 2.^(1.56:0.3:5.15); % frequency of wavelet in Hz (2.9485 : 29.0406)
            low_freq_wavelet=zeros(length(phas_freqs), wavelet_length);
            lo_freq_wavelet_fft=[];
            for f=1:length(phas_freqs)
                s  = 3/(2*pi*phas_freqs(f)); %width of wavelet
                low_freq_wavelet(f, :) = exp(2*pi*1i*phas_freqs(f).*time) .* exp(-time.^2./(2*(s)^2)); %eulers form*2pift*gaussian
                lo_freq_wavelet_fft(f, :)=fft((low_freq_wavelet(f, :))',nconv); %fft of the wavelets
                clear s
            end
            
            % Wavelets to extract hi freqs
            power_freqs = 2.^(5.3:0.15:8); % frequency of wavelet Hz
            hi_freq_wavelet=zeros(length(power_freqs), wavelet_length);
            hi_freq_wavelet_fft=zeros(length(power_freqs), nconv);
            for p=1:length(power_freqs)
                s  = 3/(2*pi*power_freqs(p));
                hi_freq_wavelet(p, :) = exp(2*pi*1i*power_freqs(p).*time) .* exp(-time.^2./(2*s^2));
                hi_freq_wavelet_fft(p, :)=fft((hi_freq_wavelet(p, :))',nconv);
                clear s
            end
            
            % FFT of Data
            data_fft_phase=zeros(size(data_for_phase, 1), nconv);
            for elec=1:size(data_for_phase,1)
                data_fft_phase(elec, :)=fft(data_for_phase(elec, :)', nconv);
            end
            
            % Phase information - conv data w/ lo freq
            A = zeros(1, (nconv));
            B = A(half_wavelet+1:nconv-half_wavelet);
            lf_phase=zeros(length(phas_freqs),size(data_for_phase, 2), size(data_for_phase, 1)); %phase freqXtimeXelec
            for elec = 1:size(data_for_phase,1) % loop thru elecs
                for f = 1:length(phas_freqs) % for each trial, loop thru freq
                    phase_conv = ifft(data_fft_phase(elec,:)' .* lo_freq_wavelet_fft(f,:)', nconv);
                    phase_conv = phase_conv(end:-1:1)';
                    phase_conv_result = phase_conv(half_wavelet:end-(half_wavelet+1));
                    lf_phase(f, :,elec) = angle(phase_conv_result);
                    clear phase_conv phase_conv_result
                end
            end
            
            % FFT of Data
            data_fft_power=zeros(size(data_for_power, 1), nconv );
            for elec=1:size(data_for_power, 1)
                data_fft_power(elec, :)=fft(data_for_power(elec, :)', nconv);
            end
            
            % Power information - conv data w/ hi freq
            hf_power=zeros(length(power_freqs),size(data_for_power, 2), size(data_for_power, 1));%hifreqXtimeXtrialXelec
            for elec = 1:size(data_for_power,1) % loop thru elecs
                for p = 1:length(power_freqs)
                    power_conv=ifft(data_fft_power(elec,:)' .* hi_freq_wavelet_fft(p,:)', nconv);
                    power_conv=power_conv(end:-1:1)';
                    power_conv_result=power_conv(half_wavelet:end-(half_wavelet+1));
                    hf_power(p, :, elec) = abs(power_conv_result.^2);
                    clear power_conv power_conv_result
                end
            end
            
            
            % Index power and phase for current time window
            hf_power_New=hf_power;
            lf_phase_New=lf_phase;
            
            lf_phase_Trans = zeros(size(lf_phase_New,2), size(lf_phase_New,1), size(lf_phase_New,3));
            
            % Transpose Phase Matrix for matrix mult
            for elec=1:size(lf_phase_New, 3) % loop thru elecs
                lf_phase_Trans(:, :,  elec)=lf_phase_New(:, :,  elec)';
            end
            
            % Multiply Phase matrix with i and get exponential: output=the second half of PAC eqn
            lf_phase_exp=exp(lf_phase_Trans*1i);
            
            % Get OBS PAC Values
            elec = 0
            OBS_PAC_New=zeros(size(hf_power_New, 1), size(lf_phase_New, 1), size(hf_power_New,4)*size(lf_phase_New,4)); %hifreqXlofreqXtrialXelec
            for FC_elec = 1:size(lf_phase_New, 3) % Loop thru FRO elec
                for HC_elec = 1:size(hf_power, 3) % Loop thru HC elec
                    elec = elec+1
                    OBS_PAC_New(:,:,elec)=abs((hf_power_New(:, :,HC_elec)*lf_phase_exp(:, :, FC_elec))/size(lf_phase_New, 2));
                end
            end
            
            
            % permutations
            tic
            elec_length      = size(hf_power_New,3)*size(lf_phase_New,3);
            num_iter         = 1000;
            random_timepoint = randsample(size(hf_power_New, 2), num_iter, 1);
            perm_hf_power    = zeros(size(hf_power_New));
            Perm_PAC         = zeros(size(hf_power_New, 1), size(lf_phase_exp, 2), num_iter,elec_length); %Power * Phase pac * iteration * elec
            
            tic
            for perm=1:num_iter
                %Shuffle power series in all trials elecs in first iteration
                perm_hf_power = hf_power_New(:, [random_timepoint(perm):end 1:real(random_timepoint(perm))-1], :,:); %powerfreqXtimesXtiralsXchans
                elec = 0;
                for FC_elec = 1:size(lf_phase_New, 3) % Loop thru FRO elec
                    for HC_elec = 1:size(hf_power, 3) % Loop thru HC elec
                        elec = elec+1;
                        Perm_PAC(:,:,perm,elec)=abs((perm_hf_power(:, :, HC_elec)*lf_phase_exp(:, :, FC_elec))/size(lf_phase_New, 2));
                    end
                end
                disp(perm)
            end
            toc
            
            mean_across_reps   = nanmean(Perm_PAC, 3);
            PAC_mn_Trial_Elec  = squeeze(mean_across_reps);
            std_across_reps    = nanstd(Perm_PAC,[],3);
            PAC_std_Trial_Elec = squeeze(std_across_reps);
            
            elec = 0;
            pacz  = zeros(size(hf_power_New, 1), size(lf_phase_exp, 2),elec_length);
            pacz_std_log_mtx = zeros(size(pacz));
            for FC_elec = 1:size(lf_phase_New, 4) % Loop thru FRO elec
                for HC_elec = 1:size(hf_power, 4) % Loop thru HC elec
                    elec = elec+1;
                    temp_mtx = zeros(size(PAC_std_Trial_Elec(:, :,elec)));
                    temp_mtx(find(PAC_std_Trial_Elec(:, :,elec)==0)) = 1;
                    pacz_std_log_mtx(:, :, elec)=temp_mtx;
                end
            end
            
            PAC_std_Trial_Elec(logical(pacz_std_log_mtx)) = nan;
            
            
            % Computing Z-scores
            elec = 0;
            for FC_elec = 1:size(lf_phase_New, 3) % Loop thru FRO elec
                for HC_elec = 1:size(hf_power, 3) % Loop thru HC elec
                    elec = elec+1
                    pacz(:, :, elec)=((OBS_PAC_New(:, :,elec))-(PAC_mn_Trial_Elec(:, :,elec)))./(PAC_std_Trial_Elec(:, :, elec)); % (obs pac - perm pac) / std of perm pac
                end
            end
            
            % overall PAC for each condition
            cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/figures/cfc'])
            h=figure; imagesc(squeeze(nanmean(pacz,3)));
            set(gca, 'xtick', 1:2:length(phas_freqs));
            x=floor(phas_freqs*10)/10;
            set(gca,'xticklabel', x(1:2:length(phas_freqs)));
            set(gca,'YDir','normal', 'ytick', 1:2:length(power_freqs)); %reverse axis and set ticks
            y=floor(power_freqs*10)/10;
            set(gca,'yticklabel', y(1:2:19));
            xlabel(['Freq for Phase in ' phase_region]);
            ylabel(['Freq for Power in ' power_region]);
            title(['cond' num2str(num2str(cond))])
            set(gca, 'FontSize', 14, 'FontWeight', 'bold')
            colorbar
            caxis([-.4 .7])
            save(['pacz_' phase_region '_' power_region '_cond' num2str(cond)], 'pacz', 'Perm_PAC','phas_freqs','power_freqs')
            saveas(gcf,[ 'allchans_' phase_region '_' power_region '_cond_' num2str(cond) '.bmp'])
            close all
            
        end
    end
end
toc