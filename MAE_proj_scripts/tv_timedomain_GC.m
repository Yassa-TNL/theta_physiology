clear all;close all;clc
tic
% develop for 1 subj
cond_num  = 2
reg1_name = 'OFC';
reg2_name = 'HC';

subj = '84';
exp_type ='tuning_correct';
lock     = 'onset';

cd(['/tmp/yassamri/iEEG/sandra/subj_' subj])
addpath('/tmp/yassamri/iEEG/sandra/analysis_pipeline_final')
[cond1,cond2,cond3,cond4, cond5, cond6] = GetCondData(subj, exp_type, lock);
[fro_chan_idx,MTL_chan_idx,temp_chan_idx,insula_chan_idx,cingulate_chan_idx,...
    OFC_chan_idx,CA3_chan_idx,CA1_chan_idx,HC_chan_idx] = get_elecs(subj);
if strcmp('OFC',reg1_name)
    reg1 = OFC_chan_idx;
elseif strcmp('FRO',reg1_name)
    reg1 = fro_chan_idx;
elseif strcmp('TEMP',reg1_name)
    reg1 = temp_chan_idx;
elseif strcmp('CING',reg1_name)
    reg1 = cingulate_chan_idx;
elseif strcmp('INS',reg1_name)
    reg1 = insula_chan_idx;
end

if strcmp('HC',reg2_name)
    reg2 = HC_chan_idx;
elseif strcmp('CA3',reg2_name)
    reg2 = CA3_chan_idx;
elseif strcmp('CA1',reg2_name)
    reg2 = CA1_chan_idx;
end

allreg = [reg1 reg2];

if cond_num==2
    cond  = 'PatternCompletion';
    cond2 = cond2(:,:,allreg);
    trial_data = cond2;
    disp('cond2-pattern comp')
elseif cond_num==3
    cond  = 'PatternSeperation';
    cond3 = cond3(:,:,allreg);
    trial_data = cond3;
    disp('cond3-pattern sep')
end


% made data stationary (mean, cov fcn of lag), remove trials with artifacts
sig_avg = squeeze(mean(trial_data, 1));
for elec = 1:size(trial_data, 3)
    for trials = 1:size(trial_data,1) % remove mean from each trial
        trial_data(trials, :,elec) = trial_data(trials, :,elec) - sig_avg(:,elec)';
    end
end

Xn = permute(trial_data, [1 3 2]); % trials X elec X time: observation mtx
c  = 1; % constant for mst and process noise estimates
order_start = 10;
order_end   = 10;
orders_tested = length(order_start:order_end);

% matrices to save:
P_order    = cell(1,orders_tested);
xhat_order = cell(1,orders_tested);
Ev         = cell(1,orders_tested);
order_counter_ii = 0;


% first get PE(n)
% dimensions
trial_length = size(trial_data,1);
elec_length  = size(trial_data,3);
time_length  = size(trial_data,2);

%init
wnbar   = nan(size(trial_data,3), size(trial_data,3), time_length-order_start);

order_counter =10; % 1:orders_tested % run the model for diff orders
order_counter_ii = order_counter_ii+1;

% init conds
xhatk = zeros(elec_length*order_counter, elec_length);
Pk = eye(elec_length*order_counter, elec_length*order_counter);
Fk = eye(elec_length*order_counter, elec_length*order_counter);
Rk_temp = eye(elec_length);
Qk = c*eye(elec_length*order_counter);

% initalize mtx
P    = nan(elec_length*order_counter, elec_length*order_counter, length(order_counter+1:time_length)); %dpxdpxtimelength once per order
xhat = zeros(elec_length*order_counter, elec_length, length(order_counter+1:time_length));

% run filter for diff orders
counter = 0;
for time = order_counter+1:time_length
    
    %time = time+1
    counter = counter+1;
    
    % get an observation
    Zn = Xn(:,:,time);
    
    % get measurement mtx
    Hk = [];
    for a = 1:order_counter % loop thru order
        Hk =[Hk Xn(:,:,time-a)]; % 140 trials X (elec*order)
    end
    
    Rk_temp = Rk_temp*(1-c)+((c*(Zn-Hk*xhatk)')*(Zn-Hk*xhatk))/trial_length-1;
    Rk = trace(Rk_temp)*eye(size(Hk,1));

    wnbar(:,:,time) = Rk_temp;
     
    % kalman filter
    % eqn 1: state prediction covariance
    Pbarkplus1= Fk*Pk*Fk' + Qk;
    
    % eqn 2: mst prection covariance
    Skplus1 = Hk*Pbarkplus1*Hk'+Rk;
    
    % eqn3: filter gain
    Wkplus1 = Pbarkplus1*Hk'/Skplus1;
    
    % eqn 4: updated covariance
    % use joseph form
    Pkplus1 = (eye(size(Hk,2))- Wkplus1*Hk)*Pbarkplus1*(eye(size(Hk,2))- Wkplus1*Hk)' + Wkplus1*Rk*Wkplus1';
    P(:,:,counter) = Pkplus1; %save estimate covariance over time
    
    % eqn 5: predicted state
    xbarkplus1 = Fk*xhatk; %predicted state
    
    % eqn 6: predicted measurement
    zbarkplus1 = Hk*xbarkplus1; %predicted measurement
    zkplus1 =Xn(:,:,time); %get measurement
    
    % eqn 7: innovation term / measurement residual
    nukplus1 = zkplus1 - zbarkplus1; %calculate residual
    
    % eqn 8: updated state estimate
    xhatkplus1 = xbarkplus1 + Wkplus1*nukplus1;%state estiamte (filtered prediction)
    xhat(:,:,counter) = xhatkplus1;
    
    % eqn 9: state predic error
    xtildekplus1 = xhatkplus1 - xbarkplus1;%state prediction error
    
    % for next iteration
    xhatk = xhatkplus1; %next prediction is previous estimate
    Pk = Pkplus1; %next state prediction error is previous state prediction error  
end



% take one elec off at a time
tvpGCI       = nan(size(trial_data,3), size(trial_data,3), time_length-order_start);
for elec = 1:size(trial_data,3)
    disp(elec)
    trial_data_new = trial_data(:,:,setdiff(1:size(trial_data,3),elec));
    Xn_new = Xn(:,setdiff(1:size(trial_data,3),elec),:);
    
    % dimensions
    elec_length  = size(trial_data_new,3);
    
    order_counter =10; % 1:orders_tested % run the model for diff orders
    order_counter_ii = order_counter_ii+1;
    
    % init conds
    xhatk = zeros(elec_length*order_counter, elec_length);
    Pk = eye(elec_length*order_counter, elec_length*order_counter);
    Fk = eye(elec_length*order_counter, elec_length*order_counter);
    Rk_temp = eye(elec_length);
    Qk = c*eye(elec_length*order_counter);
    
    % initalize mtx
    P    = nan(elec_length*order_counter, elec_length*order_counter, length(order_counter+1:time_length)); %dpxdpxtimelength once per order
    xhat = zeros(elec_length*order_counter, elec_length, length(order_counter+1:time_length));
    
    counter = 0;
    for time = order_counter+1:time_length
        
        %time = time+1
        counter = counter+1;
        % get an observation
        Zn = Xn_new(:,:,time);
        
        % get measurement mtx
        Hk = [];
        for a = 1:order_counter % loop thru order
            Hk =[Hk Xn_new(:,:,time-a)]; % 140 trials X (elec*order)
        end
        Rk_temp = Rk_temp*(1-c)+((c*(Zn-Hk*xhatk)')*(Zn-Hk*xhatk))/trial_length-1;
        Rk = trace(Rk_temp)*eye(size(Hk,1));
        
        wnbar_red = Rk_temp;
        
        tvpGCI(elec,setdiff(1:size(trial_data,3),elec),time-order_start) = log(diag(wnbar_red)./diag(wnbar(setdiff(1:size(trial_data,3),elec),setdiff(1:size(trial_data,3),elec),time)));
        
        % kalman filter
        % eqn 1: state prediction covariance
        Pbarkplus1= Fk*Pk*Fk' + Qk;
        
        % eqn 2: mst prection covariance
        Skplus1 = Hk*Pbarkplus1*Hk'+Rk;
        
        % eqn3: filter gain
        Wkplus1 = Pbarkplus1*Hk'/Skplus1;
        
        % eqn 4: updated covariance
        % use joseph form
        Pkplus1 = (eye(size(Hk,2))- Wkplus1*Hk)*Pbarkplus1*(eye(size(Hk,2))- Wkplus1*Hk)' + Wkplus1*Rk*Wkplus1';
        P(:,:,counter) = Pkplus1; %save estimate covariance over time
        
        % eqn 5: predicted state
        xbarkplus1 = Fk*xhatk; %predicted state
        
        % eqn 6: predicted measurement
        zbarkplus1 = Hk*xbarkplus1; %predicted measurement
        zkplus1 =Xn_new(:,:,time); %get measurement
        
        % eqn 7: innovation term / measurement residual
        nukplus1 = zkplus1 - zbarkplus1; %calculate residual
        
        % eqn 8: updated state estimate
        xhatkplus1 = xbarkplus1 + Wkplus1*nukplus1;%state estiamte (filtered prediction)
        xhat(:,:,counter) = xhatkplus1;
        % eqn 9: state predic error
        xtildekplus1 = xhatkplus1 - xbarkplus1;%state prediction error
        
        % for next iteration
        xhatk = xhatkplus1; %next prediction is previous estimate
        Pk = Pkplus1; %next state prediction error is previous state prediction error
        
        % save estimate and covariance as a function of order num
        P_order{order_counter_ii}    = P;
        xhat_order{order_counter_ii} = xhat;
        
    end
    
end
toc
cd /tmp/yassamri/iEEG/sandra/analysis_pipeline_final/MAE_proj_scripts/MAE_proj_data_figs
save(['subj' subj '_' cond '_' reg1_name 'vs' reg2_name],'tvpGCI' )


%%
mx = .1
mn = 0
figure
subplot(1,3,1)
imagesc(nanmean(tvpGCI(:,:,1:500), 3))
caxis([mn mx])
title('prestim')
colorbar
set(gca, 'Fontsize', 12, 'FontWeight' , 'Bold')

subplot(1,3,2)
imagesc(nanmean(tvpGCI(:,:,500:1500), 3))
title('1 sec post-stim')
caxis([mn mx])
colorbar
set(gca, 'Fontsize', 12, 'FontWeight' , 'Bold')

subplot(1,3,3)
imagesc(nanmean(tvpGCI(:,:,1500:end), 3))
title('2 sec post-stim')
caxis([mn mx])
colorbar
suptitle( [cond ' ' reg1_name ' vs. ' reg2_name])

set(gca, 'Fontsize', 12, 'FontWeight' , 'Bold')
saveas(gcf, [cond '_' reg1_name 'vs' reg2_name '.bmp'])



%% loop thru cortex elec
NC_HC = zeros(length(reg1)*length(reg2), size(tvpGCI,3))
cntr = 0
for NC_elec = 1:length(reg1)
    for HC_elec = length(reg1)+1:size(tvpGCI,2)
        cntr = cntr+1
       NC_HC(cntr,:) =tvpGCI(NC_elec, HC_elec,:);
    end
end

HC_NC = zeros(length(reg1)*length(reg2), size(tvpGCI,3))
cntr = 0
for HC_elec = length(reg1)+1:size(tvpGCI,2)
    for NC_elec =  1:length(reg1)
        cntr = cntr+1
       HC_NC(cntr,:) =tvpGCI(HC_elec, NC_elec,:);
    end
end

% plot time traces
figure
plot(nanmean(HC_NC(1:1500),1))
hold on
plot(nanmean(NC_HC(1:1500),1))
legend({'HC-NC', 'NC-HC'})

% plot imagesc for NC-HC
row_mean = [nanmean(tvpGCI(1:length(reg1),:,1:500-order_counter),1) ; nanmean(tvpGCI((length(reg1)+1):end,:,1:500-order_counter),1) ];
all_mean = [nanmean(row_mean(:,1:length(reg1),:),2) nanmean(row_mean(:,(length(reg1)+1):end,:),2) ];% OFC vs. HC

figure
subplot(3,1,1)
imagesc(nanmean(all_mean(:,:,:),3)) 
set(gca, 'XTick', 1:2, 'XTickLabel', {'NC', 'HC'},'YTick', 1:2,'YTickLabel',{'NC', 'HC'},'Fontsize',12,'FontWeight' ,'Bold')
suptitle( [cond ' ' reg1_name ' vs. ' reg2_name '0.5 sec prestim'])
colorbar

row_mean = [nanmean(tvpGCI(1:length(reg1),:,500-order_counter+1:500-order_counter+1001),1) ; nanmean(tvpGCI((length(reg1)+1):end,:,500-order_counter+1:500-order_counter+1001),1) ];
all_mean = [nanmean(row_mean(:,1:length(reg1),:),2) nanmean(row_mean(:,(length(reg1)+1):end,:),2) ];% OFC vs. HC

subplot(3,1,2)
imagesc(nanmean(all_mean(:,:,:),3)) 
set(gca, 'XTick', 1:2, 'XTickLabel', {'NC', 'HC'},'YTick', 1:2,'YTickLabel',{'NC', 'HC'},'Fontsize',12,'FontWeight' ,'Bold')
suptitle( [cond ' ' reg1_name ' vs. ' reg2_name '0.5 sec prestim'])
colorbar

row_mean = [nanmean(tvpGCI(1:length(reg1),:,500-order_counter+1:500-order_counter+1001),1) ; nanmean(tvpGCI((length(reg1)+1):end,:,500-order_counter+1:500-order_counter+1001),1) ];
all_mean = [nanmean(row_mean(:,1:length(reg1),:),2) nanmean(row_mean(:,(length(reg1)+1):end,:),2) ];% OFC vs. HC

subplot(3,1,2)
imagesc(nanmean(all_mean(:,:,:),3)) 
set(gca, 'XTick', 1:2, 'XTickLabel', {'NC', 'HC'},'YTick', 1:2,'YTickLabel',{'NC', 'HC'},'Fontsize',12,'FontWeight' ,'Bold')
suptitle( [cond ' ' reg1_name ' vs. ' reg2_name '0.5 sec prestim'])
colorbar
saveas(gcf, [cond '_' reg1_name 'vs' reg2_name 'collapsed.bmp'])


%%

find(tvpGCI(imag(tvpGCI) ~= 0)) 
tvpGCI(imag(tvpGCI) ~= 0) = NaN;

%tvpGCI(imag(tvpGCI) ~= 0) = NaN;


% figure
% imagesc(nanmean(all_mean(:,:,1:250), 3))
% caxis([0 0.036])
% title('prestim')
% colorbar
% 
% figure
% imagesc(nanmean(all_mean(:,:,250:1250), 3))
% title('first second')
% caxis([0 0.06])
% colorbar
% 
% figure
% imagesc(nanmean(all_mean(:,:,1250:2491), 3))
% title('last second')
% caxis([0 0.07])
% colorbar

% 
% 
% %%
%cd /tmp/yassamri/iEEG/sandra/analysis_pipeline_final/MAE_proj_scripts/MAE_proj_figs
%
%% 
% figure
% for time_cntr = 1:size(tvpGCI,3)
%  
%    pause(.1)
% end

