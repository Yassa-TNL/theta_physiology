% use General Linear Kalman Filter to estimate TV-MVAR model
clear all; close all; clc
tic

% develop for 1 subj
subj = '39';
cd(['/tmp/yassamri/iEEG/sandra/subj_' subj])
load(['trial_data_subj_onset_' subj '_ref__select_chan_3.mat'])
trial_data = trial_data(1:30,1:1500,10:19); % test on 30 trials, 1.5 seconds, 10 elecs

% made data stationary (mean, cov fcn of lag), remove trials with artifacts
sig_avg = squeeze(mean(trial_data, 1));
for elec = 1:size(trial_data, 3)
    for trials = 1:size(trial_data,1) % remove mean from each trial
        trial_data(trials, :,elec) = trial_data(trials, :,elec) - sig_avg(:,elec)';
    end
end

% dimensions
trial_length = size(trial_data,1);
elec_length  = size(trial_data,3);
time_length  = size(trial_data,2);

Xn = permute(trial_data, [1 3 2]); % trials X elec X time: observation mtx
c_vals  =1; % constant for mst and process noise estimates
order_start = 10;
order_end   = 10;
orders_tested = length(order_start:order_end);

% matrices to save:
P_order      = cell(1,orders_tested);
xhat_order   = cell(1,orders_tested);
P_c_order    = cell(1,length(c_vals));
xhat_c_order = cell(1,length(c_vals));
Ev           = zeros(elec_length,time_length);
det_hist     = [];


for c_counter = 1:length(c_vals)
    c = c_vals(c_counter)
    order_counter_ii = 0;
    for order_counter =order_start:order_end % 1:orders_tested % run the model for diff orders
        order_counter
        order_counter_ii = order_counter_ii+1;
        % init conds
        xhatk = zeros(elec_length*order_counter, elec_length);
        Pk = eye(elec_length*order_counter, elec_length*order_counter);
        Fk = eye(elec_length*order_counter, elec_length*order_counter);
        Rk_temp = eye(elec_length);
        % Rk_temp = eye(trial_length);
        Qk = c*eye(elec_length*order_counter);
        
        % initalize mtx
        P    = nan(elec_length*order_counter, elec_length*order_counter, length(order_counter+1:time_length)); %dpxdpxtimelength once per order
        xhat = zeros(elec_length*order_counter, elec_length, length(order_counter+1:time_length));
        bic  = zeros(length(order_counter+1:time_length),1);
        
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
            
            % Rk
            % approach 2
            Rk_temp = Rk_temp*(1-c)+((c*(Zn-Hk*xhatk)')*(Zn-Hk*xhatk))/trial_length-1;
            Rk = trace(Rk_temp)*eye(size(Hk,1));
            
            % approach 1
            %Rk = trace(Rk_temp*(1-c)+((c*(Zn-Hk*xhatk)')*(Zn-Hk*xhatk))/trial_length-1)*eye(size(Hk,1));
            
            % approach 3
            %  Rk_temp = Rk_temp*(1-c)+((c*(Zn-Hk*xhatk))*(Zn-Hk*xhatk)')/trial_length-1;
            %  Rk=Rk_temp;
            % kalman filter
            % eqn 1: state prediction covariance
            Pbarkplus1= Fk*Pk*Fk' + Qk; %+ Qk; % calc Qk
            
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
            
            
            % for next iteration
            xhatk = xhatkplus1; %next prediction is previous estimate
            Pk = Pkplus1; %next state prediction error is previous state prediction error
            
            % compute Bayes Information Criteria
            bic(counter,1) = log(det(Pk)) + (log(length(order_counter+1:time_length))*order_counter*2^2)/length(order_counter+1:time_length);
            
            % calculate Ev
            Ev_temp(:,counter) = diag(nukplus1'/Skplus1*nukplus1);

        end
        
        % save determinant and Ev
        det_hist = [det_hist; det(Pkplus1)];
%        Ev{order_counter_ii} = Ev_temp/time;
        
        % save estimate and covariance as a function of order num
        P_order{order_counter_ii}    = P;
        xhat_order{order_counter_ii} = xhat;
    end
    
    % save estimate and covariance as a function of order num
    P_c_order    {c_counter} = P_order;
    xhat_c_order {c_counter} = xhat_order;
    
end
toc
%% order vs. RMS plot
tic
order_c_yhat    = cell(1,length(c_vals));

for c_counter = 1:length(c_vals)
    order_yhat    = cell(1,orders_tested);
    order_counter_ii = 0;
    xhat_order = xhat_c_order {c_counter} ;
    for order =order_start:order_end
        order_counter_ii = order_counter_ii+1;
        Aks = xhat_order{order_counter_ii};
        yhat = zeros(elec_length, length((order+1):size(trial_data,2)), trial_length);
        for trial = 1:trial_length
            for time = 1:size(xhat_order{order_counter_ii},3) % for a given timepoint
                predic_temp=zeros(size(xhat_order{order_counter_ii},2),1);
                for ii = 1:order       % sum current point and all past (loop thru order)
                    predic_temp = predic_temp+Aks((ii-1)*elec_length+(1:elec_length),1:elec_length,time)'*squeeze(trial_data(trial, (order+time)-ii, :));
                end
                yhat(:,time, trial)= predic_temp;
            end
        end
        order_yhat{order_counter_ii}= yhat; % save predicted signal for each order
        clear AKs
    end
    order_c_yhat{c_counter} = order_yhat;
end
toc

% plot estiamted vs. real
for c_counter = 1:length(c_vals)
    get_c = order_c_yhat{c_counter} ;
    c = c_vals(c_counter);

    order_counter_ii = 0;
    for order =order_start:order_end
        order_counter_ii = order_counter_ii+1;
        get_c_order_yhat = get_c{order_counter_ii};
        figure
        elec = 1
        plot(get_c_order_yhat(elec,:, trial), 'LineWidth', 2, 'MarkerFaceColor', 'm')%+50
        hold on
        plot(trial_data(trial,(order+1):end,elec), 'LineWidth', 2, 'MarkerFaceColor', 'k')
        legend('estimated', 'observed')
        xlabel('time (ms)')
        ylabel('voltage')
        title(['cval = ' num2str(c) ' order = ' num2str(order)])
        set(gca, 'FontSize', 16, 'FontWeight', 'bold')
        saveas(gcf, ['real_predict_signals_' 'cval_' num2str(c) '_order_' num2str(order) '.bmp'])
    end
end
        xlabel('microvolts')
        ylabel('count')
saveas(gcf, 'RMS_signals_order10c1.bmp')
saveas(gcf, 'hist_count_of_data.bmp')

tic
% save mean RMS
order_rms_mn  = zeros(orders_tested, length(c_vals));
order_rms_std = zeros(orders_tested, length(c_vals));
for c_counter = 1:length(c_vals)
    get_c = order_c_yhat{c_counter} ;
     order_counter_ii = 0;
    for  order = order_start:order_end
        order_counter_ii = order_counter_ii+1;
        yhat = get_c{order_counter_ii};
        elec_trial_RMS=[];
        for elec = 1:elec_length
            for trial = 1:trial_length
                elec_trial_RMS= [elec_trial_RMS sqrt(mean((yhat(elec,:, trial)-trial_data(trial,(order+1):end,elec)).^2))];
            end
        end
        order_rms_mn  (order_counter_ii, c_counter) = mean(elec_trial_RMS);
        order_rms_std (order_counter_ii, c_counter) = std(elec_trial_RMS,0,2);
        clear yhat
    end
end
%cd /tmp/yassamri/iEEG/sandra/analysis_pipeline_final/MAE_proj_scripts/MAE_proj_figs
%save(['subj_' subj '_RMS_workspace_orders_' num2str(order_start) '_' num2str(order_end)], 'order_yhat', 'order_rms_mn', 'order_rms_std')
toc

%plot order vs. RMS

cd /tmp/yassamri/iEEG/sandra/analysis_pipeline_final/MAE_proj_scripts/MAE_proj_data_figs
figure
hold on
for c_counter = 1:length(c_vals)
    plot(order_rms_mn(:,c_counter), 'LineWidth', 2, 'MarkerFaceColor', 'm')
    xlabel('model order')
    ylabel('RMSE')
    set(gca, 'FontSize', 16, 'FontWeight', 'bold')
end
saveas(gcf, 'C_Order_RMS_signals.bmp')

figure
plot(c_vals(4:end),order_rms_mn(10,4:end), 'LineWidth', 3, 'MarkerFaceColor', 'm')
xlabel('filter paramater c')
ylabel('RMSE')
title('model order 10')
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
saveas(gcf, 'C_Order10_RMS_signals_all_3.bmp')


