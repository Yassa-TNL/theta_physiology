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
c  =1; % constant for mst and process noise estimates
order_counter = 10;
sim_legnth = 1000;
% init conds
xbar0 = zeros(elec_length*order_counter, elec_length);
Pk = eye(elec_length*order_counter, elec_length*order_counter);
Fk = eye(elec_length*order_counter, elec_length*order_counter);
Rk_temp = eye(elec_length);
Qk = c*eye(elec_length*order_counter);
Ev_temp = zeros(elec_length,time_length);
Ev=zeros(elec_length,time_length);
for sim_num = 1:sim_legnth
    
    % initalize mtx
    P    = nan(elec_length*order_counter, elec_length*order_counter, length(order_counter+1:time_length)); %dpxdpxtimelength once per order
    
    % draw and xhat
    Pk = eye(elec_length*order_counter, elec_length*order_counter);
    xhatk = xbar0 + chol(Pk)'*randn(size(xbar0));
    
    % run filter for diff orders
    counter = 0;
    for time = order_counter+1:time_length
        
        %time
        counter = counter+1;
        
        % get an observation
        Zn = Xn(:,:,time);
        
        % get measurement mtx
        Hk = [];
        for a = 1:order_counter % loop thru order
            Hk =[Hk Xn(:,:,time-a)]; % 140 trials X (elec*order)
        end
        
        % Rk
        Rk_temp = Rk_temp*(1-c)+((c*(Zn-Hk*xhatk)')*(Zn-Hk*xhatk))/trial_length-1;
        Rk = trace(Rk_temp)*eye(size(Hk,1));
        
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
        
        % calculate Ev
        Ev_temp(:,counter) = diag(nukplus1'/Skplus1*nukplus1);
    end
    
    % save Ev
    Ev= Ev+Ev_temp;
end
toc

figure;
plot(mean(Ev,1)/sim_legnth)
% calc r1 and r2
alpha = 0.05;
r1 = chi2inv(alpha/2, sim_legnth)/sim_legnth
r2 = chi2inv(1-alpha/2, sim_legnth)/sim_legnth