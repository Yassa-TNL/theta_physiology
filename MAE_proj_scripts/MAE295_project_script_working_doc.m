% use GLKF to estimate TV-MVAR model
clear all;close all;clc

% develop for 1 subj
subj = '39';
cd(['/tmp/yassamri/iEEG/sandra/subj_' subj])
load(['trial_data_subj_onset_' subj '_ref__select_chan_3.mat'])
trial_data   = trial_data(:,:, 8:10);
%
% eamine data
% figure
%  for trl = elec3trl%size(trial_data,1)
%
%     plot(trial_data(trl,:,3))
%     hold on
%    % pause
%  end
%
%  figure
% plot(squeeze(nanmean(trial_data(elec2trl,:,2),1)))
% end
% 2 3 4 5 6  7  8  9 10 11 12 13 14 15
elec1trl = [ 5 6 7 8 9 10 11  13 14 15 16 17 19 20]; % elec 1
elec2trl = [5 7 8 13 14 16 17 18 20 21 23 24 28 29]; % elec 2
elec3trl = [5 7 8 13 14 16 17 18 20 21 23 24 28 29]; % elec 2


% made data stationary (mean, cov fcn of lag), remove trials with artifacts
sig_avg = squeeze(mean(trial_data, 1));
for elec = 1:size(trial_data, 3)
    for trials = 1:size(trial_data,1) % remove mean from each trial
        trial_data(trials, :,elec) = trial_data(trials, :,elec) - sig_avg(:,elec)';
    end
end

trial_data_new(:,:,1) = trial_data(elec1trl,:,1);
trial_data_new(:,:,2) = trial_data(elec2trl,:,2);
trial_data_new(:,:,3) = trial_data(elec3trl,:,3);
clear trial_data
trial_data = trial_data_new;

%% dimensions
trial_length = size(trial_data,1);
elec_length  = size(trial_data,3);
time_length  = size(trial_data,2);

Xn = permute(trial_data, [1 3 2]); % trials X elec X time: observation mtx
c  = 0.03; % constant for mst and process noise estimates
orders_tested= 8;

% matrices to save:
P_order    = cell(1,orders_tested);
xhat_order = cell(1,orders_tested);
bic_order  = cell(1,orders_tested);
%Ev         = cell(1,orders_tested);
%% Ev = zeros(orders_tested,1);
alpha = 0.05;
test_pass = [];
det_hist = [];
for order_counter =1:15 % 1:orders_tested % run the model for diff orders
    % init conds
    xhatk = zeros(elec_length*order_counter, elec_length);
    Pk = eye(elec_length*order_counter, elec_length*order_counter);
    Fk = eye(elec_length*order_counter, elec_length*order_counter);
    Rk_temp = eye(elec_length);
    Qk = c*eye(elec_length*order_counter);
    
    % initalize mtx
    P    = nan(elec_length*order_counter, elec_length*order_counter, length(order_counter+1:time_length)); %dpxdpxtimelength once per order
    xhat = zeros(elec_length*order_counter, elec_length, length(order_counter+1:time_length));
    bic  = zeros(length(order_counter+1:time_length),1);
    % run filter for diff orders
    counter = 0;
    Ev_temp = 0;
    
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
        Rk = trace(Rk_temp*(1-c)+((c*(Zn-Hk*xhatk)')*(Zn-Hk*xhatk))/trial_length-1)*eye(size(Hk,1));
        
        % kalman filter
        % eqn 1: state prediction covariance
        Pbarkplus1= Fk*Pk*Fk' + Qk; %+ Qk; % calc Qk
        
        % eqn 2: mst prection covariance
        Skplus1 = Hk*Pbarkplus1*Hk'+Rk;
        
        % eqn3: filter gain
        Wkplus1 = Pbarkplus1*Hk'/Skplus1;
        
        % eqn 4: updated covariance
        %Pkplus1 = Pbarkplus1 - Wkplus1*Skplus1*Wkplus1';
        
        %use joseph formula
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
        
        % compute Bayes Information Criteria
        bic(counter,1) = log(det(Pk)) + (log(length(order_counter+1:time_length))*order_counter*2^2)/length(order_counter+1:time_length);
        %figure; imagesc( P_order{1}); colorbar
        Ev_temp = Ev_temp + trace(nukplus1'/Skplus1*nukplus1);
        % calculate Ev
        Ev_temp(:,:,counter) = nukplus1'/Skplus1*nukplus1;
        
        
    end
    
    % save determinant and Ev
    det_hist = [det_hist; det(Pkplus1)];
    %Ev{order_counter} = Ev_temp;
    Ev(order_counter) = Ev_temp/time;

    test_dist = Ev_temp;
    r1 = chi2inv(alpha/2, (time_length - order_counter)*numel(Zn));
    r2 = chi2inv(1 - alpha/2, (time_length - order_counter)*numel(Zn));
    test_pass = [test_pass, (r1<test_dist) && (test_dist<r2)];
    % save estimate and covariance as a function of order num
    P_order{order_counter}    = P;
    xhat_order{order_counter} = xhat;
    bic_order{order_counter}  = bic;
end


%% plot to check filter
figure; plot(Ev)
figure; plot(det_hist)
find(det_hist < 0)
%% measure filter consistencey
% run montel carlo for each order
order_counter = 8;
runs_num = 1000;
Gk = 0;
uk = 0;
gamma = eye(elec_length*order_counter);
kmax = length(order_counter+1:time_length);

P_1000 = zeros(size(Fk,1)); % covariance at diff time indices
P_2000 = zeros(size(Fk,1));
xhat0 = zeros(elec_length*order_counter, elec_length);
xtilde_avg = zeros([size(xhat0),kmax]);
P0 = eye(elec_length*order_counter, elec_length*order_counter);
for run_counter = 1:runs_num
    [xhist,zhist] = mcltisim(Fk,Gk,Hk,Qk,Rk,xhat0,P0,kmax); %
    % init conds
    xhatk = zeros(elec_length*order_counter, elec_length);
    Pk = eye(elec_length*order_counter, elec_length*order_counter);
    Fk = eye(elec_length*order_counter, elec_length*order_counter);
    Rk_temp = eye(elec_length);
    Qk = c*eye(elec_length*order_counter);
    
    % initalize mtx
    P    = nan(elec_length*order_counter, elec_length*order_counter, length(order_counter+1:time_length)); %dpxdpxtimelength once per order
    xhat = zeros(elec_length*order_counter, elec_length, length(order_counter+1:time_length));
    bic  = zeros(length(order_counter+1:time_length),1);
    % run filter for diff orders
    counter = 0;
    %Ev_temp = 0;
    
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
        Rk = trace(Rk_temp*(1-c)+((c*(Zn-Hk*xhatk)')*(Zn-Hk*xhatk))/trial_length-1)*eye(size(Hk,1));
        
        % kalman filter
        % eqn 1: state prediction covariance
        Pbarkplus1= Fk*Pk*Fk' + Qk; %+ Qk; % calc Qk
        
        % eqn 2: mst prection covariance
        Skplus1 = Hk*Pbarkplus1*Hk'+Rk;
        
        % eqn3: filter gain
        Wkplus1 = Pbarkplus1*Hk'/Skplus1;
        
        % eqn 4: updated covariance
        %Pkplus1 = Pbarkplus1 - Wkplus1*Skplus1*Wkplus1';
        
        %use joseph formula
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
        
        % compute Bayes Information Criteria
        bic(counter,1) = log(det(Pk)) + (log(length(order_counter+1:time_length))*order_counter*2^2)/length(order_counter+1:time_length);
        %figure; imagesc( P_order{1}); colorbar
        Ev_temp = Ev_temp + trace(nukplus1'/Skplus1*nukplus1);
        % calculate Ev
        Ev_temp(:,:,counter) = nukplus1'/Skplus1*nukplus1;
        
    end
    xtilde = xhist(:,:,2:end) - xhat;
    xtilde_avg = xtilde_avg + xtilde;
    P_1000 = P_1000 + xtilde(:,1000)*xtilde(:,1000)';
    P_2000 = P_2000 + xtilde(:,2000)*xtilde(:,2000)';
end

xtilde_avg = xtilde_avg/runs_num;
P_1000 = P_1000/runs_num;% average covariance at time k
P_2000 = P_2000/runs_num;% average covariance at time k

save('projMAE295', 'P_1000','P_2000','xtilde_avg')

%% check AR model

Aks = xhat_order{8};

% recover AKls
ii = 1
A1 = Aks((ii-1)*elec_length+(1:elec_length),1:elec_length,:);
ii = 2
A2 = Aks((ii-1)*elec_length+(1:elec_length),1:elec_length,:);
ii = 3
A3 = Aks((ii-1)*elec_length+(1:elec_length),1:elec_length,:);
ii = 4
A4 = Aks((ii-1)*elec_length+(1:elec_length),1:elec_length,:);
ii = 5
A5 = Aks((ii-1)*elec_length+(1:elec_length),1:elec_length,:);
ii = 6
A6 = Aks((ii-1)*elec_length+(1:elec_length),1:elec_length,:);
ii = 7
A7 = Aks((ii-1)*elec_length+(1:elec_length),1:elec_length,:);
ii = 8
A8 = Aks((ii-1)*elec_length+(1:elec_length),1:elec_length,:);

%% predic 9th point till end
trial = 4;
yhat = zeros(elec_length, length((order_counter+1):size(trial_data,2)));
for time = 1:size(A1,3)
yhat(:,time)= A1(:,:,time)'*squeeze(trial_data(trial, order_counter+time, :))+ A2(:,:,time)'*squeeze(trial_data(trial, order_counter+time, :))+...
    A3(:,:,time)'*squeeze(trial_data(trial, order_counter+time, :))+A4(:,:,time)'*squeeze(trial_data(trial, order_counter+time, :))+...
    A5(:,:,time)'*squeeze(trial_data(trial, order_counter+time, :))+A6(:,:,time)'*squeeze(trial_data(trial, order_counter+time, :))+...
   A7(:,:,time)'*squeeze(trial_data(trial, order_counter+time, :))+A8(:,:,time)'*squeeze(trial_data(trial, order_counter+time, :));
end

figure
elec = 1
plot(yhat(elec,:))
hold on
plot(trial_data(trial,(order_counter+1):end,elec))
legend('estimated', 'observed')



