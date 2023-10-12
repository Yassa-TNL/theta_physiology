runs_num = 1000;
Gk = 0;
uk = 0;
gamma = eye(elec_length*order_counter);
kmax = length(order_counter+1:time_length);
P_1000 = zeros(size(Fk,1)); % covariance at diff time indices
P_2000 = zeros(size(Fk,1));
xhat0 = zeros(elec_length*order_counter, elec_length);
P0 = eye(elec_length*order_counter, elec_length*order_counter);
 for run_counter = 1:runs_num
    [xhist,zhist] = mcltisim(Fk,Gk,Hk,Qk,Rk,xhat0,P0,kmax); %
    % init conds
    xhatk = zeros(elec_length*order_counter, elec_length);
    Pk = eye(elec_length*order_counter, elec_length*order_counter);
    Fk = eye(elec_length*order_counter, elec_length*order_counter);
    Rk_temp = eye(elec_length);
    Qk = c*eye(elec_length*order_counter);

 end