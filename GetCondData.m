function [cond1,cond2,cond3,cond4, cond5, cond6] = GetCondData(subj,exp_type, lock, DS, fn_nm, ref)
% output: trial data in time domain for each condition

% go to ptnt dir
cd(['/mnt//yassamri/iEEG/sandra/subj_' subj])

% load trial data
if strcmp('yes',DS)
    load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3' fn_nm '_fs_500.mat'])
elseif strcmp('',DS)
    load(['trial_data_subj_' lock '_' subj '_ref_' ref '_select_chan_3_NotDownsampled.mat'])
end

% remove trils w/ no resp
if strcmp('response', lock)   
    bd_trls = zeros(1,length(responses));
    bd_trls(responses<.2) = 1;
    trial_data(logical(bd_trls),:,:) = nan;
end

% load behavior
load(['behavior_subj' subj '.mat'])

% organize data into condition
study_trial_data = trial_data(1:nTrials_study,:,:);
test_trial_data  = trial_data(nTrials_study+1:end,:,:);


repeat_log_vec       = testing_behav_matrix(:,2)== 0;
new_log_vec          = testing_behav_matrix(:,2)== 0.5;
all_lure_log_vec     = testing_behav_matrix(:,2)>=1;
corr_resp_log_vec    = testing_behav_matrix(:,3)== testing_behav_matrix(:,5);


%diff from repeats --> similar to repeats 5-->1
lure_1_log_vec     =testing_behav_matrix(:,2)== 1;
lure_2_log_vec     =testing_behav_matrix(:,2)== 2;
lure_4_log_vec     =testing_behav_matrix(:,2)== 4;
lure_5_log_vec     =testing_behav_matrix(:,2)== 5;
lures              =[1 2 4 5];
diff               =[1 2 ];
lure_diff          =logical(ismember(testing_behav_matrix(:,2),diff));
lure_all           =logical(ismember(testing_behav_matrix(:,2),lures));

cond1= [];
cond2= [];
cond3= [];
cond4= [];
cond5= [];
cond6= [];

% indoor/outdoor
indoor_log_vec  = training_behav_matrix(:,3)== 1;
outdoor_log_vec = training_behav_matrix(:,3)== 2;

if strcmp('study_test',exp_type)
    cond1 = study_trial_data;
    cond2 = test_trial_data;
    cond3= [];
    cond4= [];
    cond5= [];
    cond6= [];
elseif strcmp('tuning',exp_type)
    cond1= test_trial_data(repeat_log_vec,:,:);
    cond2= test_trial_data(lure_1_log_vec,:,:);
    cond3= test_trial_data(lure_2_log_vec,:,:);
    cond4= test_trial_data(lure_4_log_vec,:,:);
    cond5= test_trial_data(lure_5_log_vec,:,:);
    cond6= test_trial_data(new_log_vec,:,:);
    
elseif strcmp('tuning_correct',exp_type)
    
    cond1= test_trial_data(repeat_log_vec&corr_resp_log_vec ==1,:,:);
    cond2= test_trial_data(lure_all==1 & corr_resp_log_vec ==0,:,:);  %pattern comp
    cond3= test_trial_data(lure_all==1 & corr_resp_log_vec ==1,:,:); %lure
    cond4= test_trial_data(new_log_vec & corr_resp_log_vec ==1,:,:);
    

elseif strcmp('tuning_incorrect',exp_type)
    cond1= test_trial_data(repeat_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond2= test_trial_data(lure_1_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond3= test_trial_data(lure_2_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond4= test_trial_data(lure_4_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond5= test_trial_data(lure_5_log_vec==1&corr_resp_log_vec ==0,:,:);
    cond6= test_trial_data(new_log_vec==1&corr_resp_log_vec ==0,:,:);
    
elseif strcmp('indoor_outdoor',exp_type)
    cond1= study_trial_data(indoor_log_vec,:,:);
    cond2= study_trial_data(outdoor_log_vec,:,:);
    cond3= [];
    cond4= [];
    cond5= [];
    cond6= [];
    
elseif strcmp('encoding',exp_type)
    all_lure_log_vec = training_behav_matrix(:,2)>=1;
    
    % get the image list of lures during encoding
   lure_images = train_images(all_lure_log_vec,:);
  
    % find trial num of each lure during test
   for img = 1:size(lure_images,1)
       temp = lure_images(img,1:3);
       img_idx_during_test(img) = find(sum(test_images(:,1:3) == temp,2)==3);
   end
   
   lure_accuracy_log_mtx = zeros(size(all_lure_log_vec,1),2);
   lure_accuracy_log_mtx(:,1) = all_lure_log_vec;
   lure_accuracy_log_mtx(all_lure_log_vec==1,2) = corr_resp_log_vec(img_idx_during_test');
   
   cond1 = study_trial_data(lure_accuracy_log_mtx(:,1)==1 & lure_accuracy_log_mtx(:,2)==1,:,:); % lure --> +
   cond2 = study_trial_data(lure_accuracy_log_mtx(:,1)==1 & lure_accuracy_log_mtx(:,2)==0,:,:);  % lure --> -
   cond3= [];
   cond4= [];
   cond5= [];
   cond6= [];
elseif strcmp('cond_spec_prestim',exp_type)
    all_lure_log_vec = training_behav_matrix(:,2)>=1;
    
    % get the image list of lures during encoding
    lure_images = train_images(all_lure_log_vec,:);
    
    % find trial num of each lure during test
    for img = 1:size(lure_images,1)
        temp = lure_images(img,1:3);
        img_idx_during_test(img) = find(sum(test_images(:,1:3) == temp,2)==3);
    end
    
    lure_accuracy_log_mtx = zeros(size(all_lure_log_vec,1),2);
    lure_accuracy_log_mtx(:,1) = all_lure_log_vec;
    lure_accuracy_log_mtx(all_lure_log_vec==1,2) = corr_resp_log_vec(img_idx_during_test');
    
    
    cond3a_log_vec = lure_accuracy_log_mtx(:,1)==1 & lure_accuracy_log_mtx(:,2)==1; %-->lure+
    cond2a_log_vec = lure_accuracy_log_mtx(:,1)==1 & lure_accuracy_log_mtx(:,2)==0; %-->lure-
    cond1_log_vec  = repeat_log_vec&corr_resp_log_vec ==1;
    cond2_log_vec  = lure_all==1 & corr_resp_log_vec ==0;
    cond3_log_vec  = lure_all==1 & corr_resp_log_vec ==1;
    cond4_log_vec  = new_log_vec & corr_resp_log_vec ==1;  
    

end

% remove nan trials

   elec = 6; 

if ~isempty(cond1)
 cond1 = cond1(~isnan(cond1(:,1,elec)),:,:);
end
if ~isempty(cond2)
    cond2 = cond2(~isnan(cond2(:,1,elec)),:,:);
end
if ~isempty(cond3)
    cond3 = cond3(~isnan(cond3(:,1,elec)),:,:);
end
if ~isempty(cond4)
    cond4 = cond4(~isnan(cond4(:,1,elec)),:,:);
end
if ~isempty(cond5)
    cond5 = cond5(~isnan(cond5(:,1,elec)),:,:);
end
if ~isempty(cond6)
    cond6 = cond6(~isnan(cond6(:,1,elec)),:,:);
end

save(['cond_data_' lock '_' exp_type], 'cond1', 'cond2', 'cond3', 'cond4', 'cond5', 'cond6')

end




%% extra tuning corr conditions
    %     cond1= test_trial_data(repeat_log_vec&corr_resp_log_vec ==1,:,:);
    %     cond2= test_trial_data(lure_1_log_vec&corr_resp_log_vec ==1,:,:);
    %     cond3= test_trial_data(lure_2_log_vec&corr_resp_log_vec ==1,:,:);
    %     cond4= test_trial_data(lure_4_log_vec&corr_resp_log_vec ==1,:,:);
    %     cond5= test_trial_data(lure_5_log_vec&corr_resp_log_vec ==1,:,:);
    %     cond6= test_trial_data(new_log_vec,:,:);