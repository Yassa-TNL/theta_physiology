function [cond3a_log_vec,cond2a_log_vec,cond1_log_vec,cond2_log_vec,cond3_log_vec,cond4_log_vec] = get_trial_condition_label(subj)
%UNTITLED3 Summary of this function goes here
% load behavior
load(['behavior_subj' subj '.mat']);

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
all_lure_log_vec   = training_behav_matrix(:,2)>=1;

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

% all conds of interest
cond3a_log_vec = lure_accuracy_log_mtx(:,1)==1 & lure_accuracy_log_mtx(:,2)==1; %-->lure+
cond2a_log_vec = lure_accuracy_log_mtx(:,1)==1 & lure_accuracy_log_mtx(:,2)==0; %-->lure-
cond1_log_vec  = repeat_log_vec&corr_resp_log_vec ==1;
cond2_log_vec  = lure_all==1 & corr_resp_log_vec ==0;
cond3_log_vec  = lure_all==1 & corr_resp_log_vec ==1;
cond4_log_vec  = new_log_vec & corr_resp_log_vec ==1;
end

