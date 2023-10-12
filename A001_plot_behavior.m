% plot behavior
clear all; close all; clc
patient_num = [39 44 57 63 66  84 85 87] ;%83


%% plot test behav data
 figure;hold on
% plot response distributions and over all accuracy for all patients
for subj = 1:length(patient_num)
   
    cd(['/mnt/yassamri/iEEG/sandra/subj_'  num2str(patient_num(subj))])
    load(['behavior_subj' num2str(patient_num(subj)) '.mat'])
    histogram((testing_behav_matrix(:,4)), 'normalization', 'probability', 'Binwidth', 0.06)
    % accuracy for entire session
    percent_corr (subj) = sum(testing_behav_matrix(:,3)==testing_behav_matrix(:,5))/length(testing_behav_matrix(:,3))
    %legend({[num2str(percent_corr(subj)*100) '% corr']})
    ylabel('probability')
    xlabel('response time')
    set(gca, 'FontSize', 18,'FontWeight','bold')
end

%% 
%% plot test behav data
 figure;hold on
% plot response distributions and over all accuracy for all patients
for subj = 1:length(patient_num)
   
    cd(['/mnt/yassamri/iEEG/sandra/subj_'  num2str(patient_num(subj))])
    load(['behavior_subj' num2str(patient_num(subj)) '.mat'])
    histogram((testing_behav_matrix(:,4)), 'normalization', 'probability', 'Binwidth', 0.06)
    % accuracy for entire session
    percent_corr (subj) = sum(testing_behav_matrix(:,3)==testing_behav_matrix(:,5))/length(testing_behav_matrix(:,3))
    %legend({[num2str(percent_corr(subj)*100) '% corr']})
    ylabel('probability')
    xlabel('response time')
    set(gca, 'FontSize', 18,'FontWeight','bold')
    clear testing_behav_matrix 
end
%% plot response distributions using pooled data
pooled_RT_data =[];
for subj =1:length(patient_num)
    cd(['/mnt/yassamri/iEEG/sandra/subj_'  num2str(patient_num(subj))])
    load(['behavior_subj' num2str(patient_num(subj))])
    pooled_RT_data =[ pooled_RT_data; testing_behav_matrix(:,4)]
end
figure;
histogram(pooled_RT_data, 'normalization', 'probability','FaceColor', 'k')
ylabel('probability')
xlabel('response time')
title('RT in Tesing Phase')
set(gca, 'FontSize', 12,'FontWeight','bold')
  


%% plot trainig response and RT histograms

subj = 84;
load(['subj' num2str(patient_num(subj)) '.mat'])

% responses
figure
subplot(2,1,1)
hist(training_behav_matrix(:,3))
xlabel('response')
ylabel('couonts')
title(['subj ' num2str(patient_num(subj)) ' response'])
set(gca, 'FontSize',14,'FontWeight','bold')

% response times
subplot(2,1,2)
hist(training_behav_matrix(:,4))
xlabel('response time')
ylabel('counts')
title(['subj ' num2str(patient_num(subj)) ' response time'])
set(gca, 'FontSize',14,'FontWeight','bold')


%% plot conditional accuracies
for subj = 1:length(patient_num)
    cd(['/mnt/yassamri/iEEG/sandra/subj_'  num2str(patient_num(subj))])
    load(['behavior_subj' num2str(patient_num(subj)) '.mat'])  
  % find repeats
  repeat_idx = find(testing_behav_matrix(:,2)==0)
  repeat_responses_patient = testing_behav_matrix(repeat_idx,3)
  repeat_responses_corr = testing_behav_matrix(repeat_idx,5)
  repeat_response_accuracy (subj) = sum(repeat_responses_patient==repeat_responses_corr)/length(repeat_idx)*100;
  
  
  
  % find lures
  lure_idx = find(testing_behav_matrix(:,2)>=1)
  lure_responses_patient = testing_behav_matrix(lure_idx,3)
  lure_responses_corr = testing_behav_matrix(lure_idx,5)
  lure_response_accuracy (subj) = sum(lure_responses_patient==lure_responses_corr)/length(lure_idx)*100;
  
  % hi sim
  lure_hi_idx = find(testing_behav_matrix(:,2) == or(1,2) )
  lure_hi_responses_patient = testing_behav_matrix(lure_hi_idx,3)
  lure_hi_responses_corr = testing_behav_matrix(lure_hi_idx,5)
  lure_hi_sim_response_accuracy (subj) = sum(lure_hi_responses_patient==lure_hi_responses_corr)/length(lure_hi_idx)*100;
  
  % mid sim
  lure_mid_idx = find(testing_behav_matrix(:,2) == 3 )
  lure_mid_responses_patient = testing_behav_matrix(lure_mid_idx,3)
  lure_mid_responses_corr = testing_behav_matrix(lure_mid_idx,5)
  lure_mid_sim_response_accuracy (subj) = sum(lure_mid_responses_patient==lure_mid_responses_corr)/length(lure_mid_idx)*100;
  
  
  %lo sim
  
  lure_lo_idx = find(round(testing_behav_matrix(:,2))==or(4,5))
  lure_lo_responses_patient = testing_behav_matrix(lure_lo_idx,3)
  lure_lo_responses_corr = testing_behav_matrix(lure_lo_idx,5)
  lure_lo_sim_response_accuracy (subj) = sum(lure_lo_responses_patient==lure_lo_responses_corr)/length(lure_lo_idx)*100;
  
  
  
  
  
  % find foil (new)
  foil_idx = find(testing_behav_matrix(:,2)==0.5)
  foil_responses_patient = testing_behav_matrix(foil_idx,3)
  foil_responses_corr = testing_behav_matrix(foil_idx,5)
  foil_response_accuracy (subj) = sum(foil_responses_patient==foil_responses_corr)/length(foil_idx)*100 ;
  
  
end


%%
figure
hold on


plot(ones(length(percent_corr),1),percent_corr*100, 'ok','MarkerSize', 10)
plot(2*ones(length(percent_corr),1),repeat_response_accuracy, 'or', 'MarkerSize', 10)
plot(3*ones(length(percent_corr),1),foil_response_accuracy, 'og', 'MarkerSize', 10)
plot(4*ones(length(percent_corr),1),lure_lo_sim_response_accuracy, 'om', 'MarkerSize', 10)
plot(5*ones(length(percent_corr),1),lure_response_accuracy, 'om', 'MarkerSize', 10)
plot(6*ones(length(percent_corr),1),lure_hi_sim_response_accuracy, 'om', 'MarkerSize', 10)
set(gca, 'XTick', 1:subj) %, 'XTickLabel')
ylabel('% accuracy')
title('Subject Performance')
xlim([0 7])
set(gca, 'FontSize',17,'FontWeight','bold')
cd('/mnt/yassamri/iEEG/sandra/GroupFigures')


%% indiv subj performance plots
labels = { 'Overall ','Repeat ','Lure ','Hi Sim Lure ','Lo Sim Lure ','Foil '}

for subj =6% 1:length(patient_num)
  figure
  bar([percent_corr(subj)*100 repeat_response_accuracy(subj) lure_response_accuracy(subj)...
    lure_hi_sim_response_accuracy(subj) lure_lo_sim_response_accuracy(subj) foil_response_accuracy(subj)] , 'm')
  set(gca,'XTickLabel', labels)
  xlabel('condition')
  ylabel('% accuracy')
  title(['Subject ' num2str(patient_num(subj)) ' Performance'])
  set(gca, 'FontSize',13,'FontWeight','bold')
  
end