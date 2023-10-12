clear all;close all;clc

% this script saves the behavior log data into matrices
% output: testing_behav_matrix, training_behav_matrix

subj = 7
patient_num = [ 39 57  66 63 44 83 84 85 87] ;
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
cd(['/mnt/yassamri/iEEG/sandra/subj_' num2str(patient_num(subj))]) % go to patient log dir
filename = [num2str(patient_num(subj)) '_MDTO_log.txt']
MDTOlog  = A000_import_log(filename)
%%
counter = 0;
for a = 1:size(MDTOlog,1) % loop thru throw
  for b = 1:size(MDTOlog,2) % loop thru columns
    if  strcmp('ImageType',MDTOlog{a,b})
      counter = counter+1
      row(counter,:) =[a b]
    end
    
    if strcmp('Scores:', MDTOlog{a,b})
      row_end = [a b]
    end
  end
  
  
end
%% start and end rows for training and testing tasks
Training_trial_row = [row(1,1)+1 row(2,1)-2];
Testing_trial_row = [row(2,1)+1 row_end(1,1)-1];


%%
%trial image imagetype responses rt
trial_length                = length(Training_trial_row(1):Training_trial_row(2));
train_behav_data            = cell(1,2);

% get trial num
train_behav_data{1}(1:trial_length,1) = MDTOlog{Training_trial_row(1),1}:MDTOlog{Training_trial_row(2),1};

% images
for a = 1:trial_length
  counter                    = Training_trial_row(1)-1+a
  
  % 2nd row: get image
  train_behav_data{2}(a,:)   = MDTOlog{counter,2}
  
  % 3rd row: get image type
  if strcmp('sR', MDTOlog{counter,3})    
    train_behav_data{3}(a,:)   = 0 ; % convert sR to 0
  else
    train_behav_data{3}(a,:)   = MDTOlog{counter,3};
  end
  
  % 4th row: get response
  if ischar(MDTOlog{counter,4})
     train_behav_data{4}(a,:) = nan;
  else
      train_behav_data{4}(a,:) = MDTOlog{counter,4};
  end
  
  % 5th row: get response time
  if isempty(MDTOlog{counter,5})
     train_behav_data{5}(a,:)  = 0;
  else
  train_behav_data{5}(a,:) = MDTOlog{counter,5};
  end
  
  
end

train_trials     =  (MDTOlog{Training_trial_row(1),1}:MDTOlog{Training_trial_row(2),1})';
train_images     = train_behav_data{2};
train_image_type = train_behav_data{3};
train_response   = train_behav_data{4};
train_RT         = train_behav_data{5};
training_behav_matrix(:,1) = train_trials;
training_behav_matrix(:,2) = train_image_type;
training_behav_matrix(:,3) = train_response;
training_behav_matrix(:,4) = train_RT;
training_behav_matrix_label = {'trials' 'image_type' 'response' 'RT'}
%%

trial_length                = length(Testing_trial_row(1):Testing_trial_row(2));
test_behav_data            = cell(1,2);

% get trial num
test_behav_data{1}(1:trial_length,1) = MDTOlog{Testing_trial_row(1),1}:MDTOlog{Testing_trial_row(2),1};

% images
for a = 1:trial_length
  counter                    = Testing_trial_row(1)-1+a
  
  % 2nd row: get image
  test_behav_data{2}(a,:)   = MDTOlog{counter,2}
  

  % 3rd row: get image type
  if strcmp('sR', MDTOlog{counter,3})    
    test_behav_data{3}(a,:)   = 0 ; % convert sR to 0
  elseif strcmp('sF',MDTOlog{counter,3})
        test_behav_data{3}(a,:)   = 0.5 ; % convert sF to 0.5

  else 
    test_behav_data{3}(a,:)   = MDTOlog{counter,3};
  end
  
  
  % 5th row: get response
  if ischar(MDTOlog{counter,5})
     test_behav_data{4}(a,:) = nan;
  else
      test_behav_data{4}(a,:) = MDTOlog{counter,5};
  end


  % 6th row: get response time
  if isempty(MDTOlog{counter,6})
     test_behav_data{5}(a,:)  = 0;
  else
  test_behav_data{5}(a,:) = str2num(MDTOlog{counter,6});
  end
  
    % 4th row: get response
  if ischar(MDTOlog{counter,4})
     test_behav_data{6}(a,:) = nan;
  else
      test_behav_data{6}(a,:) = MDTOlog{counter,4};
  end
  
end

test_trials     =  (MDTOlog{Testing_trial_row(1),1}:MDTOlog{Testing_trial_row(2),1})';
test_images     = test_behav_data{2};
test_image_type = test_behav_data{3};
test_response   = test_behav_data{4};
test_RT         = test_behav_data{5};
corr_response   = test_behav_data{6};


testing_behav_matrix(:,1) = test_trials;
testing_behav_matrix(:,2) = test_image_type;
testing_behav_matrix(:,3) = test_response;
testing_behav_matrix(:,4) = test_RT;
testing_behav_matrix(:,5) = corr_response;
testing_behav_matrix_label = {'trials' 'image_type' 'response' 'RT' 'corr_response'}


%% 
save(['behavior_subj' num2str(patient_num(subj))], 'training_behav_matrix','training_behav_matrix_label','train_images',...
      'testing_behav_matrix', 'testing_behav_matrix_label', 'test_images')
    clear all
