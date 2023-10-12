clear all;close all;clc
subj = '87'
lock = 'response'
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
load(['variables_sync_neurlynx_clinical_' lock '.mat'])

% pair research and clin labels
elec_pairs = [];
for iElec =1:length(research_chan_label)
    for iElec2 = 1:length(clinical_chan_label)
        if strcmp(clinical_chan_label{iElec2},research_chan_label{iElec})
            elec_pairs = [ elec_pairs ; iElec iElec2];
        end
    end
end

if strcmp('84',subj)
    yshift  = -16^-4;
    yscale = -10^-6;
    D = 2;
elseif strcmp('85',subj)
    yshift  = 10^-5;
    yscale = -10^-6;
    D = 1;
elseif strcmp('87',subj) % onset good
    yshift  = 10^-5;
    yscale = -10^-6;
    D = 1;
end

% plot reref for all chans
trial = 1;
for iElec = 1
    figure;plot(clinical_trial_data(trial,:,elec_pairs(iElec, 2)) - clinical_trial_data(trial,:,elec_pairs(iElec+1, 2)))
    hold on ;
    plot(decimate(double(research_trial_data(trial,:,elec_pairs(iElec, 1)) - research_trial_data(trial,:,elec_pairs(iElec+1, 1))),D)*yscale)
end


%% make data 
subj = '87'
lock = 'response'
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])

load(['trial_data_subj_' lock '_' subj '_ref__select_chan_3_research.mat'])
clearvars -except subj lock trial_data chan_label
research_trial_data = trial_data;
research_chan_label = chan_label;
clearvars -except subj lock research_trial_data research_chan_label

load(['trial_data_subj_' lock '_' subj '_ref__select_chan_3_clinical.mat'])
clinical_trial_data = trial_data;
clinical_chan_label = chan_label;
clearvars -except lock clinical_trial_data clinical_chan_label research_trial_data research_chan_label

save(['variables_sync_neurlynx_clinical_' lock '.mat'])
