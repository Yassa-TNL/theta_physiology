clear all;close all;clc
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')

exp_type      = 'tuning_correct' %encoding, tuning_correct
lock          = 'onset'   %onset response
ref           = 'LM'
fs            = 500;
baseline      = '_prestim' % '': entire recording, 'pre_stim'
desired_freq  = 5.5
reg_list      = {'HC'  'OFC' 'FRO'  'TEMP' 'CING' 'ins' 'EC'}

minfreq = 3; maxfreq = 200;
[freq] = get_freq(fs,minfreq, maxfreq);
tickmarks = 1:21:length(freq);
if strcmp('encoding', exp_type);
    cond_num=2;    titles = {'lure+','lure-'}
else
    cond_num = 4; titles = {'repeat', 'lure-', 'lure+', 'new'}
end

if strcmp('response', lock)
elseif strcmp('onset', lock)
    pre_stim = 0.5; post_stim = 2.0;
end
if strcmp('encoding', exp_type)
    time_range = [.3*fs (pre_stim+post_stim)*fs];
    current_pre_stim  = -.2;
    current_post_stim = 2;
else
    time_range = [.3*fs 1.5*fs]; % was [301 1301];
    current_pre_stim =-.2;
    current_post_stim= 1; % was 1;
end

figure;
hold on;
color_list = {'r', [1 0.5 0.4], 'k', [0 0 1], [0.3 0.4 .7], [0.6 0.7 .8], [0 1 1], [0 0.8 .2]}

for iReg = 1:length(reg_list)

reg           = reg_list{iReg}
[subj_list]   = get_subj_list(reg);

% get pooled elecs from all subj
[cond1a,cond2a,cond3a,cond4a] = get_power_subj_elecs(subj_list,exp_type,reg,ref,baseline,lock);

cond3 = cat(3,cond3a{:});
trace3 = squeeze(nanmean(cond3(freq<desired_freq, time_range(1):time_range(2), :)))';

% plot traces and stats
stdshade(trace3,.1, color_list{iReg},linspace(current_pre_stim,current_post_stim, size(trace3,2)),[] ,[], []);
h(iReg) = plot(linspace(current_pre_stim,current_post_stim, size(trace3,2)),nanmean(trace3,1),'Color',color_list{iReg}, 'LineWidth', 3);
xlim([current_pre_stim current_post_stim]);
xlabel('Time (s)'), ylabel('3-5 Hz power')
end

ylim([-0.1 0.35]);
y = ylim;
line([0 0],ylim,'color','k')
xticks([0  1])
yticks([0 y(2)])

%% test colors systematically
color_list = {'r', [1 0.5 0.4], 'k', [0 0 1], [0.3 0.4 .7], [0.6 0.7 .8], [0 1 1], [0 0.8 .2]}
figure;hold on
for iReg=1:8
plot(iReg*nanmean(trace3,1),'Color',color_list{iReg}, 'LineWidth', 2);
end
