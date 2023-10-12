clear all
regList = { 'OFC' 'FRO'  'TEMP'  'MTL'}
regList={'MTL'}
%%
% go to dir/load data
cd('/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_onset')
addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-20181130')

for reg_counter = 1:length(regList)
    reg  = regList{reg_counter}
    load([reg '_grp_spectrograms.mat'])
   
    % rearrange data
    conda=cond2(:,1:1301,:);
    condb=cond3(:,1:1301,:);


% smoothen average gamma power for each subject
fs = 1000;
minfreq = 3;
maxfreq = 200;
edge_points = 200;
dt = 1/fs;
NumVoices = 32;
a0 = 2^(1/NumVoices);
wavCenterFreq = 6/(2*pi);
minscale = wavCenterFreq/(maxfreq*dt); 
maxscale = wavCenterFreq/(minfreq*dt); 
minscale = floor(NumVoices*log2(minscale));
maxscale = ceil(NumVoices*log2(maxscale)); 
scales = a0.^(minscale:maxscale).*dt;
freq = wavCenterFreq./(fs*scales.*dt);
win= .135*fs;

for e = 1:size(conda,3)
    for gamma_freq = find(freq>40)
        conda(gamma_freq,:,e) = conv(conda(gamma_freq,:,e), ones(1,win)/win,'same');
        condb(gamma_freq,:,e) = conv(condb(gamma_freq,:,e), ones(1,win)/win,'same');

    end
end
%%rearrange data
cond2 = permute(conda, [3  1 2]);
cond3 = permute(condb, [3  1 2]);


pre_stim  =.5;
post_stim = 2;

% t-test cluster

% cond 1 definition
wDataCond2.label  = {'Region'};
wDataCond2.dimord = 'rpt_chan_freq_time';
wDataCond2.freq   = freq;
wDataCond2.time   = linspace(-.3,1,size(cond3,3));
wDataCond2.cfg    = [];
wDataCond2.powspctrm = nan(size(cond2, 1),1,size(cond2, 2),size(cond2, 3)); % subj X 1 X freq X time
wDataCond2.cumtapcnt = ones(size(cond2, 1), length(freq));
for a = 1:size(cond2,1)
    wDataCond2.powspctrm(a, 1, :, :) = cond2(a,:,:) ;
end

% cond 2 definition
wDataCond3.label  = {'Region'};
wDataCond3.dimord = 'rpt_chan_freq_time';
wDataCond3.freq   = freq;
wDataCond3.time   = linspace(-.3,1,size(cond3,3));
wDataCond3.cfg    = [];
wDataCond3.powspctrm = nan(size(cond3, 1),1,size(cond3, 2),size(cond3, 3)); % subj X cfreq X time
wDataCond3.cumtapcnt = ones(size(cond3, 1), length(freq));
for a = 1:size(cond3,1)
    wDataCond3.powspctrm(a, 1, :, :) = cond3(a,:,:) ;
end

% cfg definition
cfg = [];
cfg.channel          = {'Region'};
cfg.frequency        = 'all';
cfg.avgoverchan      = 'no';
cfg.avgovertime      = 'no';
cfg.avgoverfreq      = 'no';
cfg.parameter        = 'powspctrm';
cfg.correctm         = 'cluster';
cfg.neighbours.label = [];
cfg.neighbours.neighblabel = [];
cfg.method           = 'montecarlo';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.clusteralpha     = 0.05;
cfg.clusterstatistic = 'maxsum';
cfg.clusterthreshold = 'parametric'
cfg.tail             = 1; % 0 for both reight and left sided
cfg.clustertail      = 0;
cfg.alpha            = 0.05;
cfg.numrandomization = 1000;
design = [ones(1, size(cond3, 1)), ones(1,size(cond3, 1))*2];
cfg.design           = design;
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, wDataCond3, wDataCond2);
if ~isempty(find(stat.prob <.05))
    fprintf('\nFound significant cluster!\n')
else
    fprintf('\nno significant clusters\n')
end
sigIdx   = squeeze(stat.prob(1,:,:)) <.05;


figure
tickmarks = 1:20:length(freq);
% mn = -.3
% mx = .3
suptitle ([' Reg:' reg])
contourf(linspace(0,1,size(conda,2)), 1:length(freq),nanmean(condb,3),40,'linecolor','none')
%set(gca, 'YDir', 'reverse')
hold on
contour(linspace(0,1,size(conda,2)), 1:length(freq),flipud(sigIdx),1,'linecolor','k')
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)), 'YDir', 'reverse')
colorbar
%caxis([mn mx])
xlabel('time')
ylabel('freq') 

end