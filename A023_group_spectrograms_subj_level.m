
if strcmp('response', lock)
data_length = 1601;
cd(['/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_response_' exp_type])

elseif strcmp('onset', lock)
data_length = 1251 %2101;
cd(['/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_onset_' exp_type])
end

fs = 1000;
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

tickmarks = 1:20:length(freq);

% get start and end trial indices
if strcmp ('onset',lock)
    pre_stim  = .5;
    post_stim = 2; % changed to 1 sec. ISI is 1 sec?
   
     
elseif strcmp ('offset',lock)
    pre_stim = 1;
    post_stim = .5;
   
end
edge_points = 200
cond_num = 3

% define num of chans for each patient
grp = {'39' '57' '44' '63' '66' '83' '84'}
  grp_fro_cond1 = nan(length(freq), data_length, length(grp));
  grp_fro_cond2 = nan(length(freq), data_length, length(grp));
  grp_fro_cond3= nan(length(freq), data_length, length(grp));
  grp_fro_cond4= nan(length(freq), data_length, length(grp));
  
grp = {'39' '57' '44' '63' '66' '83' '84'}
  grp_temp_cond1= nan(length(freq), data_length,length(grp));
  grp_temp_cond2= nan(length(freq), data_length, length(grp));
  grp_temp_cond3= nan(length(freq), data_length, length(grp));
  grp_temp_cond4= nan(length(freq), data_length, length(grp));
  
grp = {'39' '57' '44' '63' '66' '84'}
  grp_MTL_cond1= nan(length(freq), data_length, length(grp));
  grp_MTL_cond2= nan(length(freq), data_length, length(grp));
  grp_MTL_cond3= nan(length(freq), data_length, length(grp));
  grp_MTL_cond4= nan(length(freq), data_length, length(grp));

grp = {'39' '57' '44' '63' '66' '84'}
  grp_OFC_cond1= nan(length(freq), data_length, length(grp));
  grp_OFC_cond2= nan(length(freq), data_length, length(grp));
  grp_OFC_cond3= nan(length(freq), data_length, length(grp));
  grp_OFC_cond4= nan(length(freq), data_length, length(grp));
  
grp = {'39' '57' '44' '63' '66' }
  grp_CING_cond1= nan(length(freq), data_length, length(grp));
  grp_CING_cond2= nan(length(freq), data_length, length(grp));
  grp_CING_cond3= nan(length(freq), data_length, length(grp));
  grp_CING_cond4= nan(length(freq), data_length, length(grp));

grp = {'39' '57' '44'  '66' }
  grp_CA3_cond1= nan(length(freq), data_length, length(grp));
  grp_CA3_cond2= nan(length(freq), data_length, length(grp));
  grp_CA3_cond3= nan(length(freq), data_length, length(grp));
  grp_CA3_cond4= nan(length(freq), data_length, length(grp));
  
grp = {'39' '57' '44' '66' '84'}
  grp_CA1_cond1= nan(length(freq), data_length, length(grp));
  grp_CA1_cond2= nan(length(freq), data_length, length(grp));
  grp_CA1_cond3= nan(length(freq), data_length, length(grp));
  grp_CA1_cond4= nan(length(freq), data_length, length(grp));

grp = {'39' '57' '44' '66' }
  grp_HC_cond1 = nan(length(freq), data_length, length(grp));
  grp_HC_cond2 = nan(length(freq), data_length, length(grp));
  grp_HC_cond3 = nan(length(freq), data_length, length(grp));
  grp_HC_cond4 = nan(length(freq), data_length, length(grp));
 %% fro 
 cond1 = [];
 cond2 = [];
 cond3 = [];
 cond4 = [];
for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat']);
      grp_fro_cond1(:,:,a) = nanmean(fro_cond1,3);
      grp_fro_cond2(:,:,a) = nanmean(fro_cond2,3);
      grp_fro_cond3(:,:,a) = nanmean(fro_cond3,3);
      grp_fro_cond4(:,:,a) = nanmean(fro_cond4,3);
end
cond1 = grp_fro_cond1;
cond2 = grp_fro_cond2;
cond3 = grp_fro_cond3;
cond4 = grp_fro_cond4;

fro_grp_1_mn = nanmean(grp_fro_cond1,3);
fro_grp_1_std = nanstd(grp_fro_cond1,0,3);

fro_grp_2_mn = nanmean(grp_fro_cond2,3);
fro_grp_2_std = nanstd(grp_fro_cond2,0,3);

fro_grp_3_mn = nanmean(grp_fro_cond3,3);
fro_grp_3_std = nanstd(grp_fro_cond3,0,3);

fro_grp_4_mn = nanmean(grp_fro_cond4,3);
fro_grp_4_std = nanstd(grp_fro_cond4,0,3);


save('FRO_grp_spectrograms','cond1','cond2','cond3','cond4')

%% temp
 cond1 = [];
 cond2 = [];
 cond3 = [];
 cond4 = [];
 grp = {'39' '57' '44' '63' '66' '83' '84'}
for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
      grp_temp_cond1(:,:,a) = nanmean(temp_cond1,3);
      grp_temp_cond2(:,:,a) = nanmean(temp_cond2,3);
      grp_temp_cond3(:,:,a) = nanmean(temp_cond3,3);
      grp_temp_cond4(:,:,a) = nanmean(temp_cond4,3);
end
cond1 = grp_temp_cond1;
cond2 = grp_temp_cond2;
cond3 = grp_temp_cond3;
cond4 = grp_temp_cond4;

temp_grp_1_mn = nanmean(grp_temp_cond1,3);
temp_grp_1_std = nanstd(grp_temp_cond1,0,3);

temp_grp_2_mn = nanmean(grp_temp_cond2,3);
temp_grp_2_std = nanstd(grp_temp_cond2,0,3);

temp_grp_3_mn = nanmean(grp_temp_cond3,3);
temp_grp_3_std = nanstd(grp_temp_cond3,0,3);

temp_grp_4_mn = nanmean(grp_temp_cond4,3);
temp_grp_4_std = nanstd(grp_temp_cond4,0,3);

save('TEMP_grp_spectrograms','cond1','cond2','cond3','cond4')

%% MTL
 cond1 = [];
 cond2 = [];
 cond3 = [];
 cond4 = [];
 grp = {'39' '57' '44' '63' '66' '84'}

for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])  
      grp_MTL_cond1(:,:,a) = nanmean(MTL_cond1,3);
      grp_MTL_cond2(:,:,a) = nanmean(MTL_cond2,3);
      grp_MTL_cond3(:,:,a) = nanmean(MTL_cond3,3);
      grp_MTL_cond4(:,:,a) = nanmean(MTL_cond4,3);
end
cond1 = grp_MTL_cond1;
cond2 = grp_MTL_cond2;
cond3 = grp_MTL_cond3;
cond4 = grp_MTL_cond4;

MTL_grp_1_mn = nanmean(grp_MTL_cond1,3);
MTL_grp_1_std = nanstd(grp_MTL_cond1,0,3);

MTL_grp_2_mn = nanmean(grp_MTL_cond2,3);
MTL_grp_2_std = nanstd(grp_MTL_cond2,0,3);

MTL_grp_3_mn = nanmean(grp_MTL_cond3,3);
MTL_grp_3_std = nanstd(grp_MTL_cond3,0,3);

MTL_grp_4_mn = nanmean(grp_MTL_cond4,3);
MTL_grp_4_std = nanstd(grp_MTL_cond4,0,3);
save('MTL_grp_spectrograms','cond1','cond2','cond3','cond4')


%% CA3
 cond1 = [];
 cond2 = [];
 cond3 = [];
 cond4 = [];
 
 grp = {'39' '57' '44' '66'}

for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
      grp_CA3_cond1(:,:,a) = nanmean(CA3_cond1,3);
      grp_CA3_cond2(:,:,a) = nanmean(CA3_cond2,3);
      grp_CA3_cond3(:,:,a) = nanmean(CA3_cond3,3);
      grp_CA3_cond4(:,:,a) = nanmean(CA3_cond4,3);
end
cond1 = grp_CA3_cond1;
cond2 = grp_CA3_cond2;
cond3 = grp_CA3_cond3;
cond4 = grp_CA3_cond4;

CA3_grp_1_mn = nanmean(grp_CA3_cond1,3);
CA3_grp_1_std = nanstd(grp_CA3_cond1,0,3);

CA3_grp_2_mn = nanmean(grp_CA3_cond2,3);
CA3_grp_2_std = nanstd(grp_CA3_cond2,0,3);

CA3_grp_3_mn = nanmean(grp_CA3_cond3,3);
CA3_grp_3_std = nanstd(grp_CA3_cond3,0,3);

CA3_grp_4_mn = nanmean(grp_CA3_cond4,3);
CA3_grp_4_std = nanstd(grp_CA3_cond4,0,3);

save('CA3_grp_spectrograms','cond1','cond2','cond3','cond4')


%% CA1
 cond1 = [];
 cond2 = [];
 cond3 = [];
 cond4 = [];
  grp = {'39' '57' '44' '66' '84'}
for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
      grp_CA1_cond1(:,:,a) = nanmean(CA1_cond1,3);
      grp_CA1_cond2(:,:,a) = nanmean(CA1_cond2,3);
      grp_CA1_cond3(:,:,a) = nanmean(CA1_cond3,3);
      grp_CA1_cond4(:,:,a) = nanmean(CA1_cond4,3);
end
cond1 = grp_CA1_cond1;
cond2 = grp_CA1_cond2;
cond3 = grp_CA1_cond3;
cond4 = grp_CA1_cond4;

CA1_grp_1_mn = nanmean(grp_CA1_cond1,3);
CA1_grp_1_std = nanstd(grp_CA1_cond1,0,3);

CA1_grp_2_mn = nanmean(grp_CA1_cond2,3);
CA1_grp_2_std = nanstd(grp_CA1_cond2,0,3);

CA1_grp_3_mn = nanmean(grp_CA1_cond3,3);
CA1_grp_3_std = nanstd(grp_CA1_cond3,0,3);

CA1_grp_4_mn = nanmean(grp_CA1_cond4,3);
CA1_grp_4_std = nanstd(grp_CA1_cond4,0,3);

save('CA1_grp_spectrograms','cond1','cond2','cond3','cond4')
%% ca1-3
 cond1 = [];
 cond2 = [];
 cond3 = [];
 cond4 = [];
 
CA13_1 = nan(size(grp_CA1_cond1));
CA13_2 = nan(size(grp_CA1_cond2));
CA13_3 = nan(size(grp_CA1_cond3));
CA13_4 = nan(size(grp_CA1_cond4));

for a = 1:length(grp)
    CA13_1(:,:,a) = nanmean(cat(3,grp_CA1_cond1(:,:,a), grp_CA3_cond1(:,:,a)),3);
    CA13_2(:,:,a) = nanmean(cat(3,grp_CA1_cond2(:,:,a), grp_CA3_cond2(:,:,a)),3);
    CA13_3(:,:,a) = nanmean(cat(3,grp_CA1_cond3(:,:,a), grp_CA3_cond3(:,:,a)),3);
    CA13_4(:,:,a) = nanmean(cat(3,grp_CA1_cond4(:,:,a), grp_CA3_cond4(:,:,a)),3);

end
cond1 = CA13_1;
cond2 = CA13_2;
cond3 = CA13_3;
cond4 = CA13_4;

CA13_grp_1_mn = nanmean(CA13_1,3);
CA13_grp_1_std = nanstd(CA13_1,0,3);

CA13_grp_2_mn = nanmean(CA13_2,3);
CA13_grp_2_std = nanstd(CA13_2,0,3);

CA13_grp_3_mn = nanmean(CA13_3,3);
CA13_grp_3_std = nanstd(CA13_3,0,3);

CA13_grp_4_mn = nanmean(CA13_4,3);
CA13_grp_4_std = nanstd(CA13_4,0,3);


save('CA13_grp_spectrograms','cond1','cond2','cond3','cond4')

%% OFC
 cond1 = [];
 cond2 = [];
 cond3 = [];
 cond4 = [];
for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
  
      grp_OFC_cond1(:,:,a) = nanmean(OFC_cond1,3);
      grp_OFC_cond2(:,:,a) = nanmean(OFC_cond2,3);
      grp_OFC_cond3(:,:,a) = nanmean(OFC_cond3,3);
      grp_OFC_cond4(:,:,a) = nanmean(OFC_cond4,3);

end
cond1 = grp_OFC_cond1;
cond2 = grp_OFC_cond2;
cond3 = grp_OFC_cond3;
cond4 = grp_OFC_cond4;

OFC_grp_1_mn  = nanmean(grp_OFC_cond1,3);
OFC_grp_1_std = nanstd(grp_OFC_cond1,0,3);

OFC_grp_2_mn  = nanmean(grp_OFC_cond2,3);
OFC_grp_2_std = nanstd(grp_OFC_cond2,0,3);

OFC_grp_3_mn  = nanmean(grp_OFC_cond3,3);
OFC_grp_3_std = nanstd(grp_OFC_cond3,0,3);

OFC_grp_4_mn  = nanmean(grp_OFC_cond4,3);
OFC_grp_4_std = nanstd(grp_OFC_cond4,0,3);

save('OFC_grp_spectrograms','cond1','cond2','cond3','cond4')

%% CING
 cond1 = [];
 cond2 = [];
 cond3 = [];
 cond4 = [];
 grp = {'39' '57' '44' '63' '66' }

for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
      grp_CING_cond1(:,:,a) = nanmean(cing_cond1,3);
      grp_CING_cond2(:,:,a) = nanmean(cing_cond2,3);
      grp_CING_cond3(:,:,a) = nanmean(cing_cond3,3);
      grp_CING_cond4(:,:,a) = nanmean(cing_cond4,3);
end
cond1 = grp_CING_cond1;
cond2 = grp_CING_cond2;
cond3 = grp_CING_cond3;
cond4 = grp_CING_cond4;

CING_grp_1_mn = nanmean(grp_CING_cond1,3);
CING_grp_1_std = nanstd(grp_CING_cond1,0,3);

CING_grp_2_mn = nanmean(grp_CING_cond2,3);
CING_grp_2_std = nanstd(grp_CING_cond2,0,3);

CING_grp_3_mn = nanmean(grp_CING_cond3,3);
CING_grp_3_std = nanstd(grp_CING_cond3,0,3);

CING_grp_4_mn = nanmean(grp_CING_cond4,3);
CING_grp_4_std = nanstd(grp_CING_cond4,0,3);

save('CING_grp_spectrograms','cond1','cond2','cond3','cond4')

%% plot
cond_num = 4
edge_points = 200
mx = .3
mn = -.3

figure
suptitle('frontal')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_1_mn./fro_grp_1_std)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_1_mn)

set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_2_mn./fro_grp_2_std)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_2_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_3_mn./fro_grp_3_std)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_3_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), fro_grp_4_mn./fro_grp_4_std)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_4_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

saveas(gcf, 'fro_onset.bmp') 





%% temporal

mx = 1.7
mn = -1

figure
suptitle('temporal')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(temp_grp_1_mn,2)), 1:length(freq), temp_grp_1_mn./temp_grp_1_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), temp_grp_2_mn./temp_grp_2_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), temp_grp_3_mn./temp_grp_3_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), temp_grp_4_mn./temp_grp_4_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

saveas(gcf, 'temp_onset_theta.bmp') 



%% MTL


mx = 1.5
mn = -.2

figure
suptitle('MTL')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_1_mn./MTL_grp_1_std)
%imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_1_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_2_mn./MTL_grp_2_std)
%imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_2_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_3_mn./MTL_grp_3_std)
%imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_3_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_4_mn./MTL_grp_4_std)
%imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_4_mn)

set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

%saveas(gcf, 'MTL.bmp') 

%% CA3

mx = 1.5
mn = -1

figure
suptitle('CA3 - 4 patients')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA3_grp_1_mn./CA3_grp_1_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq),  CA3_grp_2_mn./CA3_grp_2_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('pattern Comp: incorr lure')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA3_grp_3_mn./CA3_grp_3_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('Pattern Sep: corr lure')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA3_grp_4_mn./CA3_grp_4_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

saveas(gcf, 'ca3.bmp') 

%%  CA1


mx = 1.5
mn = -1

figure
suptitle('CA1 - 5 patients')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA1_grp_1_mn./CA1_grp_1_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq),  CA1_grp_2_mn./CA1_grp_2_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('pattern Comp: incorr lure')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA1_grp_3_mn./CA1_grp_3_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('Pattern Sep: corr lure')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA1_grp_4_mn./CA1_grp_4_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
saveas(gcf, 'ca1.bmp') 

%% CA1-CA3-DG

mx = 1.5
mn = -1.5

figure
suptitle('HC: CA1-CA3-DG')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA13_grp_1_mn./CA1_grp_4_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq),  CA13_grp_2_mn./CA1_grp_2_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure - ')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA13_grp_3_mn./CA1_grp_3_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure + ')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA13_grp_4_mn./CA1_grp_4_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

saveas(gcf, 'ca1-3dg.bmp') 

%% OFC

mx = 2
mn = -1

figure
suptitle('OFC - 6 patients')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), OFC_grp_1_mn./OFC_grp_1_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), OFC_grp_2_mn./OFC_grp_2_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), OFC_grp_3_mn./OFC_grp_3_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), OFC_grp_4_mn./OFC_grp_4_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')


saveas(gcf, 'ofc.bmp') 


%% CING

mx= 1.82
mn = -1

figure
suptitle('CING')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CING_grp_1_mn./CING_grp_1_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CING_grp_2_mn./CING_grp_2_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('pattern Comp: incorr lure')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CING_grp_3_mn./CING_grp_3_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('Pattern Sep: corr lure')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), CING_grp_4_mn./CING_grp_4_std)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

saveas(gcf, 'cing.bmp') 

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% mean only%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cond_num = 4
edge_points = 200
mx = .32
mn = -.3

figure
suptitle('frontal')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_1_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_2_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), fro_grp_3_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), fro_grp_4_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

saveas(gcf, 'fro_onset.bmp') 

%%

mx = .30
mn = -.3

figure
suptitle('MTL')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), MTL_grp_1_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('repeat')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), MTL_grp_2_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,3)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), MTL_grp_3_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
colormap jet
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,4)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(indiv_freq_chans_cond1,2)), 1:length(freq), MTL_grp_4_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('new')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

