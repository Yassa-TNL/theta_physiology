exp_type = 'encoding'
lock     = 'onset'
if strcmp('response', lock)
data_length = 1601;
cd(['/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_response_' exp_type])

elseif strcmp('onset', lock)
data_length = 2101;
cd(['/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_onset_' exp_type])
end

minfreq = 3;
maxfreq = 200;
fs=1000;
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
cond_num = 4

% define num of chans for each patient
grp = {'39' '57' '44' '63' '66' '83' '84'}
grp_fro_cond1 = nan(length(freq), data_length, length(grp));
grp_fro_cond2 = nan(length(freq), data_length, length(grp));


grp = {'39' '57' '44' '63' '66' '83' '84'}
grp_temp_cond1= nan(length(freq), data_length,length(grp));
grp_temp_cond2= nan(length(freq), data_length, length(grp));


grp = {'39' '57' '44' '63' '66' '84'}
grp_MTL_cond1= nan(length(freq), data_length, length(grp));
grp_MTL_cond2= nan(length(freq), data_length, length(grp));


grp = {'39' '57' '44' '63' '66' '84'}
grp_OFC_cond1= nan(length(freq), data_length, length(grp));
grp_OFC_cond2= nan(length(freq), data_length, length(grp));
  
grp = {'39' '57' '44' '63' '66' }
grp_CING_cond1= nan(length(freq), data_length, length(grp));
grp_CING_cond2= nan(length(freq), data_length, length(grp));


grp = {'39' '57' '44'  '66' }
grp_CA3_cond1= nan(length(freq), data_length, length(grp));
grp_CA3_cond2= nan(length(freq), data_length, length(grp));
  
grp = {'39' '57' '44' '66' '84'}
grp_CA1_cond1= nan(length(freq), data_length, length(grp));
grp_CA1_cond2= nan(length(freq), data_length, length(grp));

grp = {'39' '57' '44' '66' }
grp_HC_cond1 = nan(length(freq), data_length, length(grp));
grp_HC_cond2 = nan(length(freq), data_length, length(grp));


%%
cond1 = [];
cond2 = [];

for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat']);
      grp_fro_cond1(:,:,a) = nanmean(fro_cond1,3);
      grp_fro_cond2(:,:,a) = nanmean(fro_cond2,3);

end
cond1 = grp_fro_cond1;
cond2 = grp_fro_cond2;


fro_grp_1_mn = nanmean(grp_fro_cond1,3);
fro_grp_1_std = nanstd(grp_fro_cond1,0,3);

fro_grp_2_mn = nanmean(grp_fro_cond2,3);
fro_grp_2_std = nanstd(grp_fro_cond2,0,3);

save('FRO_grp_spectrograms','cond1','cond2')

%% temp
 cond1 = [];
 cond2 = [];

 grp = {'39' '57' '44' '63' '66' '83' '84'}
for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
      grp_temp_cond1(:,:,a) = nanmean(temp_cond1,3);
      grp_temp_cond2(:,:,a) = nanmean(temp_cond2,3);

end
cond1 = grp_temp_cond1;
cond2 = grp_temp_cond2;


temp_grp_1_mn = nanmean(grp_temp_cond1,3);
temp_grp_1_std = nanstd(grp_temp_cond1,0,3);

temp_grp_2_mn = nanmean(grp_temp_cond2,3);
temp_grp_2_std = nanstd(grp_temp_cond2,0,3);


save('TEMP_grp_spectrograms','cond1','cond2')

%% MTL
 cond1 = [];
 cond2 = [];

 grp = {'39' '57' '44' '63' '66' '84'}

for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])  
      grp_MTL_cond1(:,:,a) = nanmean(MTL_cond1,3);
      grp_MTL_cond2(:,:,a) = nanmean(MTL_cond2,3);

end
cond1 = grp_MTL_cond1;
cond2 = grp_MTL_cond2;


MTL_grp_1_mn = nanmean(grp_MTL_cond1,3);
MTL_grp_1_std = nanstd(grp_MTL_cond1,0,3);

MTL_grp_2_mn = nanmean(grp_MTL_cond2,3);
MTL_grp_2_std = nanstd(grp_MTL_cond2,0,3);


save('MTL_grp_spectrograms','cond1','cond2')


%% CA3
 cond1 = [];
 cond2 = [];

 
 grp = {'39' '57' '44' '66'}

for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
      grp_CA3_cond1(:,:,a) = nanmean(CA3_cond1,3);
      grp_CA3_cond2(:,:,a) = nanmean(CA3_cond2,3);

end
cond1 = grp_CA3_cond1;
cond2 = grp_CA3_cond2;


CA3_grp_1_mn = nanmean(grp_CA3_cond1,3);
CA3_grp_1_std = nanstd(grp_CA3_cond1,0,3);

CA3_grp_2_mn = nanmean(grp_CA3_cond2,3);
CA3_grp_2_std = nanstd(grp_CA3_cond2,0,3);



save('CA3_grp_spectrograms','cond1','cond2')


%% CA1
 cond1 = [];
 cond2 = [];

  grp = {'39' '57' '44' '66' '84'}
for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
      grp_CA1_cond1(:,:,a) = nanmean(CA1_cond1,3);
      grp_CA1_cond2(:,:,a) = nanmean(CA1_cond2,3);

end
cond1 = grp_CA1_cond1;
cond2 = grp_CA1_cond2;


CA1_grp_1_mn = nanmean(grp_CA1_cond1,3);
CA1_grp_1_std = nanstd(grp_CA1_cond1,0,3);

CA1_grp_2_mn = nanmean(grp_CA1_cond2,3);
CA1_grp_2_std = nanstd(grp_CA1_cond2,0,3);



save('CA1_grp_spectrograms','cond1','cond2')
%% ca1-3
 cond1 = [];
 cond2 = [];

 
CA13_1 = nan(size(grp_CA1_cond1));
CA13_2 = nan(size(grp_CA1_cond2));


for a = 1:length(grp)
    CA13_1(:,:,a) = nanmean(cat(3,grp_CA1_cond1(:,:,a), grp_CA3_cond1(:,:,a)),3);
    CA13_2(:,:,a) = nanmean(cat(3,grp_CA1_cond2(:,:,a), grp_CA3_cond2(:,:,a)),3);


end
cond1 = CA13_1;
cond2 = CA13_2;


CA13_grp_1_mn = nanmean(CA13_1,3);
CA13_grp_1_std = nanstd(CA13_1,0,3);

CA13_grp_2_mn = nanmean(CA13_2,3);
CA13_grp_2_std = nanstd(CA13_2,0,3);




save('CA13_grp_spectrograms','cond1','cond2')

%% OFC
 cond1 = [];
 cond2 = [];

for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
  
      grp_OFC_cond1(:,:,a) = nanmean(OFC_cond1,3);
      grp_OFC_cond2(:,:,a) = nanmean(OFC_cond2,3);


end
cond1 = grp_OFC_cond1;
cond2 = grp_OFC_cond2;


OFC_grp_1_mn  = nanmean(grp_OFC_cond1,3);
OFC_grp_1_std = nanstd(grp_OFC_cond1,0,3);

OFC_grp_2_mn  = nanmean(grp_OFC_cond2,3);
OFC_grp_2_std = nanstd(grp_OFC_cond2,0,3);



save('OFC_grp_spectrograms','cond1','cond2')

%% CING
 cond1 = [];
 cond2 = [];

 grp = {'39' '57' '44' '63' '66' }

for a = 1:length(grp)
    load(['subj' grp{a} 'spectrograms.mat'])
      grp_CING_cond1(:,:,a) = nanmean(cing_cond1,3);
      grp_CING_cond2(:,:,a) = nanmean(cing_cond2,3);

end
cond1 = grp_CING_cond1;
cond2 = grp_CING_cond2;


CING_grp_1_mn = nanmean(grp_CING_cond1,3);
CING_grp_1_std = nanstd(grp_CING_cond1,0,3);

CING_grp_2_mn = nanmean(grp_CING_cond2,3);
CING_grp_2_std = nanstd(grp_CING_cond2,0,3);



save('CING_grp_spectrograms','cond1','cond2')
%%
cond_num = 2
edge_points = 200
mx = .28
mn = -.28

figure
suptitle('OFC - 6 patients')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), OFC_grp_1_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), OFC_grp_2_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

%% fro


mx = .27
mn = -.27

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
title('lure --> +')
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
title('lure --> -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

saveas(gcf,'temp.jpg')


%% temporal
mx = .25
mn = -.25
figure
suptitle('temporal')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(temp_grp_1_mn,2)), 1:length(freq), temp_grp_1_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), temp_grp_2_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
saveas(gcf,'temp.jpg')


%% CING

figure
suptitle('CING')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CING_grp_1_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CING_grp_2_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
saveas(gcf,'temp.jpg')







%% MTL


mx = .2
mn = -.2

figure
suptitle('MTL')
subplot (cond_num,1,1)
%imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_1_mn./MTL_grp_1_std)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_1_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
%imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_2_mn./MTL_grp_2_std)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(MTL_grp_1_mn,2)), 1:length(freq), MTL_grp_2_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
saveas(gcf,'temp.jpg')

%%
mx = .23
mn = -.23
figure
suptitle('CA3 - 4 patients')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq), CA3_grp_1_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(fro_grp_1_mn,2)), 1:length(freq),  CA3_grp_2_mn)
set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
saveas(gcf,'temp.jpg')

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
title('lure -> +')
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
title('lure -> -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
saveas(gcf,'temp.jpg')


%% %% cortex
load('TEMP_grp_spectrograms.mat')
temp1 = nanmean(cond1,3);
temp2 = nanmean(cond2,3);
load('OFC_grp_spectrograms.mat')
ofc2 = nanmean(cond2,3);
ofc1 = nanmean(cond1,3);
load('FRO_grp_spectrograms.mat')
fro1 = nanmean(cond1,3);
fro2 = nanmean(cond2,3);
load('CING_grp_spectrograms.mat')
cing2 = nanmean(cond2,3);
cing1 = nanmean(cond1,3);
cortex1 = cat(3,fro1,cing1,ofc1,temp1);
cortex2 = cat(3,fro2,cing2,ofc2,temp2);


%%
mx = 3
mn = -2.5

figure
suptitle('cortex')
subplot (cond_num,1,1)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(cortex1,2)), 1:length(freq), nanmean(cortex1,3))
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(cortex1,2)), 1:length(freq), nanmean(cortex1,3)./nanstd(cortex1,0,3))

set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar 
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> +')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')

subplot (cond_num,1,2)
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(cortex1,2)), 1:length(freq),  nanmean(cortex2,3))
imagesc(linspace(-pre_stim+(edge_points/fs),post_stim-(edge_points/fs),size(cortex1,2)), 1:length(freq), nanmean(cortex2,3)./nanstd(cortex2,0,3))

set(gca, 'YTick',tickmarks,'YTickLabel',round(freq(tickmarks)))
colorbar
colormap jet
caxis([mn mx])
xlabel('time')
ylabel('freq')
hold on
plot(zeros(1,length(freq)), 1:length(freq),'Color','k','LineWidth',1)
title('lure -> -')
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
saveas(gcf,'temp.jpg')


%% generate cond2, cond 3 for permutation testing on cortex
regions_list = {'OFC' 'FRO' 'TEMP' 'CING'}
subj_list    = {'39' '57' '44' '63' '66' '84'}
total_subj   = length(subj_list)
cond2a = cell(1,4*6)
cond3a = cell(1,4*6)

% ofc
for subj = 1:total_subj
    load(['subj' subj_list{subj} 'spectrograms.mat'])
    cond2a{subj} = OFC_cond2;
    cond3a{subj} = OFC_cond3;
end

% fro
reg=2
for subj = 1:total_subj
    load(['subj' subj_list{subj} 'spectrograms.mat'])
    cond2a{subj+(total_subj*(reg-1))} = fro_cond2;
    cond3a{subj+(total_subj*(reg-1))} = fro_cond3;
end

% temp
reg=3
for subj = 1:total_subj
    load(['subj' subj_list{subj} 'spectrograms.mat'])
    cond2a{subj+(total_subj*(reg-1))} = temp_cond2;
    cond3a{subj+(total_subj*(reg-1))} = temp_cond3;
end

% cing
subj_list    = {'39' '57' '44' '63' '66' }
total_subj   = length(subj_list)

reg=4
for subj = 1:total_subj
    load(['subj' subj_list{subj} 'spectrograms.mat'])
    cond2a{subj+(total_subj*(reg-1))} = cing_cond2;
    cond3a{subj+(total_subj*(reg-1))} = cing_cond3;
end
cond2 = cat(3,cond2a{:});
cond3 = cat(3,cond3a{:});







