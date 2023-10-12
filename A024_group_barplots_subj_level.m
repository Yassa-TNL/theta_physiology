% plot bar plot

stim_onset = 301;


%% CA1-3-dg
rng = 4
dur = 750
subj_length = 5
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(CA13_1(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CA13_2(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CA13_3(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CA13_4(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(CA13_1(freq<5,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CA13_2(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CA13_3(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CA13_4(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)];
figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/sqrt(subj_length), 'rx')%ylim([mn mx])
ylabel('3-4 hz power')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
title('Hippocampus')
saveas(gcf, 'delta_barplot_HC.bmp') 
ylim([-.1 .25])
saveas(gcf, 'delta_barplot_HC.bmp') 

%% CA3-dg
rng = 4
dur = 750
subj_length = 4
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(CA3_grp_1_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CA3_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CA3_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CA3_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(CA3_grp_1_mn(freq<5,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CA3_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CA3_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CA3_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)];
figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/sqrt(subj_length), 'rx')%ylim([mn mx])
ylabel('3-4 hz power')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
title('CA3 - DG')
ylim([-.1 .25])
saveas(gcf, 'delta_barplot_CA3.bmp') 
%% CA1-dg
rng = 4
dur = 750
subj_length = 4
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(CA1_grp_1_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CA1_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CA1_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CA1_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(CA1_grp_1_mn(freq<5,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CA1_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CA1_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CA1_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)];
figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/sqrt(subj_length), 'rx')%ylim([mn mx])
ylabel('3-4 hz power')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
title('CA1')
ylim([-.1 .25])
saveas(gcf, 'delta_barplot_CA1.bmp') 

%% fro
rng = 4
dur = 1000
subj_length = 7
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(fro_grp_1_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(fro_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(fro_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(fro_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(fro_grp_1_mn(freq<5,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(fro_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(fro_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(fro_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)];
figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/sqrt(subj_length), 'rx')%ylim([mn mx])
ylabel('3-4 hz power')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
title('frontal: 0-500 ms')
ylim([-.1 .25])
saveas(gcf, 'delta_barplot_FRO.bmp') 

%% temporal
stim_onset = 301+800
dur = 300
% stim_onset = 301
% dur = 500
rng = 4
subj_length = 7
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(temp_grp_1_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(temp_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(temp_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(temp_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(temp_grp_1_mn(freq<5,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(temp_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(temp_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(temp_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)];
figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/sqrt(subj_length), 'rx')%ylim([mn mx])
ylabel('3-4 hz power')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
title('Temp: 800-110 ms')
%title('Temp: 0-500 ms')

ylim([-.1 .25])
saveas(gcf, 'delta_barplot_temp.bmp') 

%% MTL

stim_onset = 301
dur = 750
rng = 4
subj_length = 6
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(MTL_grp_1_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(MTL_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(MTL_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(MTL_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(MTL_grp_1_mn(freq<5,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(MTL_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(MTL_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(MTL_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)];
figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/sqrt(subj_length), 'rx')%ylim([mn mx])
ylabel('3-4 hz power')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
title('MTL: 0-750 ms')
ylim([-.1 .25])
saveas(gcf, 'delta_barplot_MTL.bmp') 


%% OFC
stim_onset = 301
dur = 400
rng = 4
subj_length = 6
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(OFC_grp_1_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(OFC_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(OFC_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(OFC_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(OFC_grp_1_mn(freq<5,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(OFC_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(OFC_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(OFC_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)];
figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/sqrt(subj_length), 'rx')%ylim([mn mx])
ylabel('3-4 hz power')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
title('OFC: 0-400 ms')
ylim([-.1 .25])
saveas(gcf, 'delta_barplot_OFC.bmp') 


%% CING
stim_onset = 301
dur = 750
rng = 4
subj_length = 5
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(CING_grp_1_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CING_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CING_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(CING_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(OFC_grp_1_mn(freq<5,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CING_grp_2_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CING_grp_3_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(CING_grp_4_mn(freq<rng,stim_onset:stim_onset+dur,: ),1))', 2),0,1)];
figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/sqrt(subj_length), 'rx')%ylim([mn mx])
ylabel('3-4 hz power')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
title('OFC: 0-750 ms')
ylim([-.1 .25])
saveas(gcf, 'delta_barplot_CING.bmp')


%% fro high theta: t-value
rng         = [6 9]
dur1        = 200
dur2        = 470

cond1 = nanmean(nanmean(temp_grp_1_mn(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: )./...
            temp_grp_1_std(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,:)));
        
        
cond2 = nanmean(nanmean(temp_grp_2_mn(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: )./...
            temp_grp_2_std(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,:)));
        
        
cond3 = nanmean(nanmean(temp_grp_3_mn(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: )./...
            temp_grp_3_std(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,:)));
        
        
cond4 = nanmean(nanmean(temp_grp_4_mn(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: )./...
            temp_grp_4_std(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,:)));
bar_vector_mn1 = [cond1 cond2 cond3 cond4];

figure
hold on
bar(bar_vector_mn1)   
ylabel(['reliable increase in 6-9 hz power'])
title('Temporal Cortex')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')

saveas(gcf, 'temp_hightheta_barplot_HC.bmp') 

%% fro high theta: t-value
rng         = [6 9]
subj_length = 5
dur1        = 200
dur2        = 470

cond1 = nanmean(nanmean(fro_grp_1_mn(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: )./...
            fro_grp_1_std(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,:)));
        
        
cond2 = nanmean(nanmean(fro_grp_2_mn(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: )./...
            fro_grp_2_std(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,:)));
        
        
cond3 = nanmean(nanmean(fro_grp_3_mn(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: )./...
            fro_grp_3_std(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,:)));
        
        
cond4 = nanmean(nanmean(fro_grp_4_mn(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: )./...
            fro_grp_4_std(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,:)));
bar_vector_mn1 = [cond1 cond2 cond3 cond4];

figure
hold on
bar(bar_vector_mn1)         

%% MTL Gamma reponse locked
rng = 33
stim_onset = 1301
dur1=950
dur2=420
subj_length = 5
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(grp_MTL_cond1(freq>rng,stim_onset-dur1:stim_onset-dur2,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(grp_MTL_cond2(freq>rng,stim_onset-dur1:stim_onset-dur2,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(grp_MTL_cond3(freq>rng,stim_onset-dur1:stim_onset-dur2,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(grp_MTL_cond4(freq>rng,stim_onset-dur1:stim_onset-dur2,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(grp_MTL_cond1(freq>rng,stim_onset-dur1:stim_onset-dur2,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(grp_MTL_cond2(freq>rng,stim_onset-dur1:stim_onset-dur2,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(grp_MTL_cond3(freq>rng,stim_onset-dur1:stim_onset-dur2,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(grp_MTL_cond4(freq>rng,stim_onset-dur1:stim_onset-dur2,: ),1))', 2),0,1)];
%% 
rng = [33 84]

bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(grp_MTL_cond1(freq>rng(1) & freq<rng(2),stim_onset-dur1:stim_onset-dur2,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(grp_MTL_cond2(freq>rng(1) & freq<rng(2),stim_onset-dur1:stim_onset-dur2,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(grp_MTL_cond3(freq>rng(1) & freq<rng(2),stim_onset-dur1:stim_onset-dur2,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(grp_MTL_cond4(freq>rng(1) & freq<rng(2),stim_onset-dur1:stim_onset-dur2,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(grp_MTL_cond1(freq>rng(1) & freq<rng(2),stim_onset-dur1:stim_onset-dur2,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(grp_MTL_cond2(freq>rng(1) & freq<rng(2),stim_onset-dur1:stim_onset-dur2,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(grp_MTL_cond3(freq>rng(1) & freq<rng(2),stim_onset-dur1:stim_onset-dur2,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(grp_MTL_cond4(freq>rng(1) & freq<rng(2),stim_onset-dur1:stim_onset-dur2,: ),1))', 2),0,1)];


figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/sqrt(subj_length), 'rx')%ylim([mn mx])
ylabel('33 - 84 hz power')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
title('MTL')

saveas(gcf, 'gamma_barplot_MTL.bmp') 




%% fro high theta
% rng = [6 10]
% subj_length = 5
% dur1 = 100
% dur2 = 450
% bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(grp_fro_cond1(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2))...
%     nanmean(nanmean(squeeze(nanmean(grp_fro_cond2(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2))...
%     nanmean(nanmean(squeeze(nanmean(grp_fro_cond3(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2))...
%     nanmean(nanmean(squeeze(nanmean(grp_fro_cond4(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2))];
%          
% bar_vector_std = [nanstd(nanmean(squeeze(nanmean(grp_fro_cond1(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2),0,1)...
%     nanstd(nanmean(squeeze(nanmean(grp_fro_cond2(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2),0,1)...
%     nanstd(nanmean(squeeze(nanmean(grp_fro_cond3(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2),0,1)...
%     nanstd(nanmean(squeeze(nanmean(grp_fro_cond4(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2),0,1)];
% figure
% hold on
% bar(bar_vector_mn1)         
% errorbar(1:4,bar_vector_mn1,bar_vector_std/subj_length, 'rx')%ylim([mn mx])
% ylabel('average 6-10 hz power')
% title('Frontal Cortex')
% set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
% set(gca, 'FontSize', 16, 'FontWeight', 'bold')
% saveas(gcf, 'fro_hightehta_barplot_HC.bmp') 
% 
% %%
rng = [6 9]
subj_length = 5
dur1 = 200
dur2 = 470
bar_vector_mn1 = [nanmean(nanmean(squeeze(nanmean(grp_temp_cond1(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(grp_temp_cond2(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(grp_temp_cond3(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2))...
    nanmean(nanmean(squeeze(nanmean(grp_temp_cond4(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2))];
         
bar_vector_std = [nanstd(nanmean(squeeze(nanmean(grp_temp_cond1(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(grp_temp_cond2(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(grp_temp_cond3(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2),0,1)...
    nanstd(nanmean(squeeze(nanmean(grp_temp_cond4(freq>rng(1) & freq<rng(2),stim_onset+dur1:stim_onset+dur2,: ),1))', 2),0,1)];
figure
hold on
bar(bar_vector_mn1)         
errorbar(1:4,bar_vector_mn1,bar_vector_std/subj_length, 'rx')%ylim([mn mx])
ylabel('average 6-10 hz power')
title('Temporal Cortex')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
set(gca, 'FontSize', 16, 'FontWeight', 'bold')
saveas(gcf, 'temp_hightehta_barplot_HC.bmp') 