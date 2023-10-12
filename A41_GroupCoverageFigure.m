clear all;close all;clc
addpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-20181130'))

load('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-20181130/template/anatomy/surface_pial_left.mat');
left_mesh = mesh; clear mesh
load('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-20181130/template/anatomy/surface_pial_right.mat');
righ_mesh = mesh;
%% lateral view: right hemi
figure;hold on
ft_plot_mesh(left_mesh, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.4,'EdgeColor', 'none');
ft_plot_mesh(righ_mesh, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.4,'EdgeColor', 'none');

subj_list = {'39' '44'  '57'  '63' '66' '84'}
size_vector=15
for iSubj = 1:length(subj_list)             
subj=subj_list{iSubj}
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/FT_Pipeline/Electrodes'])
load(['IR' subj '_elec_mni_v.mat'])
hs = ft_plot_sens(elec_mni_v, 'elecsize', size_vector, 'facecolor', 'b', 'shape', 'sphere');
end

view([90 0]); % left to right, vs. up and down
material dull;
lighting gouraud;
camlight;
print('-clipboard','-dbitmap')
%% lateral view: left hemi
figure;hold on
ft_plot_mesh(left_mesh, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.4,'EdgeColor', 'none');
ft_plot_mesh(righ_mesh, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.4,'EdgeColor', 'none');

subj_list = {'39' '44'  '57'  '63' '66' '84'}
size_vector=15
for iSubj = 1:length(subj_list)             
subj=subj_list{iSubj}
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/FT_Pipeline/Electrodes'])
load(['IR' subj '_elec_mni_v.mat'])
hs = ft_plot_sens(elec_mni_v, 'elecsize', size_vector, 'facecolor', 'b', 'shape', 'sphere');
end
view([-90 0]); % left to right, vs. up and down
material dull;
lighting gouraud;
camlight;
print('-clipboard','-dbitmap')

%% anterior view
figure;hold on
ft_plot_mesh(left_mesh, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.4,'EdgeColor', 'none');
ft_plot_mesh(righ_mesh, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.4,'EdgeColor', 'none');

subj_list = {'39' '44'  '57'  '63' '66' '84'}
size_vector=15
for iSubj = 1:length(subj_list)             
subj=subj_list{iSubj}
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/FT_Pipeline/Electrodes'])
load(['IR' subj '_elec_mni_v.mat'])
hs = ft_plot_sens(elec_mni_v, 'elecsize', size_vector, 'facecolor', 'b', 'shape', 'sphere');
end
view([-180 0]);
material dull;
lighting gouraud;
camlight;
print('-clipboard','-dbitmap')
%% ventral view
figure;hold on
ft_plot_mesh(left_mesh, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.4,'EdgeColor', 'none');
ft_plot_mesh(righ_mesh, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.4,'EdgeColor', 'none');

subj_list = {'39' '44'  '57'  '63' '66' '84'}
size_vector=15
for iSubj = 1:length(subj_list)             
subj=subj_list{iSubj}
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/FT_Pipeline/Electrodes'])
load(['IR' subj '_elec_mni_v.mat'])
hs = ft_plot_sens(elec_mni_v, 'elecsize', size_vector, 'facecolor', 'b', 'shape', 'sphere');
end
view([-180 -70]);
material dull;
lighting gouraud;
camlight;
print('-clipboard','-dbitmap')

