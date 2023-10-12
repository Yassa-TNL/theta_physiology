% run A003, then A011, then this

subj = '57'
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj])
if strcmp(subj, '55')
    load('elec_loc_wrkspc.mat')
end
addpath(genpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final/fieldtrip-20181130'))
load(['/mnt/yassamri/iEEG/sandra/subj_' subj '/FT_Pipeline/Electrodes/IR'  subj '_elec_acpc_f.mat'])
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/FT_Pipeline/Scans'])
if subj == '66'
    mri= ft_read_mri(['IR' subj '_fsMR_post_acpc.nii'])
else
    mri = ft_read_mri(['IR' subj '_fsMR_pre_acpc.nii']);
end
cd(['/mnt/yassamri/iEEG/sandra/subj_' subj '/FT_Pipeline/Surfaces'])
cortex1 = ft_read_headshape(['IR' subj '_cortex_lh.mat']);
cortex2 = ft_read_headshape(['IR' subj '_cortex_rh.mat']);
cfg = [];
cfg.elec = elec_acpc_f;
%ft_electrodeplacement(cfg, mri)

%% get desired chan index
original_label = cfg.elec.label
chan_label_new = cell(length(chan_label),1)
if  ~strcmp('84',subj)
    for chan = 1:length(chan_label)
        
        chan_label_new{chan} = chan_label{chan}(4:end-3)
    end
end
if  strcmp('63',subj)
    chan_label_new{117} = 'RLES1'
elseif strcmp('84',subj)
    chan_label_new = chan_label;
end


if strcmp('83',subj) | strcmp('63',subj)
    % get a vector of channels inside brain
    for chan = 1:length(original_label)
        temp = original_label{chan}
        for chan_or = 1:chan_counter
            if strcmp(temp, chan_label_new{chan_or})
                desired_chan_idx (chan) = chan_or
            end
        end
    end
    
else
    
    
    for chan = 1:chan_counter
        temp = chan_label_new{chan}
        for chan_or = 1:length(original_label)
            if strcmp(temp, original_label{chan_or})
                desired_chan_idx (chan) = chan_or;
            end
        end
    end
end

% indicate the size
size_vector = zeros(length(original_label),1);
size_vector(:) = 20;

if strcmp('83',subj) | strcmp('63',subj)
    size_vector(:) = (chan_peak_freq_disc_magnitude(desired_chan_idx,2))*150;
    peak_freq = chan_peak_freq_disc_magnitude(desired_chan_idx,1);
    size_vector(find(size_vector==0)) = 20;

else
    size_vector(desired_chan_idx) = (chan_peak_freq_disc_magnitude(:,2))*150;
    peak_freq = chan_peak_freq_disc_magnitude(:,1);

    % add 20 where for non disc chans
    size_vector(find(size_vector==0)) = 20;
end

% color vector for peak freq, black with white outline for zero discrinimination
my_colormap = jet(length(all_poss_unique_freq));
unique_vals = sort(unique(all_poss_unique_freq));

freq_color_map = zeros(length(peak_freq),3);
white = ones(1,3)

for chan = 1:length(peak_freq) % loop thru chans
    % if a peak discriminatory freq exists
    if peak_freq(chan)==0
        freq_color_map(chan,:) = white(1,:);
    else
        idx = find(peak_freq(chan)==unique_vals);
        freq_color_map(chan,:) = my_colormap(idx,:);
    end
end


freq_color_map_all = zeros(length(original_label),3);
if strcmp('83',subj) | strcmp('63',subj)

        freq_color_map_all = freq_color_map;

else
    for a = 1:length(desired_chan_idx)
        freq_color_map_all(desired_chan_idx(a),:) = freq_color_map(a,:);
    end
end


%% plot
figure
ft_plot_mesh(cortex1, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.5,'EdgeColor', 'none');
ft_plot_mesh(cortex2, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.5,'EdgeColor', 'none');
view([-90 10]); lighting gouraud; camlight;
%view([90 -10]); lighting gouraud; camlight;

hs = ft_plot_sens(elec_acpc_f, 'elecsize', size_vector, 'facecolor', freq_color_map_all, 'shape', 'sphere');
hs = ft_plot_sens(elec_acpc_f, 'elecsize', 1.5*size_vector, 'facecolor', 'b', 'shape', 'sphere');

colormap(my_colormap)
unique_vals = round(unique_vals,2)
colorbar('Ticks',linspace(0,1,length(unique_vals)), 'TickLabels', {unique_vals})
%title(['IR' subj  ' ' num2str(strt_time) '-' num2str(end_time) ' ms'])
set(gca, 'FontSize', 14, 'FontWeight', 'bold')
saveas(gcf, 'ventralbrain.bmp') 

%% to identify elec loc
freq_color_map_all(1:64,:) = 0; 
freq_color_map_all(65:68,:) =0; 

size_vector(:) = 30


%%
% colormap(my_colormap)
% colorbar
% figure
% imagesc(peak_freq)
% colormap(my_colormap)
% colorbar
% figure;
% for a = 1:2%size(my_colormap,1)
% hold on
% plot(1:10, [a:a+9], 'color', my_colormap(a,:), 'LineWidth', 10)
% end


% color vector for peak freq
my_colormap = jet(3*length(all_poss_unique_freq));
unique_vals = sort(unique(all_poss_unique_freq));

freq_color_map = zeros(length(peak_freq),3);
for chan = 1:length(peak_freq) % loop thru chans
    idx = find(peak_freq(chan)==unique_vals);
    freq_color_map(chan,:) = my_colormap(idx+3,:);
end
for a = 1:length(desired_chan_idx)
    freq_color_map_all(desired_chan_idx(a),:) = freq_color_map(a,:);
end

%% plot coverage only
figure
ft_plot_mesh(cortex1, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.5,'EdgeColor', 'none');
ft_plot_mesh(cortex2, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.5,'EdgeColor', 'none');
view([-90 10]); lighting gouraud; camlight;
hs = ft_plot_sens(elec_acpc_f, 'elecsize', 1.5*size_vector, 'facecolor', 'b', 'shape', 'sphere');

figure
ft_plot_mesh(cortex1, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.5,'EdgeColor', 'none');
ft_plot_mesh(cortex2, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.5,'EdgeColor', 'none');
view([90 10]); lighting gouraud; camlight;
hs = ft_plot_sens(elec_acpc_f, 'elecsize', 1.5*size_vector, 'facecolor', 'b', 'shape', 'sphere');

figure
ft_plot_mesh(cortex1, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.5,'EdgeColor', 'none');
ft_plot_mesh(cortex2, 'facecolor', [0.781 0.762 0.664], 'facealpha', 0.5,'EdgeColor', 'none');
view([180 -90]); lighting gouraud; camlight;
hs = ft_plot_sens(elec_acpc_f, 'elecsize', 1.5*size_vector, 'facecolor', 'b', 'shape', 'sphere');

print('-clipboard','-dbitmap')