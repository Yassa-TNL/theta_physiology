% plot cluster determined group bar plots
clear all;close all;clc
reg           = 'NC'
% NC HC CA3 'OFC_FRO_TEMP' 'OFC_FRO_TEMP_significant' 'OFC_FRO_TEMP_CING_INS_EC_significant'
fn_ext        = '_cue_responseiveyes'
baseline      = [ '_cond_spec_prestim' fn_ext]
lock          = 'response' % onset % response
exp_type      = 'tuning_correct'% 'encoding' 'tuning_correct'
ref           = 'LM' 
plot_4_conds  = 'yes'

addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock ])

subj_list = {'39' '44' '57' '63' '66' '84' '85' '87'}  
fs=500;

if strcmp('onset', lock)
    pre_stim = 0.5; post_stim = 2.0;
    if strcmp('encoding',exp_type)
        time_idx = [.3*fs 2.5*fs]
    else
        time_idx = [.3*fs 1.5*fs]
    end
elseif strcmp('response', lock)
    time_idx = [0.1*fs 1.2*fs];
    pre_stim = 0.5; post_stim = 2.0;

elseif strcmp('onset2', lock)
    pre_stim = 0.5; post_stim = 2.6;
end

% get pooled elecs from all subj
[cond1a,cond2a,cond3a,cond4a] = get_power_subj_elecs(subj_list,exp_type,reg,ref,baseline,lock);

if strcmp('tuning_correct',exp_type)
    cond1 = cat(3,cond1a{:});
    cond4 = cat(3,cond4a{:});
end

cond2 = cat(3,cond2a{:});
cond3 = cat(3,cond3a{:});
if ~strcmp('response', lock)
reset = nanmean(nanmean(cond2(:,0.3*fs:pre_stim*fs,:),3),2);
cond2 = cond2 - reset;
reset = nanmean(nanmean(cond3(:,0.3*fs:pre_stim*fs,:),3),2);
cond3 = cond3 - reset;
reset = nanmean(nanmean(cond1(:,0.3*fs:pre_stim*fs,:),3),2);
cond1 = cond1 - reset;
reset = nanmean(nanmean(cond4(:,0.3*fs:pre_stim*fs,:),3),2);
cond4 = cond4 - reset;
end

% load cluster matrix and preview
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/cluster_matrices'])
range = 'deltatheta'   %theta delta deltatheta alpha gamma
if strcmp('OFC_FRO_TEMP',reg) ||  strcmp('OFC_FRO_TEMP_significant',reg) ...
        || strcmp('OFC_FRO_TEMP_CING_INS_EC_significant',reg)
    regprime = 'NC';
    load([regprime '_' lock '_' exp_type '_groupcluster_' range '.mat'])
else
 load([reg '_' lock '_' exp_type '_groupcluster_' range '.mat'])   
end
figure;imagesc(zmapthresh)
zmapthresh_temp = zeros(size(cond2,1),size(cond2,2));
zmapthresh_temp(:,time_idx(1):time_idx(2)) = logical(zmapthresh);


if strcmp('gamma',range)
    minfreq       = 3;
    maxfreq       = 200;
    [freq] = get_freq(fs,minfreq, maxfreq);
    zmapthresh_temp(:)=0;
    zmapthresh_temp(freq>35,1.5*fs:2*fs)=1;
    figure;imagesc(zmapthresh_temp)
end
%% group plot
% make group mean bar plots
for elec = 1:size(cond2,3)
    
    if strcmp('encoding', exp_type)
        temp2 = [];
        temp2 = (cond2(:,:,elec));
        cond2_grp_val(elec) = [nanmean(temp2(logical(zmapthresh_temp)) )];
        
        temp3 = [];
        temp3 = (cond3(:,:,elec));
        cond3_grp_val(elec) = [nanmean(temp3(logical(zmapthresh_temp)) )];
    else strcmp('tuning_correct', exp_type)
        
        temp1 = [];
        temp1 = (cond1(:,:,elec));
        cond1_grp_val(elec) = [nanmean(temp1(logical(zmapthresh_temp)) )];
        
        temp2 = [];
        temp2 = (cond2(:,:,elec));
        cond2_grp_val(elec) = [nanmean(temp2(logical(zmapthresh_temp)) )];
        
        temp3 = [];
        temp3 = (cond3(:,:,elec));
        cond3_grp_val(elec) = [nanmean(temp3(logical(zmapthresh_temp)) )];
        
        temp4 = [];
        temp4 = (cond4(:,:,elec));
        cond4_grp_val(elec) = [nanmean(temp4(logical(zmapthresh_temp)) )];
    end
end
if strcmp('encoding', exp_type)
    bar_vector_mn1  = [nanmean(cond3_grp_val)  nanmean(cond2_grp_val) ];
    bar_vector_std1 = [nanstd(cond3_grp_val,0,2) nanstd(cond2_grp_val,0,2) ];
    labels =  {'lure +' 'lure -' };
    
    % encoding stats
    reps = 1000;
    adata = cond3_grp_val;
    bdata = cond2_grp_val;
    p2    = permutation_paired(adata, bdata, reps)
    
else
    bar_vector_mn1 = [nanmean(cond1_grp_val) nanmean(cond2_grp_val) nanmean(cond3_grp_val) nanmean(cond4_grp_val)];
    bar_vector_std1 = [nanstd(cond1_grp_val,0,2) nanstd(cond2_grp_val,0,2) nanstd(cond3_grp_val,0,2) nanstd(cond4_grp_val,0,2)];
    labels = {'repeat' 'lure -' 'lure +' 'new'};

    % retrieval stats
reps  = 1000;
adata = cond3_grp_val; 
bdata = cond1_grp_val;
[p1,real_meandiff(1) ]    = permutation_paired(adata, bdata, reps)

adata = cond3_grp_val; 
bdata = cond2_grp_val;
[p2,real_meandiff(2) ]    = permutation_paired(adata, bdata, reps)

adata = cond3_grp_val; 
bdata = cond4_grp_val;
[p3,real_meandiff(3) ]    = permutation_paired(adata, bdata, reps)

alpha = .05
pvals = [p1 p3];
[p_fdr, p_masked] = fdr( pvals, alpha)
[pvals p_fdr]

end

% plot
figure
hold on
if strcmp('yes',plot_4_conds) && strcmp('tuning_correct',exp_type) || strcmp('encoding',exp_type)
bar([bar_vector_mn1])         
errorbar(1:length(bar_vector_mn1),bar_vector_mn1,bar_vector_std1/sqrt(size(cond2,3)), 'rx')
set(gca, 'XTick', 1:length(bar_vector_mn1), 'XTickLabel', labels,'XTickLabelRotation',45)
else
bar([bar_vector_mn1(2:3)])         
errorbar(1:length(bar_vector_mn1(2:3)),bar_vector_mn1(2:3),bar_vector_std1(2:3)/sqrt(size(cond2,3)), 'rx')
set(gca, 'XTick', 1:length(bar_vector_mn1(2:3)), 'XTickLabel', labels(2:3),'XTickLabelRotation',45)
end
%ylim([-.1 .5])
y=ylim;
if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
end
print('-clipboard','-dbitmap')
cd('/mnt/yassamri/iEEG/sandra/GroupFigures')
set(gca, 'FontSize', 16, 'FontWeight', 'bold') 
%title([reg ' ' range ])
%ylabel([reg ' power'])


%% Compute and overlay plot indiv sub vals - a single val per trial across chans
clear cond1 cond2 cond3 cond4 cond2_grp_val cond3_grp_val
cd(['/mnt/yassamri/iEEG/sandra/group_data/groupdata_spectrograms/' ref '_reref/' exp_type '_' lock  ])
cond1_grp_mn= nan(length(subj_list),1);
cond1_grp_sem=nan(length(subj_list),1);

cond2_grp_mn= nan(length(subj_list),1);
cond2_grp_sem= nan(length(subj_list),1);

cond3_grp_mn= nan(length(subj_list),1);
cond3_grp_sem= nan(length(subj_list),1);

cond4_grp_mn= nan(length(subj_list),1);
cond4_grp_sem= nan(length(subj_list),1);
for iSubj = 1:length(subj_list)

    %  use for retrieval
    load(['subj' subj_list{iSubj}  'spectrograms_all_trials_across_channels_cond_spec_prestim_cue_responseiveyes'])

    if  strcmp('OFC',reg)
        data1 = OFC_cond1;
        data2 = OFC_cond2;
        data3 = OFC_cond3;
        data4 = OFC_cond4;
        
    elseif strcmp('FRO',reg)
        data1 = fro_cond1;
        data2 = fro_cond2;
        data3 = fro_cond3;
        data4 = fro_cond4;
        
    elseif strcmp('TEMP',reg)
        data1 = temp_cond1;
        data2 = temp_cond2;
        data3 = temp_cond3;
        data4 = temp_cond4;  
        
    elseif strcmp('CING',reg)
        data1 = cing_cond1;
        data2 = cing_cond2;
        data3 = cing_cond3;
        data4 = cing_cond4; 
        
    elseif strcmp('ins',reg)
        data1 = ins_cond1;
        data2 = ins_cond2;
        data3 = ins_cond3;
        data4 = ins_cond4;  
        
    elseif strcmp('EC',reg)
        data1 = EC_cond1;
        data2 = EC_cond2;
        data3 = EC_cond3;
        data4 = EC_cond4;   
        
    elseif strcmp('HC',reg)
        data1 = HC_cond1;
        data2 = HC_cond2;
        data3 = HC_cond3;
        data4 = HC_cond4;
        
    elseif strcmp('CA1',reg)
        data1 = CA1_cond1;
        data2 = CA1_cond2;
        data3 = CA1_cond3;
        data4 = CA1_cond4;
        
    elseif strcmp('CA3',reg)
        data1 = CA3_cond1;
        data2 = CA3_cond2;
        data3 = CA3_cond3;
        data4 = CA3_cond4;
        
    elseif strcmp('NC',reg)
        data1 = NC_cond1;
        data2 = NC_cond2;
        data3 = NC_cond3;
        data4 = NC_cond4;
    end
      if sum(sum(sum(isnan(data2))))== length(data2(:)) % if subj doesnt have coverage 
      continue
      end
      
      if ~strcmp('response', lock)
            % reset to zero
      reset1 = nanmean(nanmean(data1(:,0.3*fs:pre_stim*fs,:),3),2);
      data1 = data1-reset1;
      reset2 = nanmean(nanmean(data2(:,0.3*fs:pre_stim*fs,:),3),2);
      data2 = data2-reset2;
      reset3 = nanmean(nanmean(data3(:,0.3*fs:pre_stim*fs,:),3),2);
      data3 = data3-reset3;
      reset4 = nanmean(nanmean(data4(:,0.3*fs:pre_stim*fs,:),3),2);
      data4 = data4-reset4;
      end
      
  
    if strcmp('encoding', exp_type)
        cond1 = []; cond2 = []; cond3 = []; cond4 = [];
        for iTrl = 1:size(data1,3)
            temp1 = data1(:,:,iTrl);
            temp2 = nanmean(temp1(logical(zmapthresh)),1);
            if ~isnan(temp2)
                cond1 = [cond1 temp2];
            end
            clear temp1 temp2
        end
        for iTrl = 1:size(data2,3)
            temp1 = data2(:,:,iTrl);
            temp2 = nanmean(temp1(logical(zmapthresh)),1);
            if ~isnan(temp2)
                cond2 = [cond2 temp2];
            end
            clear temp1 temp2
        end

        % store indiv subj data
        cond1_grp_mn(iSubj)  = nanmean(cond1);
        cond1_grp_sem(iSubj) = nanstd(cond1,0,2)/sqrt(length(cond1));
        
        cond2_grp_mn(iSubj)  = nanmean(cond2);
        cond2_grp_sem(iSubj) = nanstd(cond2,0,2)/sqrt(length(cond2));
        
        % encoding stats
        reps = 1000;
        adata = cond1;
        bdata = cond2;
        p2    = permutation_unpaired(adata, bdata, reps)
        pvals_subjs(iSubj,:) = [p2];
    else strcmp('tuning_correct', exp_type)
        cond1 = []; cond2 = []; cond3 = []; cond4 = [];
        for iTrl = 1:size(data1,3)
            temp1 = data1(:,:,iTrl);
            temp2 = nanmean(temp1(logical(zmapthresh)),1);
            if ~isnan(temp2)
                cond1 = [cond1 temp2];
            end
            clear temp1 temp2
        end
        for iTrl = 1:size(data2,3)
            temp1 = data2(:,:,iTrl);
            temp2 = nanmean(temp1(logical(zmapthresh)),1);
            if ~isnan(temp2)
                cond2 = [cond2 temp2];
            end
            clear temp1 temp2
        end
        for iTrl = 1:size(data3,3)
            temp1 = data3(:,:,iTrl);
            temp2 = nanmean(temp1(logical(zmapthresh)),1);
            if ~isnan(temp2)
                cond3 = [cond3 temp2];
            end
            clear temp1 temp2
        end
        for iTrl = 1:size(data4,3)
            temp1 = data4(:,:,iTrl);
            temp2 = nanmean(temp1(logical(zmapthresh)),1);
            if ~isnan(temp2)
                cond4 = [cond4 temp2];
            end
            clear temp1 temp2
        end  
        
        cond1_grp_mn(iSubj)  = nanmean(cond1);
        cond1_grp_sem(iSubj) = nanstd(cond1,0,2)/sqrt(length(cond1));
        
        cond2_grp_mn(iSubj)  = nanmean(cond2);
        cond2_grp_sem(iSubj) = nanstd(cond2,0,2)/sqrt(length(cond2));
        
        cond3_grp_mn(iSubj)  = nanmean(cond3);
        cond3_grp_sem(iSubj) = nanstd(cond3,0,2)/sqrt(length(cond3));
        
        cond4_grp_mn(iSubj)  = nanmean(cond4)
        cond4_grp_sem(iSubj) = nanstd(cond4,0,2)/sqrt(length(cond4));
      
        % save for source data
           cond1_save{iSubj}=cond1;
           cond2_save{iSubj}=cond2;
           cond3_save{iSubj}=cond3;
           cond4_save{iSubj}=cond4;
           
            % retrieval stats
            reps  = 500;
            adata = cond3;
            bdata = cond1;
            p1     = permutation_unpaired(adata, bdata, reps)
            
            adata = cond3;
            bdata = cond2;
            p2    = permutation_unpaired(adata, bdata, reps)
            
            adata = cond3;
            bdata = cond4;
            p3    = permutation_unpaired(adata, bdata, reps)
            
            alpha = .05;
            [p_fdr, p_masked] = fdr( [p1 p2 p3], alpha);
            pvals_subjs(iSubj,1:4) = [p1 p2 p3 p_fdr];
    end
clear data1 data2 data3 data4
end
%%
plot_4_conds =  ''
figure
hold on
colorList = {[0 0 0],[1 0 0],[0 1 0],[1 0.4 0],[.7 0.5 0.25],[0.7 0.2 0.1],[0.3 0.7 0.4],[0 0.5 0.5]};

if strcmp('encoding', exp_type)
    labels =  {'lure +' 'lure -'};
    hold on
    for iSubj = 1:length(cond2_grp_mn)
        plot( [cond1_grp_mn(iSubj) cond2_grp_mn(iSubj)], colorList{iSubj}, 'MarkerFaceColor', color{iSubj}, 'LineWidth', 1)
        errorbar(1:2,[cond1_grp_mn(iSubj) cond2_grp_mn(iSubj)],[cond1_grp_sem(iSubj) cond2_grp_sem(iSubj)],colorList{iSubj},'LineWidth',1)
    end
else

    for iSubj = 1:length(cond2_grp_mn)
        if isnan(cond2_grp_mn(iSubj) )
            continue
        end
        if strcmp('yes',plot_4_conds)
        plot( [cond1_grp_mn(iSubj) cond2_grp_mn(iSubj)  cond3_grp_mn(iSubj)  cond4_grp_mn(iSubj)], 'Color', colorList{iSubj}, 'LineWidth', 2)
         errorbar(1:4,[cond1_grp_mn(iSubj) cond2_grp_mn(iSubj) cond3_grp_mn(iSubj)  cond4_grp_mn(iSubj)],[cond1_grp_sem(iSubj) cond2_grp_sem(iSubj) cond3_grp_sem(iSubj) cond4_grp_sem(iSubj)],'Color',colorList{iSubj},'LineWidth',2)
            labels = {'repeat' 'lure -' 'lure +' 'new'};
        else
            plot( [ cond2_grp_mn(iSubj)  cond3_grp_mn(iSubj) ], 'Color',colorList{iSubj},'LineWidth',2)
            errorbar(1:2,[cond2_grp_mn(iSubj) cond3_grp_mn(iSubj) ],[ cond2_grp_sem(iSubj) cond3_grp_sem(iSubj)], 'Color',colorList{iSubj},'LineWidth',2)
            labels = { 'lure -' 'lure +'};
        end
    end
end
cd('/mnt/yassamri/iEEG/sandra/GroupFigures')
xlim([0 length(labels)+1])
y = ylim;
%ylim([-0.5 0.5])
yticks([y(1) 0 y(2)])
print('-clipboard','-dbitmap')

%% plot two conds at a time
% repeat vs. lure+
figure;hold on
for iSubj = 1:length(cond2_grp_mn)
    if isnan(cond2_grp_mn(iSubj) )
        continue
    end
    
    plot( [cond1_grp_mn(iSubj) cond3_grp_mn(iSubj)  ], 'Color', colorList{iSubj}, 'LineWidth', 2)
    errorbar(1:2,[cond1_grp_mn(iSubj) cond3_grp_mn(iSubj) ],[cond1_grp_sem(iSubj) cond3_grp_sem(iSubj) ],'Color',colorList{iSubj},'LineWidth',2)
     
end
labels = {'repeat+', 'lure+'}
set(gca, 'XTick', 1:2, 'XTickLabel', labels(1:2),'XTickLabelRotation',45, 'FontSize', 12, 'FontWeight', 'bold')
xlim([0 3])
title(reg)
print('-clipboard','-dbitmap')

% lure- vs. lure+
figure;hold on
for iSubj = 1:length(cond2_grp_mn)
    if isnan(cond2_grp_mn(iSubj) )
        continue
    end
    
    plot( [cond2_grp_mn(iSubj) cond3_grp_mn(iSubj)  ], 'Color', colorList{iSubj}, 'LineWidth', 2)
    errorbar(1:2,[cond2_grp_mn(iSubj) cond3_grp_mn(iSubj) ],[cond2_grp_sem(iSubj) cond3_grp_sem(iSubj) ],'Color',colorList{iSubj},'LineWidth',2)
     
end
labels = {'lure-', 'lure+'}
set(gca, 'XTick', 1:2, 'XTickLabel', labels(1:2),'XTickLabelRotation',45, 'FontSize', 12, 'FontWeight', 'bold')
xlim([0 3])
title(reg)
print('-clipboard','-dbitmap')

% new+ vs. lure+
figure;hold on
for iSubj = 1:length(cond2_grp_mn)
    if isnan(cond2_grp_mn(iSubj) )
        continue
    end
    
    plot( [cond4_grp_mn(iSubj) cond3_grp_mn(iSubj)  ], 'Color', colorList{iSubj}, 'LineWidth', 2)
    errorbar(1:2,[cond4_grp_mn(iSubj) cond3_grp_mn(iSubj) ],[cond4_grp_sem(iSubj) cond3_grp_sem(iSubj) ],'Color',colorList{iSubj},'LineWidth',2)
     
end
labels = {'new+', 'lure+'}
set(gca, 'XTick', 1:2, 'XTickLabel', labels(1:2),'XTickLabelRotation',45, 'FontSize', 12, 'FontWeight', 'bold')
xlim([0 3])
title(reg)
print('-clipboard','-dbitmap')
%% Indiv subj average across trials in a given chan
[cond1a,cond2a,cond3a,cond4a] = get_power_subj_elecs(subj_list,exp_type,reg,ref,baseline,lock);
cond2 = cat(3,cond2a{:});
cond3 = cat(3,cond3a{:});
reset = nanmean(nanmean(cond2(:,0.3*fs:pre_stim*fs,:),3),2);
cond2 = cond2 - reset;
reset = nanmean(nanmean(cond3(:,0.3*fs:pre_stim*fs,:),3),2);
cond3 = cond3 - reset;
if strcmp('tuning_correct',exp_type)
    cond1 = cat(3,cond1a{:});
    cond4 = cat(3,cond4a{:});
    reset = nanmean(nanmean(cond1(:,0.3*fs:pre_stim*fs,:),3),2);
    cond1 = cond1 - reset;
    reset = nanmean(nanmean(cond4(:,0.3*fs:pre_stim*fs,:),3),2);
    cond4 = cond4 - reset;
end

cond1_grp_mn= nan(length(subj_list),1);
cond2_grp_mn= nan(length(subj_list),1);
cond3_grp_mn= nan(length(subj_list),1);
cond4_grp_mn= nan(length(subj_list),1);




for iSubj = 1:length(subj_list) %loop thru subj
    temp_data_per_subj_1 = nanmean(cond1a{iSubj},3);
    cond1_grp_mn(iSubj) = nanmean(temp_data_per_subj_1(logical(zmapthresh_temp)));
    
    temp_data_per_subj_2 = nanmean(cond2a{iSubj},3);
    cond2_grp_mn(iSubj) = nanmean(temp_data_per_subj_2(logical(zmapthresh_temp)));
    
    temp_data_per_subj_3 = nanmean(cond3a{iSubj},3);
    cond3_grp_mn(iSubj) = nanmean(temp_data_per_subj_3(logical(zmapthresh_temp)));
    
    temp_data_per_subj_4 = nanmean(cond4a{iSubj},3);
    cond4_grp_mn(iSubj) = nanmean(temp_data_per_subj_4(logical(zmapthresh_temp)));
    
    clear temp_data_per_subj_1 temp_data_per_subj_2 temp_data_per_subj_3 temp_data_per_subj_4
end

figure;hold on
colorList = {[0 0 0],[1 0 0],[0 1 0],[1 0.4 0],[.7 0.5 0.25],[0.7 0.2 0.1],[0.3 0.7 0.4],[0 0.5 0.5]};
for iSubj = 1:length(cond2_grp_mn)
    if isnan(cond2_grp_mn(iSubj) )
        continue
    end
   plot( [cond1_grp_mn(iSubj) cond3_grp_mn(iSubj) ], 'Color', colorList{iSubj}, 'LineWidth', 2)
end
labels = {'repeat+', 'lure+'}
set(gca, 'XTick', 1:2, 'XTickLabel', labels(1:2),'XTickLabelRotation',45, 'FontSize', 12, 'FontWeight', 'bold')
xlim([0 3])
title(reg)
print('-clipboard','-dbitmap')

figure;hold on
colorList = {[0 0 0],[1 0 0],[0 1 0],[1 0.4 0],[.7 0.5 0.25],[0.7 0.2 0.1],[0.3 0.7 0.4],[0 0.5 0.5]};
for iSubj = 1:length(cond2_grp_mn)
    if isnan(cond2_grp_mn(iSubj) )
        continue
    end
   plot( [cond2_grp_mn(iSubj) cond3_grp_mn(iSubj) ], 'Color', colorList{iSubj}, 'LineWidth', 2)
end
labels = {'lure-', 'lure+'}
set(gca, 'XTick', 1:2, 'XTickLabel', labels(1:2),'XTickLabelRotation',45, 'FontSize', 12, 'FontWeight', 'bold')
xlim([0 3])
title(reg)
print('-clipboard','-dbitmap')

figure;hold on
colorList = {[0 0 0],[1 0 0],[0 1 0],[1 0.4 0],[.7 0.5 0.25],[0.7 0.2 0.1],[0.3 0.7 0.4],[0 0.5 0.5]};
for iSubj = 1:length(cond2_grp_mn)
    if isnan(cond2_grp_mn(iSubj) )
        continue
    end
   plot( [cond4_grp_mn(iSubj) cond3_grp_mn(iSubj) ], 'Color', colorList{iSubj}, 'LineWidth', 2)
end
labels = {'new+', 'lure+'}
set(gca, 'XTick', 1:2, 'XTickLabel', labels(1:2),'XTickLabelRotation',45, 'FontSize', 12, 'FontWeight', 'bold')
xlim([0 3])
title(reg)
print('-clipboard','-dbitmap')



%% using indiv subj cluster - MTL stim onset
cd('/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_onset')
load(['subj39spectrograms.mat'])

for subj = 1:length(subj_list)
    % load zmap
    cd('/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_onset/cluster_matrices')
    load(['subj_' subj_list{subj} '_reg_MTLzmapthresh.mat'])
    
    if strcmp('84', subj_list{subj})
        zmap_weight = zeros(size(MTL_cond1,1),size(MTL_cond1,2));
        zmap_weight(:,301:1301) = zmapthresh;

    else
        zmapthresh_temp = zeros(size(MTL_cond1,1),size(MTL_cond1,2));
        zmapthresh_temp(:,301:1301) = logical(zmapthresh);
    end
    
    cd('/mnt/yassamri/iEEG/sandra/groupdata_spectrograms_onset')
    load(['subj' subj_list{subj} 'spectrograms.mat'])
    
    temp1 = [];
    temp1 = nanmean(MTL_cond1,3);
    
    temp2 = [];
    temp2 = nanmean(MTL_cond2,3);
    
    temp3 = [];
    temp3 = nanmean(MTL_cond3,3); 
    
    temp4 = [];
    temp4 = nanmean(MTL_cond4,3);  
    
    if strcmp('84',subj_list{subj})
    % cond 1
    temp1_a = temp1.*zmap_weight;
    cond1_grp_val(subj) = [nanmean(temp1_a(logical(zmapthresh_temp)))];    
    % cond 2  
    temp2_a = temp2.*zmap_weight;
    cond2_grp_val(subj) = [nanmean(temp2_a(logical(zmapthresh_temp)))];    
    % cond 3
    temp3_a = temp3.*zmap_weight;
    cond3_grp_val(subj) = [nanmean(temp3_a(logical(zmapthresh_temp)))];    
    % cond 4
    temp4_a = temp4.*zmap_weight;
    cond4_grp_val(subj) = [nanmean(temp4_a(logical(zmapthresh_temp)))]; 
    else
    % cond 1
    cond1_grp_val(subj) = [nanmean(temp1(logical(zmapthresh_temp)))];    
    % cond 2  
    cond2_grp_val(subj) = [nanmean(temp2(logical(zmapthresh_temp)))];    
    % cond 3
    cond3_grp_val(subj) = [nanmean(temp3(logical(zmapthresh_temp)))];    
    % cond 4
    cond4_grp_val(subj) = [nanmean(temp4(logical(zmapthresh_temp)))];   
    end
end


bar_vector_mn1  = [nanmean(cond1_grp_val) nanmean(cond2_grp_val) nanmean(cond3_grp_val) nanmean(cond4_grp_val)];
bar_vector_std1 = [nanstd(cond1_grp_val,0,2) nanstd(cond2_grp_val,0,2) nanstd(cond3_grp_val,0,2) nanstd(cond4_grp_val,0,2)];

% plot
figure
hold on
bar([bar_vector_mn1])         
errorbar(1:4,bar_vector_mn1,bar_vector_std1/sqrt(length(subj_list)), 'rx')
set(gca, 'XTick', 1:4, 'XTickLabel', {'repeat' 'lure -' 'lure +' 'new'},'XTickLabelRotation',45)
ylabel('MTL power')
set(gca, 'FontSize', 16, 'FontWeight', 'bold') 
title([[]])
ylim([])
ylim([])
 

% run stats
reps  = 10000;
adata = cond3_grp_val; 
bdata = cond1_grp_val;
p1    = permutation_paired(adata, bdata, reps)

adata = cond3_grp_val; 
bdata = cond2_grp_val;
p2    = permutation_paired(adata, bdata, reps)

adata = cond3_grp_val; 
bdata = cond4_grp_val;
p3    = permutation_paired(adata, bdata, reps)

alpha = .05
pvals = [p1 p2 p3]
[p_fdr, p_masked] = fdr(pvals, alpha)

