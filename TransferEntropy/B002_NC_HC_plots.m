figure
group_plots      = 'yes'
if strcmp('yes',group_plots)
    subjLength=1;
else
    subjLength=8;
end
for iSubj        = 1:subjLength% change this for indiv subj line color
    clearvars -except conda condb condc cond1 cond2 cond3 cond4 cond5 cond5 iSubj group_plots ...
        reg1_reg2_LureMinus reg1_reg2_LurePlus;
subj_list_original = {'39' '44'  '57'  '63' '66' '84' '85' '87' }
subj             = subj_list_original{iSubj}

%;
hold on
sig_chans        =''
%'39' '44' '57'  '63' '66' '84' '85' '87' }
lock             = 'onset'
phase            = 'retrieval'
%phase            = 'encoding'
freq_range       = 'deltatheta'
fpass            = [4 5]
fn_ext           = 'cue_resp' %
plot_4_conds     = ''
lock             = 'onset'
times1           = {'0' '0' '0.5'}
times2           = {'2' '1' '1.5'}
plot_type        = 'main' % main phase
bidirec          = 0.5
colorList        = {'r-o', 'c-o', 'g-o', 'k-o', 'm-o',  'y-o','b-o', 'k-o'};
color            = {'r', 'c', 'g', 'b', 'm','k', 'y', 'k'};
condCounterList  =2; % default is only compare lure+ vs. lure- unles plot 4 conds

if strcmp('encoding',phase)
    time1 = times1{1}
    time2 = times2{1}
    lock = 'onset'
elseif strcmp('retrieval',phase) && strcmp('onset', lock)
    time1 = times1{2}
    time2 = times2{2}
elseif strcmp('retrieval',phase) && strcmp('response', lock)
    time1 = times1{3}
    time2 = times2{3}
end
reg1_list   = {'TEMP' };%'FRO' 'TEMP'
reg2_list   = {'HC'};



addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
cd (['/mnt/yassamri/iEEG/sandra/PTE_results/' fn_ext '/' freq_range '/' phase '/' lock '/' num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1) '_' num2str(time2) 'sec'])

mtx_sz = 40
% cond 2: lure -
reg1_reg2_conda = cell(1,mtx_sz);
% cond 3: lure+
reg1_reg2_condb = cell(1,mtx_sz);
% cond1: repeat
reg1_reg2_condc = cell(1,mtx_sz);
% cond 4: new
reg1_reg2_condd = cell(1,mtx_sz);

i = 0;
for iReg1 = 1:length(reg1_list)
    
    % get NC site
    reg1_name = reg1_list{iReg1};
    reg2_name = reg2_list{1};
    
    if strcmp('yes',group_plots)
                     subj_list = {'39' '44' '57'  '63' '66' '84' '85' '87' }
                    % subj_list = { '44' '57'  '63'  '85'  }

    else
        subj_list{1} = subj;
    end
    % get all data from each subj
    for sub_counter = 1:length(subj_list)
        i = i+1;

        if isfile([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2' sig_chans '.mat' ])
            if strcmp('retrieval',phase)
                % load lure -
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2' sig_chans ])
                reg1_reg2_conda{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                
                % load lure +
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3' sig_chans ])
                reg1_reg2_condb{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                
                if strcmp('yes',plot_4_conds)
                    condCounterList = [1 2 3];

                    % load repeat
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' sig_chans ])
                    reg1_reg2_condc{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                    
                    % load new
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond4' sig_chans ])
                    reg1_reg2_condd{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                end
            else
                % load lure -
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2' sig_chans ])
                reg1_reg2_conda{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                
                % load lure +
                load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' sig_chans ])
                reg1_reg2_condb{i} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
            end
        end
    end
end

reg1_reg2_conda = cat(1,reg1_reg2_conda{:}); % lure-
reg1_reg2_condb = cat(1,reg1_reg2_condb{:}); % lure+
reg1_reg2_condc = cat(1,reg1_reg2_condc{:}); % repeat
reg1_reg2_condd = cat(1,reg1_reg2_condd{:}); % new

% save data for sourcedata
reg1_reg2_LureMinus{iSubj} = reg1_reg2_conda; % lure-
reg1_reg2_LurePlus {iSubj} = reg1_reg2_condb; % lure+

% run stats for main figure (comparing conditions withnin phases), and plot
%if strcmp('main',plot_type)
reps = 1000;
alpha = 0.05;
cntr = 0
for iCond = condCounterList 
    cntr=cntr+1
    
    if iCond ==1
        % lure+ vs. repeat
        adata =reg1_reg2_condb;% lure+
        bdata =reg1_reg2_condc;% repeat
        
    elseif iCond ==2
        % lure+ vs. lure-
        adata =reg1_reg2_condb;% lure+
        bdata =reg1_reg2_conda;% lure-
        
    elseif iCond ==3
        % lure+ vs. new+
        adata =reg1_reg2_condb;% lure+
        bdata =reg1_reg2_condd;% new+
    end
    [p,effect] = permutation_paired(adata, bdata, reps);
    p_val(cntr) =  p;
     effect(cntr) =  effect;
end
[p_fdr, p_masked] = fdr( p_val, alpha);
[ p_val p_fdr]

    % plot
    mn = -.012;
    mx = .012;
    
    
    if strcmp('encoding', phase)
        if strcmp('yes', group_plots)
            %subplot(1,3,1)
            hold on
            bar([nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec])
            errorbar(1:2,[nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec],[nanstd(reg1_reg2_conda)/sqrt(length(reg1_reg2_conda)) nanstd(reg1_reg2_condb)/sqrt(length(reg1_reg2_condb))], 'rx')
            set(gca, 'XTick', 1:2, 'XTickLabel', {'lure -' 'lure +'},'XTickLabelRotation',45)
            title([phase ' ' lock])
            set(gca, 'FontSize', 16, 'FontWeight', 'bold')
            % ylim([mn mx])
        else
            hold on
            plot([nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec], colorList{iSubj}, 'MarkerFaceColor', color{iSubj}, 'LineWidth', 1)
            errorbar(1:2,[nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec],[nanstd(reg1_reg2_conda)/sqrt(length(reg1_reg2_conda)) nanstd(reg1_reg2_condb)/sqrt(length(reg1_reg2_condb))],colorList{iSubj},'LineWidth',1)
            set(gca, 'XTick', 1:2, 'XTickLabel', {'lure -' 'lure +'},'XTickLabelRotation',45)
            title([phase ' ' lock])
            set(gca, 'FontSize', 16, 'FontWeight', 'bold')
            
        end
    elseif strcmp('retrieval', phase)
        
        if strcmp('onset', lock)
            if strcmp('yes', group_plots)
                % subplot(1,3,2)
                hold on
                if strcmp('yes',plot_4_conds)
                    mn_vector  = [nanmean(reg1_reg2_condc)-bidirec nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec nanmean(reg1_reg2_condd)-bidirec];
                    std_vector = [nanstd(reg1_reg2_condc)/sqrt(length(reg1_reg2_condc)) nanstd(reg1_reg2_conda)/sqrt(length(reg1_reg2_conda))...
                        nanstd(reg1_reg2_condb)/sqrt(length(reg1_reg2_condb)) nanstd(reg1_reg2_condd)/sqrt(length(reg1_reg2_condd))];
                    iCond =4;label = {'repeat +' 'lure -' 'lure +' 'new +'};
                else
                    mn_vector  = [ nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec ];
                    std_vector = [ nanstd(reg1_reg2_conda)/sqrt(length(reg1_reg2_conda)) nanstd(reg1_reg2_condb)/sqrt(length(reg1_reg2_condb)) ];
                    iCond =2; label = { 'lure -' 'lure +' };
                end
                bar(mn_vector)
                errorbar(1:iCond,mn_vector,std_vector, 'rx')
                set(gca, 'XTick', 1:iCond, 'XTickLabel', label ,'XTickLabelRotation',45)
                title([phase ' ' lock])
                set(gca, 'FontSize', 16, 'FontWeight', 'bold')
                %    ylim([mn mx])
            else
                hold on
                plot([nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec], colorList{iSubj}, 'MarkerFaceColor', color{iSubj}, 'LineWidth', 1)
                errorbar(1:2,[nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec],[nanstd(reg1_reg2_conda)/sqrt(length(reg1_reg2_conda)) nanstd(reg1_reg2_condb)/sqrt(length(reg1_reg2_condb))],colorList{iSubj},'LineWidth',1)
                set(gca, 'XTick', 1:2, 'XTickLabel', {'lure -' 'lure +'},'XTickLabelRotation',45)
                title([phase ' ' lock])
                set(gca, 'FontSize', 16, 'FontWeight', 'bold')
            end
            
        elseif strcmp('response', lock)
           
            if strcmp('yes',plot_4_conds)
                mn_vector  = [nanmean(reg1_reg2_condc)-bidirec nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec nanmean(reg1_reg2_condd)-bidirec];
                std_vector = [nanstd(reg1_reg2_condc)/sqrt(length(reg1_reg2_condc)) nanstd(reg1_reg2_conda)/sqrt(length(reg1_reg2_conda))...
                    nanstd(reg1_reg2_condb)/sqrt(length(reg1_reg2_condb)) nanstd(reg1_reg2_condd)/sqrt(length(reg1_reg2_condd))];
                iCond =4;
            else
                mn_vector  = [ nanmean(reg1_reg2_conda)-bidirec nanmean(reg1_reg2_condb)-bidirec ];
                std_vector = [ nanstd(reg1_reg2_conda)/sqrt(length(reg1_reg2_conda)) nanstd(reg1_reg2_condb)/sqrt(length(reg1_reg2_condb)) ];
                iCond =2;
            end
            bar(mn_vector)
            errorbar(1:iCond,mn_vector,std_vector, 'rx')
            set(gca, 'XTick', 1:iCond, 'XTickLabel', { 'lure -' 'lure +' },'XTickLabelRotation',45)
            title([phase ' ' lock])
            set(gca, 'FontSize', 16, 'FontWeight', 'bold')
            %   ylim([mn mx])
        end
    end
end
cd('/mnt/yassamri/iEEG/sandra/GroupFigures')

%% plot encoding, retrieval onset, & retrieval response (all lure +), and run stats

mn = -0.015;
mx = 0.01;
if strcmp('phase',plot_type)
    if strcmp('encoding', phase)
        conda = reg1_reg2_condb;
    elseif strcmp('retrieval', phase)
        if strcmp('onset', lock)
            condb = reg1_reg2_condb;condc=[];
        elseif strcmp('response', lock)
            condc = reg1_reg2_condb;
        end
    end
    figure
    hold on
    bar([nanmean(conda) nanmean(condb) nanmean(condc)]-bidirec)
    errorbar(1:3,[nanmean(conda) nanmean(condb) nanmean(condc)]-bidirec,[nanstd(conda)/sqrt(length(conda)) nanstd(condb)/sqrt(length(condb)) nanstd(condc)/sqrt(length(condc))], 'rx')
    set(gca, 'XTick', 1:3, 'XTickLabel', {'encode lure+' 'ret onst lure+' 'ret resp lure+'},'XTickLabelRotation',45)
    set(gca, 'FontSize', 16, 'FontWeight', 'bold')
    ylim([mn mx])
    
    %% stats for phase plots
    cond1_stat = condc;
    cond2_stat = condb;
    
    real_diff  = mean(cond1_stat - cond2_stat);
    both_conds = cat(1, cond1_stat, cond2_stat);
    n_permutes = 1000;
    null_dist_of_diff = nan(1,n_permutes);
    for permi = 1:n_permutes
        shuf_cntr = [];
        for ii = 1:size(cond1_stat,1)
            shuf_cntr = [shuf_cntr randi([0 1],1)];
        end
        % fake conds
        cond1 =[cond1_stat(logical(shuf_cntr))' cond2_stat(logical(shuf_cntr==0))'];
        cond2 =[cond2_stat(logical(shuf_cntr))' cond1_stat(logical(shuf_cntr==0))'];
        
        null_dist_of_diff (permi) = mean(cond1-cond2);
    end
    p1a = length(find(null_dist_of_diff>real_diff))/n_permutes;
    p1b = length(find(null_dist_of_diff<real_diff))/n_permutes;
    
    
    if p1a < 0.05 || p1b<0.05
        if real_diff<0
            disp(['sig p = ' num2str(p1b) ' '  '!!!' ])
        else
            disp(['sig p = ' num2str(p1a) ' '  '!!!' ])
        end
    end
    
end
%ylim([-.01 .04])
y=ylim;
if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
end
print('-clipboard','-dbitmap')

