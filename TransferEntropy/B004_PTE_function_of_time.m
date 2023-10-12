clear all;clc

fpass            = [4 5];
lock             = 'onset';
phase            = 'retrieval';% retrieval encoding
group_plots      = 'yes';
subj             = '';
plot_4_conds     = '';
freq_range       = 'deltatheta';
bidirec          = 0.5;
fs               = 500;

% define time windows
if strcmp('encoding',phase)
%     time1_list = round([0:0.25:1.5],2)
%     time2_list = round([0.5:0.25:2],2)
%     
%             time1_list = round([0:0.01:1.5],2);
%         time2_list = round([0.5:0.01:2],2);
%         
%         
%     labels           = {'0-.5' '0.05-0.55' '0.1-0.60' '0.15-0.65' '0.20-0.70'...
%         '0.25-0.75' '0.30-0.80' '0.35-0.85' '0.40-0.90' '0.45-0.95'...
%         '0.50-1' '0.55-1.05' '0.60-1.10' '0.65-1.15' '0.70-1.20'...
%         '0.75-1.2' '0.80-1.30' '0.85-1.35' '0.90-1.40' '0.95-1.45'...
%         '1-.5' '1.05-1.55' '1.10-1.60' '1.15-1.65' '1.20-1.70'...
%         '1.25-1.75' '1.30-1.80' '1.35-1.85' '1.40-1.90' '1.45-1.95' '1.50-2'}
%     
%     
%     time1_list = round([0:0.05:1.5],2);
%     time2_list = round([0.5:0.05:2],2);
%     
%     
% used for supplement fig S2e
         time1_list = [0:0.01:1.5];
         time2_list = [0.5:0.01:2];
         labels           = {'0-.5' '0.1-0.6'  '0.2-0.7' '0.3-0.8' '0.4-0.9' '0.5-1' '0.6-1' '0.7-1.2' ...
           '0.8-1.3' '0.9-1.4' '1-1.5' '1.1-1.6' '1.2-1.7' '1.3-1.8' '1.4-1.9' '1.5-2'}
       
       %%
    %  labels           = {'0-.5' '0.25-0.75' '0.5-1' '0.75-1.25' '1-1.5' '1.25-1.75' '1.5-2' '1.5-2'}
elseif strcmp('retrieval', phase) && strcmp('onset', lock)
    %     time1_list = round([0:0.25:0.5],2)
    %     time2_list = round([0.5:0.25:1],2)
    
    %       time1_list = round([0:26/fs:0.5],2) % 90% overlap
    %       time2_list = round([0.5:26/fs:1],2)
     
%     time1_list = round([0:0.05:0.5],2) % 90% overlap
%     time2_list = round([0.5:0.05:1],2)
     labels           = {'0-.5' '0.1-0.6'  '0.2-0.7' '0.3-0.8' '0.4-0.9' '0.5-1'}
    
       time1_list = 0:0.01:0.5; % 90% overlap
       time2_list = 0.5:0.01:1;
    % labels           = {'0-.5' '0.25-0.75' '0.5-1'}
elseif strcmp('retrieval', phase) && strcmp('response', lock)
    time1_list = round([0.5:0.25:1],2)
    time2_list = round([1:0.25:1.5],2)
    labels           = {'0 - -0.5' '-0.75 - -0.25' '-0.5 - 0'}
end

reg1_list   = { 'OFC' 'FRO' 'TEMP'}; % 'OFC' 'FRO'
reg2_list   = {'HC'};
fn_ext      = 'cue_resp';
table       =  []; % for mixed effects model in R

for iEpoch = 1:length(time1_list)
    
    addpath('/mnt/yassamri/iEEG/sandra/analysis_pipeline_final')
    cd (['/mnt/yassamri/iEEG/sandra/PTE_results/' fn_ext '/' freq_range '/' phase '/' lock '/' num2str(fpass(1)) '_' num2str(fpass(2)) 'Hz_' num2str(time1_list(iEpoch)) '_' num2str(time2_list(iEpoch)) 'sec'])
    
    mtx_sz = 40;
    % cond 2: lure -
    reg1_reg2_conda = cell(1,mtx_sz);
    % cond 3: lure+
    reg1_reg2_condb = cell(1,mtx_sz);
    % cond1: repeat
    reg1_reg2_condc = cell(1,mtx_sz);
    % cond 4: new
    reg1_reg2_condd = cell(1,mtx_sz);
    
    iCounter = 0;
    % making table for mixed effects model
    subjIDCol = []; elecPairCol=[]; condCol=[]; reg1CodeCol = [];reg2CodeCol = []; timeCol=[];
    
    for iReg1 = 1:length(reg1_list)
        
        % get NC site
        reg1_name = reg1_list{iReg1};
        reg2_name = reg2_list{1};
        
        % group vs. indiv subj
        if strcmp('yes',group_plots)
             subj_list = {'39' '44' '57'  '63' '66' '84' '85' '87'};
        else
            subj_list{1} =subj;
        end
        
        % get all data from each subj
        for sub_counter = 1:length(subj_list)
            if isfile([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2.mat' ])
                iCounter = iCounter+1;
                if strcmp('retrieval',phase)
                    % load lure -
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2' ])
                    reg1_reg2_conda{iCounter} = PTE_ch1_to_ch2_norm';
%                     subjIDCol   = [subjIDCol repmat(sub_counter,[1 length(PTE_ch1_to_ch2_norm)])];
%                     elecPairCol = [elecPairCol elecPair]; %NC-HC appended
%                     condCol     = [condCol repmat(2,[1 length(PTE_ch1_to_ch2_norm)])];
%                     reg2CodeCol = [reg2CodeCol reg2Code-1];%HC: subtract 1 to start at 0
%                     reg1CodeCol = [reg1CodeCol repmat(iReg1-1,[1 length(PTE_ch1_to_ch2_norm)])];%NC
%                     timeCol     = [timeCol repmat(iEpoch,[1 length(PTE_ch1_to_ch2_norm)])];
                    
                    clear PTE_ch1_to_ch2_norm
                    
                    % load lure +
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3' ])
                    reg1_reg2_condb{iCounter} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                    
                    if strcmp('yes',plot_4_conds)
                    % load repeat
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' ])
                    reg1_reg2_condc{iCounter} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                    
                    % load new
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond4' ])
                    reg1_reg2_condd{iCounter} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                        
                    end
                else
                    % load lure -
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond2' ])
                    reg1_reg2_conda{iCounter} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                    
                    % load lure +
                    load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' ])
                    reg1_reg2_condb{iCounter} = PTE_ch1_to_ch2_norm';clear PTE_ch1_to_ch2_norm
                end
            else
                continue
            end
        end
    end
    
    reg1_reg2_conda = cat(1,reg1_reg2_conda{:}); % lure-
    reg1_reg2_condb = cat(1,reg1_reg2_condb{:}); % lure+


    conda_val_time_mn_store (iEpoch,:) = reg1_reg2_conda;
    condb_val_time_mn_store (iEpoch,:) = reg1_reg2_condb;

    conda_val_time_mn (iEpoch) = nanmean(reg1_reg2_conda);
    condb_val_time_mn (iEpoch) = nanmean(reg1_reg2_condb);
    
    conda_val_time_std (iEpoch) = nanstd(reg1_reg2_conda)/sqrt(length(reg1_reg2_conda));
    condb_val_time_std (iEpoch) = nanstd(reg1_reg2_condb)/sqrt(length(reg1_reg2_condb));

    if strcmp('yes',plot_4_conds)
        reg1_reg2_condc = cat(1,reg1_reg2_condc{:}); % repeat
        reg1_reg2_condd = cat(1,reg1_reg2_condd{:}); % new
        condc_val_time_mn_store (iEpoch,:) = reg1_reg2_condc;
        condd_val_time_mn_store (iEpoch,:) = reg1_reg2_condd;
        condc_val_time_mn (iEpoch) = nanmean(reg1_reg2_condc);
        condd_val_time_mn (iEpoch) = nanmean(reg1_reg2_condd);
        condc_val_time_std (iEpoch) = nanstd(reg1_reg2_condc)/sqrt(length(reg1_reg2_condc));
        condd_val_time_std (iEpoch) = nanstd(reg1_reg2_condd)/sqrt(length(reg1_reg2_condd));
    end
    
    %table = [table; subjIDCol' condCol' reg1CodeCol' reg2CodeCol' elecPairCol' timeCol'];
    

end

%% plot

figure;hold on
% lure-
plot( conda_val_time_mn-bidirec, '-o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b','LineWidth', 1)
h=errorbar(1:length(time1_list),conda_val_time_mn-bidirec,conda_val_time_std, 'b','LineWidth',2);
h.CapSize = 0;
% lure+
plot( condb_val_time_mn-bidirec, '-o','MarkerFaceColor', 'm', 'MarkerEdgeColor','m','LineWidth', 1)
h=errorbar(1:length(time1_list),condb_val_time_mn-bidirec,condb_val_time_std, 'm','LineWidth',2);
h.CapSize = 0;
if strcmp('yes', plot_4_conds)
    figure;hold on
    plot( condb_val_time_mn-bidirec, '-o','MarkerFaceColor', 'm', 'MarkerEdgeColor','m','LineWidth', 1)
    h=errorbar(1:length(time1_list),condb_val_time_mn-bidirec,condb_val_time_std, 'm','LineWidth',2);
    h.CapSize = 0;
    % repeat+
    plot( condc_val_time_mn-bidirec, 'MarkerFaceColor', 'r', 'LineWidth', 1)
    h=errorbar(1:length(time1_list),condc_val_time_mn-bidirec,condc_val_time_std, 'r','LineWidth',2);
    h.CapSize = 0;
    
    figure;hold on
    plot( condb_val_time_mn-bidirec, '-o','MarkerFaceColor', 'm', 'MarkerEdgeColor','m','LineWidth', 1)
    h=errorbar(1:length(time1_list),condb_val_time_mn-bidirec,condb_val_time_std, 'm','LineWidth',2);
    h.CapSize = 0;
    % new+
    plot( condd_val_time_mn-bidirec, 'MarkerFaceColor', 'g', 'LineWidth', 1)
    h= errorbar(1:length(time1_list),condd_val_time_mn-bidirec,condd_val_time_std, 'k','LineWidth',2);
    h.CapSize = 0;
end
cd('/mnt/yassamri/iEEG/sandra/GroupFigures')
xticks(1:10:length(time1_list))
xticklabels(labels(1:1:end))
%ylim([-.015 .015])
y=ylim;
if y(1)<0 && y(2)>0; yticks([y(1) 0 y(2)])
elseif y(1)<0 && y(2)<0 ||  y(1)>0 && y(2)>0; yticks([y(1) y(2)])
end
print('-clipboard','-dbitmap')
cd('/mnt/yassamri/iEEG/sandra/GroupFigures')
%line([0 length(time1_list)+1], [0 0], 'color', 'k')
%% run anova on observed data
cond1_vals = condb_val_time_mn_store; % should always be b: lure+
cond2_vals = conda_val_time_mn_store; % a: lure-, c: repeat, d: new
cond1 = ones(size(cond1_vals,2),1);
cond2 = 2*ones(size(cond2_vals,2),1);
labels = [cond1' cond2'];
tabel_lable = cell(1,length(time1_list)+1);
tabel_lable{1} = 'cond';
for iCounter = 2:length(time1_list)+1
    tabel_lable{iCounter} = ['t' num2str(iCounter-1)];
end
Time = [1:length(time1_list)]';
x = [cond1_vals' ; cond2_vals'];
t = array2table([labels' x], 'VariableNames',tabel_lable);
rm = fitrm(t,['t1-t' num2str(Time(end)) ' ~ cond'],'WithinDesign',Time); %used for orig ms
sphericityAsump = mauchly(rm) % if significant, use one of the last 3 p-values in the table
% which are adjusted for sphericity
StatsOutputs   = ranova(rm)
etaSquaredtaTime = StatsOutputs{1,1}/sum(StatsOutputs{1:3,1})
etaSquaredtaTimeCondnrx = StatsOutputs{2,1}/sum(StatsOutputs{1:3,1})

% normality assumption
% normplot(conda_val_time_mn)
% for iTime = 1:length(time1_list)
% h(iTime) = adtest(conda_val_time_mn_store(iTime,:)) 
% end

%% follow by independent samples
reps = 1000;alpha = 0.05;
for iCounter = 1:length(time1_list)
adata = cond1_vals (iCounter,:);
bdata = cond2_vals (iCounter,:);
[p_val(iCounter) ~] = permutation_paired(adata, bdata, reps);
end
[p_fdr, p_masked] = fdr( p_val, alpha);
m=p_val(p_masked)'


%% make table to un ICC in R - original
tabel_lable = cell(1,length(time1_list)+2)
tabel_lable{1} = 'id'
tabel_lable{2} = 'cond'
for iCounter = 3:length(time1_list)+2
    tabel_lable{iCounter} = ['t' num2str(iCounter-2)]
end
tForR=array2table([[subjIDCol subjIDCol]' labels' x], 'VariableNames',tabel_lable);

filename = 'PTEFOTdataForR.csv';
writetable(tForR,filename)

%% run anova on shuffled data for reviews
fakeLabels = ones(1,length(labels));
fakeCond2 = randi( length(labels),[1 length(cond1_vals)]);
fakeLabels(fakeCond2)=2;

cond1_valsFakeTimeEpoch = nan(size(cond1_vals));
cond2_valsFakeTimeEpoch = nan(size(cond2_vals));

for iObsv = 1:size(cond1_vals,2) % for each observation, shuffle the label it belongs to
    randTimeEpoch = randi(size(cond1_vals,1),[1 size(cond1_vals,1)] );
    cond1_valsFakeTimeEpoch (:,iObsv) = cond1_vals(randTimeEpoch,iObsv);
    cond2_valsFakeTimeEpoch (:,iObsv) = cond2_vals(randTimeEpoch,iObsv);
    
end
tabel_lable = cell(1,length(time1_list)+1)
tabel_lable{1} = 'cond'
for iCounter = 2:length(time1_list)+1
    tabel_lable{iCounter} = ['t' num2str(iCounter-1)]
end
Time = [1:length(time1_list)]';
x = [cond1_valsFakeTimeEpoch' ; cond2_valsFakeTimeEpoch'];
t = array2table([fakeLabels' x], 'VariableNames',tabel_lable);
rm = fitrm(t,['t1-t' num2str(Time(end)) ' ~ cond'],'WithinDesign',Time);
ranova(rm)
%% make table for mixed effects model in R
cond2Data = conda_val_time_mn_store';
cond3Data = condb_val_time_mn_store';

AllTableCond2 = [table cond2Data(:) ];
AllTableCond3 = [table cond3Data(:) ];
AllTableCond3(:,2) = 3; % change condition to 3

AllTable = [AllTableCond2; AllTableCond3];

tabel_lable = cell(1,4)
tabel_lable{1} = 'id'
tabel_lable{2} = 'cond'
tabel_lable{3} = 'HCreg'
tabel_lable{4} = 'NCreg'
tabel_lable{5} = 'elecPair'
tabel_lable{6} = 'time'
tabel_lable{7} = 'couplingValue'

tForR=array2table(AllTable, 'VariableNames',tabel_lable);

filename = 'PTEFOTdataForLMEinR.csv';
writetable(tForR,filename)
%%
% smoothen average gamma power for each subject
win= round(.01*fs)
conda_val_time_mn_conv = conv(conda_val_time_mn-bidirec, ones(1,win)/win,'same');
condb_val_time_mn_conv = conv(condb_val_time_mn-bidirec, ones(1,win)/win,'same');

figure;hold on
plot( conda_val_time_mn_conv,  'b','LineWidth', 2)
plot( condb_val_time_mn_conv,  'm','LineWidth', 2)
h=errorbar(1:length(time1_list),conda_val_time_mn-bidirec,conda_val_time_std, 'b','LineWidth',2)
h.CapSize = 0
plot( condb_val_time_mn_conv-bidirec, 'o','MarkerFaceColor', 'm', 'MarkerEdgeColor','m','LineWidth', 1)
h=errorbar(1:length(time1_list),condb_val_time_mn-bidirec,condb_val_time_std, 'm','LineWidth',2)
h.CapSize = 0

xtickangle(45)
ylim([-0.012        0.012])
%yticks([-0.01 0 0.005])
xlim([0 length(time1_list)+1])
title(subj)
set(gca, 'FontSize', 14, 'FontWeight', 'bold', 'FontName', 'Arial')
