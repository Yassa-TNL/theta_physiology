% plot condition traces as a function of delay (either normalize or subtraction plots)

clear all
lock  = 'onset'
%lock = 'response'

phase = 'encoding'
%phase = 'retrieval'

reg1_list = {'OFC' 'FRO'  'TEMP' 'CING' 'INS'}; %' 
reg2_list = {'HC'};

addpath('/tmp/yassamri/iEEG/sandra/analysis_pipeline_final')
cd (['/tmp/yassamri/iEEG/sandra/PTE_results/PTE_w_lags/' phase '/' lock])


mtx_sz = 40

reg1 = cell(1,mtx_sz);
reg2 = cell(1,mtx_sz);
reg3 = cell(1,mtx_sz);
reg4 = cell(1,mtx_sz);
reg5 = cell(1,mtx_sz);

for iReg1 = 1:length(reg1_list)
    
    % get NC site
    reg1_name = reg1_list{iReg1};
    reg2_name = reg2_list{1};
    
    % get subj list
    if (strcmp('TEMP',reg1_name) || strcmp('FRO',reg1_name) ||  strcmp('OFC',reg1_name) || strcmp('CING',reg1_name)) && strcmp('CA3',reg2_name)
        subj_list = { '39' '57' '66' '44'};
    elseif (strcmp('TEMP',reg1_name) || strcmp('FRO',reg1_name) ||  strcmp('OFC',reg1_name)) && strcmp('HC',reg2_name)
        subj_list = { '39' '57' '66' '44' '84' '85'};
    elseif strcmp('INS',reg1_name) || strcmp('INS',reg2_name)
        subj_list = {'39' '57' '85'};
    elseif strcmp('CING',reg1_name)  && strcmp('HC',reg2_name)
        subj_list = { '39' '57' '66' '44' '85'};
    elseif (strcmp('OFC',reg1_name) && (strcmp('FRO',reg2_name) || strcmp('TEMP',reg2_name)))
        subj_list = { '39' '57' '66' '44' '84' '63' '85'};
    elseif (strcmp('OFC',reg1_name) || strcmp('FRO',reg1_name)|| strcmp('TEMP',reg1_name))  && (strcmp('CING',reg2_name))
        subj_list = { '39' '57' '66' '44' '63' '85'};
    elseif (strcmp('FRO',reg1_name)  &&  strcmp('TEMP',reg2_name))
        subj_list = { '39' '57' '66' '44' '84' '83' '63' '85'};
    end
    
    
    % get all data from each subj
    i = 0;
    for sub_counter = 1:length(subj_list)
        i = i+1;   
        % load cond3
        if strcmp('encoding',phase)
            load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond1' ])
        else
            load([reg1_name 'vs' reg2_name '_subj'  subj_list{sub_counter} '_cond3' ])
        end
   
    
    if strcmp('retrieval',phase)
        
        if iReg1 ==1
            % reg1{i} = ch1_to_ch2;
            reg1{i} = PTE_ch1_to_ch2_norm;
            reg1{i} = ch1_to_ch2 - ch2_to_ch1;
        elseif iReg1 ==2
            reg2{i} = PTE_ch1_to_ch2_norm;
            reg2{i} = ch1_to_ch2 - ch2_to_ch1;
            
        elseif iReg1 ==3
            reg3{i} = PTE_ch1_to_ch2_norm;
            reg3{i} = ch1_to_ch2 - ch2_to_ch1;
            
        elseif iReg1 ==4
            reg4{i} = PTE_ch1_to_ch2_norm;
            reg4{i} = ch1_to_ch2 - ch2_to_ch1;
            
        elseif iReg1 ==5
            reg5{i} = PTE_ch1_to_ch2_norm;
            reg5{i} = ch1_to_ch2 - ch2_to_ch1;
            
        end
        
    elseif strcmp('encoding',phase)
        
        if iReg1 ==1
            %reg1{i} = PTE_ch1_to_ch2_norm;
            reg1{i} = ch2_to_ch1 - ch1_to_ch2;
            
        elseif iReg1 ==2
            %reg2{i} = PTE_ch1_to_ch2_norm;
            reg2{i} = ch2_to_ch1 - ch1_to_ch2;
            
        elseif iReg1 ==3
           % reg3{i} = PTE_ch1_to_ch2_norm;
            reg3{i} = ch2_to_ch1 - ch1_to_ch2;
            
        elseif iReg1 ==4
            %reg4{i} = PTE_ch1_to_ch2_norm;
            reg4{i} = ch2_to_ch1 - ch1_to_ch2;
            
        elseif iReg1 ==5
            %reg5{i} = PTE_ch1_to_ch2_norm;
            reg5{i} = ch2_to_ch1 - ch1_to_ch2;
            
        end
        
    end
     end
end

reg1=cat(1,reg1{:});
reg2=cat(1,reg2{:});
reg3=cat(1,reg3{:});
reg4=cat(1,reg4{:});
reg5=cat(1,reg5{:});

figure
hold on
x_stps = (0.005:.005:.2)*1000;
stdshade(reg1,.1,'m',x_stps, [],[])
h1 = plot(x_stps, nanmean(reg1,1), '-m', 'LineWidth', 3)

stdshade(reg2,.1,'r',x_stps, [],[])
h2 = plot(x_stps, nanmean(reg2,1), '-r', 'LineWidth', 3)

stdshade(reg3,.1,'b',x_stps, [],[])
h3 = plot(x_stps, nanmean(reg3,1), '-b', 'LineWidth', 3)

stdshade(reg4,.1,'y',x_stps, [],[])
h4 = plot(x_stps, nanmean(reg4,1), '-y', 'LineWidth', 3)

stdshade(reg5,.1,'g',x_stps, [],[])
h5 = plot(x_stps, nanmean(reg5,1), '-g', 'LineWidth', 3)

%line([x_stps(1) x_stps(end)],[0.5 0.5], 'Color', 'k', 'LineWidth',3)

legend([h1, h2, h3, h4, h5],{'OFC' 'FRO' 'TEMP' 'CING' 'INS'})

title([phase ' ' lock]) 
xlabel(['delay (msec)']);
ylabel(['PTE'])
set(gca, 'FontSize', 14, 'FontWeight', 'bold')

suptitle('subtraction')