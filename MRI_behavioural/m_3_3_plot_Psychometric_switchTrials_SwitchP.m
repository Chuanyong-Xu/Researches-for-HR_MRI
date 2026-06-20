    clear;close all
    c    = [[255/255 51/255 51/255];[051/255 153/255 255/255]]; % set colors
    set(0, 'DefaultAxesFontName', 'Helvetica');   % Fonttype for axis
    set(0, 'DefaultTextFontName', 'Helvetica');   % Fonttype for text
    
    load('psychometric_fmri_exp_behvaior_29subj.mat');%******
    load('model_res_psychometric_fmri_60sub_pAlphErr1_test3.mat')
    thrs = thrs; %Difficulty levels
    isSave=1; %save figures

%% Part1: Plot group psychometric
    hit1 = mean([respVector{1,1} ./ countVector{1,1},respVector{1,2} ./ countVector{1,2}]);
    hit2 = mean([respVector{2,1} ./ countVector{2,1},respVector{2,2} ./ countVector{2,2}]);
    hit3 = mean([respVector{3,1} ./ countVector{3,1},respVector{3,2} ./ countVector{3,2}]);
    ste1 = std([respVector{1,1} ./ countVector{1,1},respVector{1,2} ./ countVector{1,2}])/sqrt(length(id));
    ste2 = std([respVector{2,1} ./ countVector{2,1},respVector{2,2} ./ countVector{2,2}])/sqrt(length(id));
    ste3 = std([respVector{3,1} ./ countVector{3,1},respVector{3,2} ./ countVector{3,2}])/sqrt(length(id));
    clear temp;
% figure layout
    figure;
    set(gcf,'units','centimeters','position',[2 5 45 18]);
%prepare data for plotting
%Inferred: subjective rules
    subjData = [model.psychometric.subj.p_anti];
    temp(:,:) = [squeeze(mean(subjData(2,1:2:end,:))),squeeze(mean(subjData(2,2:2:end,:)))];
    x_axis = model.psychometric.subj(1).x_axis;

    subplot(2,4,1); 
    plot(x_axis, temp','LineWidth',2,'Color','k'); hold on
    ylabel('Pr(Right | Rule, Anti)','fontsize',14); 
    xlabel('Probability of Anti-Clockwise','fontsize',14);
    errorbar(thrs, hit2(:,1:5), ste2(:,1:5), 'o','LineWidth',2,'MarkerFaceColor',c(1,:),'MarkerEdgeColor','k','color','k','MarkerSize',10); hold on
    errorbar(thrs, hit2(:,6:10), ste2(:,6:10), 'o','LineWidth',2,'MarkerFaceColor',c(2,:),'MarkerEdgeColor','k','color','k','MarkerSize',10)
    set(gca,'box','off', 'TitleHorizontalAlignment','center')
    axis([0,1,0,1]); set(gca,'XTiCk',0:.5:1); set(gca,'yTiCk',0:.5:1); set(gca,'XTickLabel',{'0','0.5','1'});
    set(gca,'FontSize',14,'tickdir','out')
    title('Inferred: subjective rules','fontsize',18,'FontWeight','normal');
    %Inferred: objective rules
    subjData = [model.psychometric.obj.p_anti];
    temp(:,:) = [squeeze(mean(subjData(2,1:2:end,:))),squeeze(mean(subjData(2,2:2:end,:)))];
    x_axis = model.psychometric.obj(1).x_axis;


    subplot(2,4,2); 
    plot(x_axis, temp','LineWidth',2,'Color','k');
    xlabel('Probability of Anti-Clockwise');
    hold on
    errorbar(thrs, hit3(:,1:5), ste3(:,1:5), 'o','LineWidth',2,'MarkerFaceColor',c(1,:),'MarkerEdgeColor','k','color','k','MarkerSize',10); hold on
    errorbar(thrs, hit3(:,6:10), ste3(:,6:10), 'o','LineWidth',2,'MarkerFaceColor',c(2,:),'MarkerEdgeColor','k','color','k','MarkerSize',10)
    % set(gca,'FontSize',12,'LineWidth',2.5)
    set(gca,'box','off', 'TitleHorizontalAlignment','center')
    axis([0,1,0,1]); set(gca,'XTiCk',0:.5:1); set(gca,'yTiCk',0:.5:1); set(gca,'XTickLabel',{'0','0.5','1'});
    set(gca,'FontSize',14,'tickdir','out') 
    title('Inferred: objective rules','fontsize',18,'FontWeight','normal');
    %Instructed: objective rules
    subjData = [model.psychometric.obj.p_anti];
    temp(:,:) = [squeeze(mean(subjData(1,1:2:end,:))),squeeze(mean(subjData(1,2:2:end,:)))];
    x_axis = model.psychometric.obj(1).x_axis;

    subplot(2,4,5); 
    plot(x_axis, temp','LineWidth',2,'Color','k');
    xlabel('Probability of Anti-Clockwise');  
    ylabel('Pr(Right | Rule, Anti)'); 
    title('Instructed','fontsize',14,'FontWeight','normal');hold on
    errorbar(thrs, hit1(:,1:5), ste1(:,1:5), 'o','LineWidth',2,'MarkerFaceColor',c(1,:),'MarkerEdgeColor','k','color','k','MarkerSize',10); hold on
    errorbar(thrs, hit1(:,6:10), ste1(:,6:10), 'o','LineWidth',2,'MarkerFaceColor',c(2,:),'MarkerEdgeColor','k','color','k','MarkerSize',10)
    % set(gca,'FontSize',12,'LineWidth',2.5)
    set(gca,'box','off', 'TitleHorizontalAlignment','center')
    axis([0,1,0,1]); set(gca,'XTiCk',0:.5:1); set(gca,'yTiCk',0:.5:1); set(gca,'XTickLabel',{'0','0.5','1'});
        set(gca,'FontSize',14,'tickdir','out') 
    title('Instructed','fontsize',18,'FontWeight','normal');
%% Part2: trials/time process switch to new rules
    subplot(2,4,6); 

    matdir = uigetdir('','select data folder');%select the location of time series 
    inferFile = dir([matdir,filesep,'*inferred_exp*.mat']);
    exclud1=[]%[2 3 6 10 32 33 35 37 44     9 16 21 22 23 34 39 40];
    inferFile(exclud1,:)=[];
    instrFile = dir([matdir,filesep,'*instruct_exp*.mat']);
    exclud2=[]%[2 3 6 10 32 33 35 37 44     9 16 21 22 23 34 39 40];% excluding the subjects with low accuracy or bad quality
    instrFile(exclud2,:)=[];

for id = 1:length(instrFile)
    load(fullfile(instrFile(id).folder,instrFile(id).name));
    ruleList     = exp_setting.rulelist;
    objSwitchIdx = find(xor(ruleList(1:end-1),ruleList(2:end)))+1+[0:3];
    subSwitch    = double((ruleResp==ruleList));
    
    instr_all(id,:)=[mean(double(ruleResp(objSwitchIdx(:,1)-1)==ruleList(objSwitchIdx(:,1)))),mean(subSwitch(objSwitchIdx))];
    %******
    load(fullfile(inferFile(id).folder,inferFile(id).name));
    ruleList  = exp_setting.rulelist;
    objSwitchIdx = find(xor(ruleList(1:end-1),ruleList(2:end)))+1+[0:3];
    subSwitch = double((ruleResp==ruleList));
    infer_all(id,:)=[mean(double(ruleResp(objSwitchIdx(:,1)-1)==ruleList(objSwitchIdx(:,1)))),mean(subSwitch(objSwitchIdx))];
    %******
end
    
for num =1:5 %******
    instr_ste(:,num) = nanstd(instr_all(:,num))./sqrt((length(id)-sum(isnan(instr_all(:,num)))));
    infer_ste(:,num) = nanstd(infer_all(:,num))./sqrt((length(id)-sum(isnan(infer_all(:,num)))));
end

    instr=plot([-1:3],[mean(instr_all)], '-',...
    'Color',[1 0.84 0],'LineWidth',2,'MarkerFaceColor',[1 0.84 0],'MarkerEdgeColor','k'); hold on
    errorbar([-1:3], nanmean(instr_all),instr_ste, 'ko','LineWidth',2,'MarkerFaceColor',[1 0.84 0],'MarkerSize',10)
    hold on
    
    infer=plot([-1:3],[mean(infer_all)], '-',...
    'Color',[0 0.4470 0.7410],'LineWidth',2,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerEdgeColor','k');
    errorbar([-1:3], nanmean(infer_all),infer_ste, 'ko','LineWidth',2,'MarkerFaceColor',[0 0.4470 0.7410],'MarkerSize',10)
    hold on
    
    stem(0,1,'--','Marker','none','Color','k','LineWidth',1)
    axis([-1.5 3.5 -0.05 1.05]);
    xticks([-1 0 1 2 3 4]);
    yticks([0 0.25 0.5 0.75 1]);
    % set(gca,'FontSize', 12,'LineWidth',2.5)
    set(gca,'tickdir','out')
    title('Switch to new rules','fontsize',14,'FontWeight','normal')
    legend([instr,infer],{'instr','infer'},'Location','east','fontsize',12); legend boxoff;%'location', 'SouthEast' 

    xlabel('Trials after rule switch')
    ylabel('Pr(switch to new rule)')
    set(gca,'box','off', 'TitleHorizontalAlignment','center')
    set(gca,'FontSize',14,'tickdir','out') 
    title('Switch to new rules','fontsize',18,'FontWeight','normal')

%% Part3: rule switch probability
    addpath(genpath('Hiearchica_Reasoning_model/depend/Functions'));
    load('model_res_switch_fmri_m_4_29sub.mat')
for iT = 1:3
    % notes added by xu 20231013
    %# res_o* from original data; res_m* from machine simulation; res_inp from model fitting
    %# Pr_Sw_o_5 is the switch rate from original data (5 points)
    %# Pr_Sw_m_5 is the switch rate from machine simulation data (5 points)
    %### Pr_Sw_inp_5 is the switch rate calculated from model fitting data (5 points)
    %### Inp_Pr_Sw_5 is the switch rate from model fitting data (5 points,model Optimization)
    %# M_Pr_Sw_5 is the switch rate from model fitting of machine simulation data (5 points,model Optimization, machine)
    res_o_switch_5 = [res_o{iT}.res_switch(:,1),res_o{iT}.res_switch(:,2),res_o{iT}.res_switch(:,3),res_o{iT}.res_switch(:,4),res_o{iT}.res_switch(:,5)];
    res_o_count_5 = [res_o{iT}.res_count(:,1),res_o{iT}.res_count(:,2),res_o{iT}.res_count(:,3),res_o{iT}.res_count(:,4),res_o{iT}.res_count(:,5)];
    Pr_Sw_o_5{iT} = res_o_switch_5 ./ res_o_count_5;%******
    sprintf('Check res_o{%d}: %d',iT,sum(sum(res_o{iT}.res_switch > res_o{iT}.res_count)))
    
    res_m_switch_5 = [res_m{iT}.res_switch(:,1),res_m{iT}.res_switch(:,2),res_m{iT}.res_switch(:,3),res_m{iT}.res_switch(:,4),res_m{iT}.res_switch(:,5)];
    res_m_count_5 = [res_m{iT}.res_count(:,1),res_m{iT}.res_count(:,2),res_m{iT}.res_count(:,3),res_m{iT}.res_count(:,4),res_m{iT}.res_count(:,5)];
    Pr_Sw_m_5{iT} = res_m_switch_5 ./ res_m_count_5;
    sprintf('Check res_m{%d}: %d',iT,sum(sum(res_m{iT}.res_switch > res_m{iT}.res_count)))
    
    res_inp_switch_5 = [res_inp{iT}.res_switch(:,1),res_inp{iT}.res_switch(:,2),res_inp{iT}.res_switch(:,3),res_inp{iT}.res_switch(:,4),res_inp{iT}.res_switch(:,5)]; 
    res_inp_count_5 = [res_inp{iT}.res_count(:,1),res_inp{iT}.res_count(:,2),res_inp{iT}.res_count(:,3),res_inp{iT}.res_count(:,4),res_inp{iT}.res_count(:,5)];
    Pr_Sw_inp_5{iT} = res_inp_switch_5 ./ res_inp_count_5;
    sprintf('Check res_inp{%d}: %d',iT,sum(sum(res_inp{iT}.res_switch > res_inp{iT}.res_count)))

    Inp_Pr_Sw_5{iT} = [Inp_Pr_Sw{iT}(:,1),Inp_Pr_Sw{iT}(:,2),Inp_Pr_Sw{iT}(:,3),Inp_Pr_Sw{iT}(:,4),Inp_Pr_Sw{iT}(:,5)];
    M_Pr_Sw_5{iT}   = [M_Pr_Sw{iT}(:,1),M_Pr_Sw{iT}(:,2),M_Pr_Sw{iT}(:,3),M_Pr_Sw{iT}(:,4),M_Pr_Sw{iT}(:,5)];
end

    res_pr.Pr_Sw_o_5   = Pr_Sw_o_5;
    res_pr.Pr_Sw_inp_5 = Pr_Sw_inp_5;
    res_pr.Pr_Sw_m_5   = Pr_Sw_m_5;
    res_pr.Inp_Pr_Sw_5 = Inp_Pr_Sw_5;
    res_pr.M_Pr_Sw_5   = M_Pr_Sw_5;
    
    res_pr.Pr_Sw       = Pr_Sw;
    res_pr.Inp_Sw      = Inp_Sw;
    res_pr.Inp_Pr_Sw   = Inp_Pr_Sw;
    res_pr.M_Sw        = M_Sw;
    res_pr.M_Pr_Sw     = M_Pr_Sw;

% figure_1: 5 point
    c = {[0.2784    0.4471    0.451],  [0.7412    0.3569    0.0235],[0.4660 0.6740 0.1880]};%{[0.9290 0.6940 0.1250],[0.6350 0.0780 0.1840]};

    subplot(2,4,3)
for ki = 1:3
    ste = nanstd(Pr_Sw_o_5{ki})./sqrt((length(id)-sum(isnan(Pr_Sw_o_5{ki}))));
    errorbar(1:5, nanmean(Pr_Sw_o_5{ki}),ste, 'ko','LineWidth',2,'MarkerFaceColor',c{ki},'MarkerSize',10); hold on
end
    fillsteplottusc(Pr_Sw_inp_5{1},2); hold on
    fillsteplottust(Pr_Sw_inp_5{2},2); hold on
    fillsteplotgreen(Pr_Sw_inp_5{3},2); hold on
    set(gca,'XTiCk',1:1:5);
    set(gca,'XTickLabel',{'0.1' '0.3' '0.5' '0.7' '0.9'});
    axis([0.5 5.5 0 1]); yticks([0:.25:1])
    set(gca,'tickdir','out')
    set(gca,'box','off', 'TitleHorizontalAlignment','center')
    title('Behavioral data','fontsize',12,'FontWeight','normal')
    ylabel('Pr(Switch)','fontsize',12)
    set(gca,'FontSize',14) 
    title('Behavioral data','fontsize',18,'FontWeight','normal')
    
    subplot(2,4,4)
    for ki = 1:3
        ste = nanstd(Pr_Sw_o_5{ki})./sqrt((length(id)-sum(isnan(Pr_Sw_o_5{ki}))));
        errorbar(1:5, nanmean(Pr_Sw_o_5{ki}),ste, 'ko','LineWidth',2,'MarkerFaceColor',c{ki},'MarkerSize',10)
        hold on
    end
    h3 = fillsteplottust(Inp_Pr_Sw_5{2},2); hold on
    h2 = fillsteplottusc(Inp_Pr_Sw_5{1},2); hold on
    h1 = fillsteplotgreen(Inp_Pr_Sw_5{3},2); hold on
    set(gca,'XTiCk',1:1:5);
    set(gca,'XTickLabel',{'0.1' '0.3' '0.5' '0.7' '0.9'});
    axis([0.5 5.5 0 1]); yticks([0:.25:1])
    set(gca,'tickdir','out')
    set(gca,'box','off', 'TitleHorizontalAlignment','center')
    legend([h3,h2,h1],{'E2','E1','E0'},'Position',[0.88 0.8 0.1 0.1],'Orientation','vertical'); legend boxoff;%'location', 'SouthEast' 
    set(gca,'FontSize',14) 
    title('CBM fitting behavioral data','fontsize',18,'FontWeight','normal')

    subplot(2,4,7)
    for ki = 1:3
        ste = nanstd(Pr_Sw_m_5{ki})./sqrt((length(id)-sum(isnan(Pr_Sw_m_5{ki}))));
        errorbar(1:5, nanmean(Pr_Sw_m_5{ki}),ste, 'ko','LineWidth',2,'MarkerFaceColor',c{ki},'MarkerSize',10)
        hold on
    end
    fillsteplottusc(Pr_Sw_m_5{1},2); hold on
    fillsteplottust(Pr_Sw_m_5{2},2); hold on
    fillsteplotgreen(Pr_Sw_m_5{3},2); hold on
    set(gca,'XTiCk',1:1:5);
    set(gca,'XTickLabel',{'0.1' '0.3' '0.5' '0.7' '0.9'});
    axis([0.5 5.5 0 1]); yticks([0:.25:1])
    set(gca,'tickdir','out')
    set(gca,'box','off', 'TitleHorizontalAlignment','center')
        ylabel('Pr(Switch)','fontsize',12)
        xlabel('Difficulties (Anti probability)','fontsize',12)
    set(gca,'FontSize',14) 
    title('Simulated data','fontsize',18,'FontWeight','normal')

    subplot(2,4,8)
    for ki = 1:3
        ste = nanstd(Pr_Sw_m_5{ki})./sqrt((length(id)-sum(isnan(Pr_Sw_m_5{ki}))));
        errorbar(1:5, nanmean(Pr_Sw_m_5{ki}),ste, 'ko','LineWidth',2,'MarkerFaceColor',c{ki},'MarkerSize',10)
        hold on
    end
    fillsteplottusc(M_Pr_Sw_5{1},2); hold on
    fillsteplottust(M_Pr_Sw_5{2},2); hold on
    fillsteplotgreen(M_Pr_Sw_5{3},2); hold on
    set(gca,'XTiCk',1:1:5);
    set(gca,'XTickLabel',{'0.1' '0.3' '0.5' '0.7' '0.9'});
    axis([0.5 5.5 0 1]); yticks([0:.25:1])
    set(gca,'tickdir','out')
    set(gca,'box','off', 'TitleHorizontalAlignment','center')
        xlabel('Difficulties (Anti probability)','fontsize',12)
    set(gca,'FontSize',14) 
    title('CBM fitting simulated data','fontsize',18,'FontWeight','normal')

%
if isSave
% set(gcf, 'PaperUnits', 'centimeters','PaperSize', [18 6]);
PATH = pwd; %**************************************************
print('HR_PSYCHOMETRIC_SwitchTrils_SwitchP','-dsvg','-painters')
end


%% model comparasion
% ============================================================
% compute log-likelihood, AIC, BIC

    load model_res_switch_fmri_m_4_29sub.mat
    load model_res_psychometric_fmri_60sub_pAlphErr1_test3.mat
    for subs=1:length(Input)

        ntri_subs1(subs) = sum(Input(subs).StateExp==1);
        n_total1 = ntri_subs1(subs);
        ntri_subs2(subs) = sum(Input(subs).StateExp==1 & Input(subs).Nback_act>0 & Input(subs).Nback_act<=2);
        n_total2 = ntri_subs2(subs);
        
        k_total1 = 4;
        k_total2 = 3;

        loglikelihood1(subs) = (model.psychometric.subj(subs).MLE);
        loglikelihood2(subs) = MLE(subs);
        
        loglikelihood(subs)  = (loglikelihood1(subs) .* (-1) + loglikelihood2(subs) .* (-1));

        AIC_2nd(subs) = (2*(k_total1+k_total2)-2*(loglikelihood1(subs)+loglikelihood2(subs)));
        BIC_2nd(subs) = (k_total1*log(n_total1) + k_total2*log(n_total2) -2*(loglikelihood1(subs)+loglikelihood2(subs)));

    end

    modcomp_m(1,1) = mean(loglikelihood)%CBM
    modcomp_m(1,2) = mean(AIC_2nd)%CBM
    modcomp_m(1,3) = mean(BIC_2nd)%CBM


    Sync_output = readtable('Beha_ModCom/Sync_data/Sync_output.csv','VariableNamingRule','preserve');
    Sync_output.LogLik = Sync_output.LogLik .* (-1);
    modcomp_m(2,1) = mean(Sync_output.LogLik)%Sync
    modcomp_m(2,2) = mean(Sync_output.AIC)%Sync
    modcomp_m(2,3) = mean(Sync_output.BIC)%Sync

    Hybrid_output = readtable('Beha_ModCom/Hybrid_data/Hybrid_output.csv','VariableNamingRule','preserve');
    Hybrid_output.LogLik = Hybrid_output.LogLik .* (-1);
    modcomp_m(3,1) = mean(Hybrid_output.LogLik)%Hybrid
    modcomp_m(3,2) = mean(Hybrid_output.AIC)%Hybrid
    modcomp_m(3,3) = mean(Hybrid_output.BIC)%Hybrid

    RW_output = readtable('Beha_ModCom/RW_data/RW_output.csv','VariableNamingRule','preserve');
    RW_output.LogLik = RW_output.LogLik .* (-1);
    modcomp_m(4,1) = mean(RW_output.LogLik)%RW
    modcomp_m(4,2) = mean(RW_output.AIC)%RW
    modcomp_m(4,3) = mean(RW_output.BIC)%RW

    n = 29; % number of subjects

    modcomp_e(1,1) = std(loglikelihood)/sqrt(n);
    modcomp_e(1,2) = std(AIC_2nd)/sqrt(n);
    modcomp_e(1,3) = std(BIC_2nd)/sqrt(n);
    
    modcomp_e(2,1) = std(Sync_output.LogLik)/sqrt(n);
    modcomp_e(2,2) = std(Sync_output.AIC)/sqrt(n);
    modcomp_e(2,3) = std(Sync_output.BIC)/sqrt(n);
    
    modcomp_e(3,1) = std(Hybrid_output.LogLik)/sqrt(n);
    modcomp_e(3,2) = std(Hybrid_output.AIC)/sqrt(n);
    modcomp_e(3,3) = std(Hybrid_output.BIC)/sqrt(n);
    
    modcomp_e(4,1) = std(RW_output.LogLik)/sqrt(n);
    modcomp_e(4,2) = std(RW_output.AIC)/sqrt(n);
    modcomp_e(4,3) = std(RW_output.BIC)/sqrt(n);
    
% ============================================================
% model recovery
% subject-level model ranking robustness
% Added to address:
% "it is still difficult to judge robustness of the model ranking across participants"

modelNames = {'CBM','Sync','ALR','RW'};  % Hybrid_output is labelled as ALR in the figure

% subject x model matrices
NLL_sub  = [loglikelihood(:), Sync_output.LogLik(:), Hybrid_output.LogLik(:), RW_output.LogLik(:)];
AIC_sub = [AIC_2nd(:),       Sync_output.AIC(:),    Hybrid_output.AIC(:),    RW_output.AIC(:)];
BIC_sub = [BIC_2nd(:),       Sync_output.BIC(:),    Hybrid_output.BIC(:),    RW_output.BIC(:)];

% Make sure all models have the same subject number
nSub = min([size(NLL_sub,1), size(AIC_sub,1), size(BIC_sub,1)]);
NLL_sub  = NLL_sub(1:nSub,:);
AIC_sub = AIC_sub(1:nSub,:);
BIC_sub = BIC_sub(1:nSub,:);

% Ranking:
% NLL: smaller is better
% AIC/BIC: smaller is better
[~, rank_NLL_order]  = sort(NLL_sub,  2, 'ascend');
[~, rank_AIC_order] = sort(AIC_sub, 2, 'ascend');
[~, rank_BIC_order] = sort(BIC_sub, 2, 'ascend');

best_NLL  = rank_NLL_order(:,1);
best_AIC = rank_AIC_order(:,1);
best_BIC = rank_BIC_order(:,1);

% Count how many participants each model wins
bestCount_NLL  = accumarray(best_NLL,  1, [numel(modelNames), 1])';
bestCount_AIC = accumarray(best_AIC, 1, [numel(modelNames), 1])';
bestCount_BIC = accumarray(best_BIC, 1, [numel(modelNames), 1])';

bestProp_NLL  = bestCount_NLL  ./ nSub;
bestProp_AIC = bestCount_AIC ./ nSub;
bestProp_BIC = bestCount_BIC ./ nSub;

rankSummary = table(modelNames(:), ...
    bestCount_NLL(:),  bestProp_NLL(:), ...
    bestCount_AIC(:), bestProp_AIC(:), ...
    bestCount_BIC(:), bestProp_BIC(:), ...
    'VariableNames', {'Model', ...
    'BestN_NLL','BestProp_NLL', ...
    'BestN_AIC','BestProp_AIC', ...
    'BestN_BIC','BestProp_BIC'});

disp('Subject-level model ranking robustness:')
disp(rankSummary)

% Paired subject-level comparison: CBM vs alternative models
% For NLL: CBM - alternative > 0 means CBM is better
% For AIC/BIC: alternative - CBM > 0 means CBM is better
compNames = {'CBM_vs_Sync','CBM_vs_ALR','CBM_vs_RW'}';

p_NLL  = nan(3,1);
p_AIC = nan(3,1);
p_BIC = nan(3,1);

medianDiff_NLL  = nan(3,1);
medianDiff_AIC = nan(3,1);
medianDiff_BIC = nan(3,1);

for mi = 2:4
    idx = mi - 1;

    diff_NLL  = NLL_sub(:,mi) - NLL_sub(:,1);
    diff_AIC = AIC_sub(:,mi) - AIC_sub(:,1);
    diff_BIC = BIC_sub(:,mi) - BIC_sub(:,1);

    medianDiff_NLL(idx)  = median(diff_NLL,  'omitnan');
    medianDiff_AIC(idx) = median(diff_AIC, 'omitnan');
    medianDiff_BIC(idx) = median(diff_BIC, 'omitnan');
    p_NLL(idx)  = signrank(diff_NLL,  0, 'tail', 'right');
    p_AIC(idx) = signrank(diff_AIC, 0, 'tail', 'right');
    p_BIC(idx) = signrank(diff_BIC, 0, 'tail', 'right');
end

pairedSummary = table(compNames, ...
    medianDiff_NLL,  p_NLL, ...
    medianDiff_AIC, p_AIC, ...
    medianDiff_BIC, p_BIC, ...
    'VariableNames', {'Comparison', ...
    'MedianDelta_NLL','P_NLL', ...
    'MedianDelta_AIC','P_AIC', ...
    'MedianDelta_BIC','P_BIC'});

disp('Paired subject-level model comparison:')
disp(pairedSummary)

% Save numerical results for reporting
writetable(rankSummary,  'HR_subjectLevel_modelRanking_summary.csv');
writetable(pairedSummary,'HR_subjectLevel_modelRanking_pairedTests.csv');





%% ------------------------------------------------------------
% plot model compare, model recovery, parameter recovery, leave-one-out
% prediction
%------------------------------------------------------------
% plot
addpath(genpath('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/matlab-special-heatmap'));

figure
set(gcf,'unit','pixels','position',[300 50 1350 760]);
set(gcf,'Color','w');

% unified cold color palette
col_NLL   = [0.40 0.56 0.78];
col_AIC  = [0.40 0.74 0.84];
col_BIC  = [0.58 0.82 0.78];

col_CBM  = [0.28 0.52 0.78];
col_Sync = [0.42 0.72 0.86];
col_ALR  = [0.55 0.78 0.76];
col_RW   = [0.62 0.64 0.78];

col_1B   = [0.28 0.52 0.78];
col_2B   = [0.55 0.78 0.76];

col_scatter = [0.32 0.60 0.78];
col_diag    = [0.25 0.25 0.25];

% ============================================================
% (A) Model comparison
% ============================================================
subplot('position',[0.06 0.58 0.50 0.34])
set(gca,'Fontsize',18);

hbar = bar(modcomp_m);

hbar(1).FaceColor = col_NLL;   hbar(1).EdgeColor = col_NLL;
hbar(2).FaceColor = col_AIC;  hbar(2).EdgeColor = col_AIC;
hbar(3).FaceColor = col_BIC;  hbar(3).EdgeColor = col_BIC;

hbar(1).FaceAlpha = 0.55; hbar(2).FaceAlpha = 0.55; hbar(3).FaceAlpha = 0.55;
hbar(1).EdgeAlpha = 0.85; hbar(2).EdgeAlpha = 0.85; hbar(3).EdgeAlpha = 0.85;
hbar(1).LineWidth = 2;    hbar(2).LineWidth = 2;    hbar(3).LineWidth = 2;

ylabel('Values');
xlabel('Models');
xticklabels({'CBM','Sync','ALR','RW'});

hold on;
ngroups = size(modcomp_m,1);
nbars   = size(modcomp_m,2);
groupwidth = min(0.8, nbars/(nbars+1.5));

for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);

    errorbar(x, modcomp_m(:,i), modcomp_e(:,i), ...
        'k', 'linestyle','none','LineWidth',1.5);
end

legend([hbar(1), hbar(2), hbar(3)],{'NLL','AIC','BIC'},'Location','northwest','fontsize',20);
legend boxoff
title('Model comparison','fontsize',20,'FontWeight','normal',...
    'Units','normalized','Position',[0.5 1.08 0]);

text(-0.10, 1.2, 'A','Units','normalized', ...
    'FontSize',22, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');

set(gca,'TickDir','out','LineWidth',1.2);
box off; grid off;


% ============================================================
% (B) Subject-level winning model counts
% ============================================================
subplot('position',[0.66 0.58 0.27 0.34])

rankMat = [bestCount_NLL; bestCount_AIC; bestCount_BIC];

modelColors = [
    col_CBM;
    col_Sync;
    col_ALR;
    col_RW
];

h = bar(rankMat, 'stacked');

for ii = 1:numel(h)
    h(ii).FaceColor = modelColors(ii,:);
    h(ii).EdgeColor = modelColors(ii,:);
    h(ii).LineWidth = 2;
    h(ii).FaceAlpha = 0.55;
    h(ii).EdgeAlpha = 0.85;
end

set(gca,'XTickLabel',{'NLL','AIC','BIC'}, ...
    'FontSize',18, ...
    'TickDir','out', ...
    'LineWidth',1.2);

ylabel('Number of subjects');
ylim([0 30]);

legend(modelNames, 'Location','eastoutside', 'Orientation','vertical');
legend boxoff

title('Subject-level winning model counts','fontsize',20,'FontWeight','normal',...
    'Units','normalized','Position',[0.5 1.08 0]);
box off

text(-0.22, 1.2, 'B','Units','normalized', ...
    'FontSize',22, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');


% ============================================================
% load parameter recovery and CV results
% ============================================================
load model_res_switch_parameter_recovery_allsubj.mat
load model_res_switch_CV.mat


% ============================================================
% (C) Parameter recovery: sigma
% ============================================================
subplot('position',[0.06 0.10 0.18 0.32])
hold on;

x = true_sigma(:);
y = rec_sigma(:);
valid = ~isnan(x) & ~isnan(y);

[r_sigma, p_sigma] = corr(x(valid), y(valid), 'Type', 'Pearson');

scatter(x(valid), y(valid), 28, ...
    'filled', ...
    'MarkerFaceColor', col_scatter, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.55);

plot([0.1 10], [0.1 10], '--', ...
    'Color', col_diag, ...
    'LineWidth', 1.5);

xlabel('True \sigma');
ylabel('Recovered \sigma');
title(sprintf('\\sigma, r = %.2f, p = %.2g', r_sigma, p_sigma), ...
    'FontSize',18,'FontWeight','normal',...
    'Units','normalized','Position',[0.5 1.08 0]);

xlim([0 10]);
ylim([0 10]);
axis square
box off

set(gca,'FontSize',18,'TickDir','out','LineWidth',1.2);

text(-0.26, 1.25, 'C','Units','normalized', ...
    'FontSize',22, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');


% ============================================================
% (D) Parameter recovery: alpha
% ============================================================
subplot('position',[0.29 0.10 0.18 0.32])
hold on;

x = true_alpha(:);
y = rec_alpha(:);
valid = ~isnan(x) & ~isnan(y);

[r_alpha, p_alpha] = corr(x(valid), y(valid), 'Type', 'Pearson');

scatter(x(valid), y(valid), 28, ...
    'filled', ...
    'MarkerFaceColor', col_scatter, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.55);

plot([0 1], [0 1], '--', ...
    'Color', col_diag, ...
    'LineWidth', 1.5);

xlabel('True \alpha');
ylabel('Recovered \alpha');
title(sprintf('\\alpha, r = %.2f, p = %.2g', r_alpha, p_alpha), ...
    'FontSize',18,'FontWeight','normal',...
    'Units','normalized','Position',[0.5 1.08 0]);

xlim([0 1]);
ylim([0 1]);
axis square
box off

set(gca,'FontSize',18,'TickDir','out','LineWidth',1.2);

text(-0.2, 1.25, 'D','Units','normalized', ...
    'FontSize',22, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');


% ============================================================
% (E) Parameter recovery: pam3
% ============================================================
subplot('position',[0.52 0.10 0.18 0.32])
hold on;

x = true_pam3(:);
y = rec_pam3(:);
valid = ~isnan(x) & ~isnan(y);

[r_pam3, p_pam3] = corr(x(valid), y(valid), 'Type', 'Pearson');

scatter(x(valid), y(valid), 28, ...
    'filled', ...
    'MarkerFaceColor', col_scatter, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceAlpha', 0.55);

plot([0 1], [0 1], '--', ...
    'Color', col_diag, ...
    'LineWidth', 1.5);

xlabel('True \lambda');
ylabel('Recovered \lambda');
title(sprintf('pam3, r = %.2f, p = %.2g', r_pam3, p_pam3), ...
    'FontSize',18,'FontWeight','normal',...
    'Units','normalized','Position',[0.5 1.08 0]);

xlim([0 1]);
ylim([0 1]);
axis square
box off

set(gca,'FontSize',18,'TickDir','out','LineWidth',1.2);

text(-0.2, 1.25, 'E','Units','normalized', ...
    'FontSize',22, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');


% ============================================================
% (F) 2-fold CV: CBM prediction
% ============================================================
subplot('position',[0.76 0.10 0.20 0.32])
hold on;

CV_obs_subj   = squeeze(nanmean(CV_Pr_Sw_obs, 2));     % subj x T x difficulty
CV_model_subj = squeeze(nanmean(CV_Pr_Sw_model, 2));   % subj x T x difficulty

xvalue = 1:5;
clear h_line h_fill

for iT = 1:2

    y_obs = squeeze(CV_obs_subj(:,iT,:));
    y_mod = squeeze(CV_model_subj(:,iT,:));

    obs_mean = nanmean(y_obs, 1);
    obs_sem  = nanstd(y_obs, [], 1) ./ sqrt(sum(~isnan(y_obs), 1));

    mod_mean = nanmean(y_mod, 1);
    mod_sem  = nanstd(y_mod, [], 1) ./ sqrt(sum(~isnan(y_mod), 1));

    if iT == 1
        thisColor = col_1B;
    else
        thisColor = col_2B;
    end

    % observed held-out data
    errorbar(xvalue, obs_mean, obs_sem, 'o', ...
        'LineWidth', 1.8, ...
        'Color', [0.15 0.15 0.15], ...
        'MarkerFaceColor', thisColor, ...
        'MarkerEdgeColor', [0.15 0.15 0.15], ...
        'MarkerSize', 7);
    hold on;

    % model held-out prediction
    xx = [xvalue, fliplr(xvalue)];
    yy = [mod_mean + mod_sem, fliplr(mod_mean - mod_sem)];

    h_fill(iT) = fill(xx, yy, thisColor, ...
        'FaceAlpha', 0.22, ...
        'EdgeColor', 'none');
    hold on;

    h_line(iT) = plot(xvalue, mod_mean, '-', ...
        'Color', thisColor, ...
        'LineWidth', 2.5);
    hold on;

end

set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.1','0.3','0.5','0.7','0.9'});
axis([0.5 5.5 0 1]);
yticks(0:0.25:1);

xlabel('Difficulties');
ylabel('Pr(Switch)');
title('2-fold CV: CBM prediction','fontsize',18,'FontWeight','normal',...
    'Units','normalized','Position',[0.5 1.08 0]);

legend([h_line(1), h_line(2)], {'E1','E2'}, ...
    'Location','best');
legend boxoff

set(gca,'FontSize',18,'TickDir','out','LineWidth',1.2);
box off

text(-0.2, 1.2, 'F','Units','normalized', ...
    'FontSize',22, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');


% ============================================================
% save
% ============================================================
print('HR_modelComparison_ranking_recovery_CV','-dsvg','-painters')











%% ============================================================
% Minimal BioRender-style cartoon (no custom functions)
% x-axis = perceptual difficulty
% y-axis = errors
% ============================================================
clear; close all; clc;

% -------------------- style --------------------
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');
set(0,'DefaultAxesFontSize',11);
set(0,'DefaultTextFontSize',11);

fig = figure('Color','w','Position',[100 100 700 350]);

green       = [92 137 112]/255;
orange      = [191 103 72]/255;
lightGreen  = [230 230 230]/255;
lightOrange = [245 245 245]/255;
gray        = [0.65 0.65 0.65];
gray_edge   = [0.85 0.85 0.85];

lims = [-0.45 2.45];
b = 2.2;                            % decision boundary: D + E = b
% 
[xbg,ybg] = meshgrid(linspace(lims(1),lims(2),300), ...
                     linspace(lims(1),lims(2),300));
[Dg,Eg] = meshgrid(0:2,0:2);
D = Dg(:);   % x = difficulty
E = Eg(:);   % y = errors

% ---- define the plane / boundary by two endpoints ----
x1 = lims(1);   y1 = 2.20;   % left endpoint: move DOWN a bit
x2 = lims(2);   y2 = -0.2;   % right endpoint: move UP a bit
% line equation: y = m*x + c
m = (y2 - y1) / (x2 - x1);
c = y1 - m * x1;
% line for plotting
xline = linspace(lims(1), lims(2), 200);
yline = m * xline + c;
% background classification
zbg = ybg - (m * xbg + c);
% point classification
isSwitch = E > (m * D + c);
% ============================================================
% A. Near-orthogonal geometry
% ============================================================
axA = axes('Position',[0.06 0.22 0.27 0.63]); hold(axA,'on');

ha = pcolor(axA,xbg,ybg,double(zbg>0));
set(ha,'EdgeColor','none','FaceAlpha',0.3);
shading(axA,'flat'); 
colormap([lightGreen; lightOrange]);
plot(axA,xline,yline,'k-','LineWidth',2);

scatter(axA,D(~isSwitch),E(~isSwitch),250,green,...
    'filled','MarkerEdgeColor',gray_edge,'LineWidth',2);
scatter(axA,D(isSwitch),E(isSwitch),250,orange,...
    'filled','MarkerEdgeColor',gray_edge,'LineWidth',2);

axis(axA,[lims lims]); axis(axA,'square'); box(axA,'off');
set(axA,'XTick',0:2,'XTickLabel',{'High','Medium','Low'},...
        'YTick',0:2,'YTickLabel',{'0E','1E','2E'},...
        'TickDir','out','LineWidth',1.0);

title(axA,'Near-orthogonal geometry','FontSize',18,'FontWeight','bold');
xlabel(axA,'Perceptual difficulty','Color',green,'FontSize',15,'FontWeight','bold');
ylabel(axA,'Errors','Color',orange,'FontSize',15,'FontWeight','bold');

text(axA,-0.25,-0.2,'No switch','Color',green,'FontSize',15,'FontWeight','bold');
text(axA,1.70,2.2,'Switch','Color',orange,'FontSize',15,'FontWeight','bold');

plot(axA,[1.03 1.38],[1.03 1.38],'k--','LineWidth',1.1);
quiver(axA,1.20,1.20, 0.16, 0.16,0,'Color','k','LineWidth',1.1,'MaxHeadSize',1.2);
quiver(axA,1.20,1.20,-0.16,-0.16,0,'Color','k','LineWidth',1.1,'MaxHeadSize',1.2);
text(axA,1.52,1.35,{'large','margin'},'Color',gray,'FontSize',16);

annotation(fig,'textbox',[0.02 0.92 0.04 0.06],...
    'String','A','EdgeColor','none','FontSize',30,'FontWeight','bold');
annotation(fig,'textbox',[0.05 0.10 0.30 0.05],...
    'String','Orthogonal axes support reliable separation.',...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',14,'Color',gray);

% ============================================================
% B. Simple readout
% ============================================================
annotation(fig,'textbox',[0.39 0.92 0.04 0.06],...
    'String','B','EdgeColor','none','FontSize',30,'FontWeight','bold');
% annotation(fig,'textbox',[0.45 0.85 0.17 0.05],...
%     'String','Simple readout','EdgeColor','none',...
%     'HorizontalAlignment','center','FontSize',18,'FontWeight','bold');


% Flow chart
annotation(fig,'textbox',[0.38 0.70 0.1 0.09],...
    'String',{'D','(difficulty)'},...
    'EdgeColor',green,'Color',green,'BackgroundColor',[0.97 1.00 0.97],...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'FontSize',12,'FontWeight','bold');

annotation(fig,'textbox',[0.52 0.70 0.1 0.09],...
    'String',{'E','(errors)'},...
    'EdgeColor',orange,'Color',orange,'BackgroundColor',[1.00 0.97 0.94],...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'FontSize',12,'FontWeight','bold');

annotation(fig,'textbox',[0.42 0.57 0.15 0.06],...
    'String','Linear readout',...
    'EdgeColor',[0.6 0.6 0.6],'Color',[0.2 0.2 0.2],...
    'BackgroundColor',[0.96 0.96 0.96],...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'FontSize',12,'FontWeight','bold');

annotation(fig,'textbox',[0.39 0.46 0.22 0.06],...
    'String','z = w_D D + w_E E + b',...
    'EdgeColor',[0.3 0.3 0.3],'BackgroundColor','w',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'FontSize',12,'FontWeight','bold');

annotation(fig,'textbox',[0.39 0.35 0.22 0.06],...
    'String','P(switch) = sigmoid(z)',...
    'EdgeColor',[0.3 0.3 0.3],'BackgroundColor','w',...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'FontSize',12,'FontWeight','bold');

annotation(fig,'textbox',[0.38 0.20 0.11 0.05],...
    'String','No switch',...
    'EdgeColor',lightGreen,'Color','k','BackgroundColor',lightGreen,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'FontSize',12,'FontWeight','bold');

annotation(fig,'textbox',[0.52 0.20 0.11 0.05],...
    'String','Switch',...
    'EdgeColor',lightOrange,'Color','k','BackgroundColor',lightOrange,...
    'HorizontalAlignment','center','VerticalAlignment','middle',...
    'FontSize',12,'FontWeight','bold');

annotation(fig,'arrow',[0.42 0.465],[0.70 0.63],'Color',green,'LineWidth',1.2);
annotation(fig,'arrow',[0.57 0.525],[0.70 0.63],'Color',orange,'LineWidth',1.2);
annotation(fig,'arrow',[0.49 0.49],[0.57 0.53],'Color',gray,'LineWidth',1.2);
annotation(fig,'arrow',[0.49 0.49],[0.46 0.41],'Color',gray,'LineWidth',1.2);
annotation(fig,'arrow',[0.49 0.425],[0.35 0.25],'Color',gray,'LineWidth',1.2);
annotation(fig,'arrow',[0.49 0.565],[0.35 0.25],'Color',gray,'LineWidth',1.2);

annotation(fig,'textbox',[0.35 0.10 0.31 0.06],...
    'String','Independent weighting enables simple switch inference.',...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',14,'Color',gray);

% ============================================================
% C. Compressed geometry
% ============================================================
theta = pi/4;
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
S = [1.15 0; 0 0.18];

pts = [D E];
ctr = mean(pts,1);
ptsC = (R * S * (pts - ctr)')' + ctr;

rng(1);
DC = ptsC(:,1) + 0.03*randn(size(D));
EC = ptsC(:,2) + 0.03*randn(size(E));
isSwitchC = (DC + EC) > b;

axC = axes('Position',[0.72 0.22 0.25 0.63]); hold(axC,'on');

hc = pcolor(axC,xbg,ybg,double(zbg>0));
set(hc,'EdgeColor','none','FaceAlpha',0.3);
shading(axA,'flat'); 
colormap([lightGreen; lightOrange]);
plot(axC,xline,yline,'k-','LineWidth',2);

scatter(axC,DC(~isSwitchC),EC(~isSwitchC),250,green,...
    'filled','MarkerEdgeColor',gray_edge,'LineWidth',2);
scatter(axC,DC(isSwitchC),EC(isSwitchC),250,orange,...
    'filled','MarkerEdgeColor',gray_edge,'LineWidth',2);

axis(axC,[lims lims]); axis(axC,'square'); box(axC,'off');
set(axC,'XTick',0:2,'XTickLabel',{'High','Medium','Low'},...
        'YTick',0:2,'YTickLabel',{'0E','1E','2E'},...
        'TickDir','out','LineWidth',1.0);

title(axC,'Compressed geometry','FontSize',18,'FontWeight','bold');
xlabel(axC,'Perceptual difficulty','Color',green,'FontSize',15,'FontWeight','bold');
ylabel(axC,'Errors','Color',orange,'FontSize',15,'FontWeight','bold');

text(axC,-0.25,-0.2,'No switch','Color',green,'FontSize',15,'FontWeight','bold');
text(axC,1.70,2.2,'Switch','Color',orange,'FontSize',15,'FontWeight','bold');

plot(axC,[1.22 1.42],[1.22 1.42],'k--','LineWidth',1.1);
quiver(axC,1.32,1.32,0.10,0.10,0,'Color','k','LineWidth',1.1,'MaxHeadSize',1.2);
text(axC,1.63,1.30,{'small','margin'},'Color',gray,'FontSize',14);

annotation(fig,'textbox',[0.69 0.92 0.04 0.06],...
    'String','C','EdgeColor','none','FontSize',30,'FontWeight','bold');
annotation(fig,'textbox',[0.70 0.10 0.28 0.05],...
    'String','Compression reduces separability and reliability.',...
    'EdgeColor','none','HorizontalAlignment','center',...
    'FontSize',14,'Color',gray);

% -------------------- save --------------------
outName = 'minimal_biorender_geometry_readout_noFunction';

print(outName, '-dsvg','-painters');














%% Plot the behaviour for each alternative models (only for rule, not consider perceptual accuracy)
% matdir = uigetdir('','select data folder');%select the location of time series 
% dir1 = dir([matdir,filesep,'*Behavioral_data_subject*.csv']);
% load model_res_switch_fmri_m_4_29sub.mat
% for subs =1:length(dir1)
%     subs_data = readtable([matdir,filesep,dir1(subs).name]);
% 
%     pr_of_swi{subs,1}  = double(Input(subs).tDev( Input(subs).StateExp==1) )';
%     pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.9)=5; pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.1)=1;
%     pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.7)=4; pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.3)=2;
%     pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.5)=3;
% 
%     pr_of_swi{subs,2}  = double(Input(subs).Nback_act( Input(subs).StateExp==1))';
%     pr_of_swi{subs,3}  = 0;
%     pr_of_swi{subs,4}  = double(subs_data.Response_likelihood)';
%     pr_of_swi{subs,5}  = double(subs_data.Reward)';
%     % accuracy
%         pr_of_swi{subs,4}(pr_of_swi{subs,5}==0) = 1-pr_of_swi{subs,4}(pr_of_swi{subs,5}==0);
%         %sync model
%         Pr_Sw_inp_3_sync{1}(subs,1) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==1));
%         Pr_Sw_inp_3_sync{1}(subs,2) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==2));
%         Pr_Sw_inp_3_sync{1}(subs,3) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==3));
%         Pr_Sw_inp_3_sync{1}(subs,4) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==4));
%         Pr_Sw_inp_3_sync{1}(subs,5) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==5));
%         %subs data
%         Pr_Sw_inp_3_sync{2}(subs,1) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==1)/sum(pr_of_swi{subs,1}==1);
%         Pr_Sw_inp_3_sync{2}(subs,2) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==2)/sum(pr_of_swi{subs,1}==2);
%         Pr_Sw_inp_3_sync{2}(subs,3) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==3)/sum(pr_of_swi{subs,1}==3);
%         Pr_Sw_inp_3_sync{2}(subs,4) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==4)/sum(pr_of_swi{subs,1}==4);
%         Pr_Sw_inp_3_sync{2}(subs,5) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==5)/sum(pr_of_swi{subs,1}==5);
%         %Bayesian simulation data
%         Pr_Sw_inp_3_sync{3}(subs,1) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.1))/sum(MachineSimulation(subs).tDev==0.1 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.1));
%         Pr_Sw_inp_3_sync{3}(subs,2) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.3))/sum(MachineSimulation(subs).tDev==0.3 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.3));
%         Pr_Sw_inp_3_sync{3}(subs,3) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.5))/sum(MachineSimulation(subs).tDev==0.5 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.5));
%         Pr_Sw_inp_3_sync{3}(subs,4) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.7))/sum(MachineSimulation(subs).tDev==0.7 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.7));
%         Pr_Sw_inp_3_sync{3}(subs,5) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.9))/sum(MachineSimulation(subs).tDev==0.9 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.9));
% 
% end
% 
% set(0, 'DefaultAxesFontName', 'Helvetica');   % Fonttype for axis
% set(0, 'DefaultTextFontName', 'Helvetica');   % Fonttype for text
% addpath(genpath('/Users/vincentxu/Desktop/MRI_behavioural/Hiearchica_Reasoning_model/depend/Functions'));
% load('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/model_res_switch_fmri_m_4_29sub.mat')
% %
% rn=[];
% figure
% set(gcf,'unit','pixels','position',[200 50 450 400]);
% subplot(1,1,1)
%     set(gca,'FontSize',14,'tickdir','out')
% 
% fillsteplotgreen(Pr_Sw_inp_3_sync{1},4); hold on
% fillsteplotblue(Pr_Sw_inp_3_sync{2},4); hold on
% fillsteplotred(Pr_Sw_inp_3_sync{3},4); hold on
% 
% title('Model Comparison','fontsize',18,'FontWeight','Normal')
% 
% % print('Pr_switch_bayesian_rnn0', '-dsvg','-painters')
