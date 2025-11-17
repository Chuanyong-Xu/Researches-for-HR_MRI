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
    xlabel('Probability of Anti-Colokwise','fontsize',14);
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
    xlabel('Probability of Anti-Colokwise');  
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
    xlabel('Probability of Anti-Colokwise');  
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

%% logistic regression (pr of only 1-er, 2-er)
Sw0   = reshape(Pr_Sw_o_5{1},[],1);
Sw1   = reshape(Pr_Sw_o_5{2},[],1);
Sw2   = reshape(Pr_Sw_o_5{3},[],1);
Sw    = [Sw0; Sw1; Sw2];
Ners  = [repelem([1 2 0]', 29*5)];
Diffs = [repelem([3 2 1 2 3]',29); repelem([3 2 1 2 3]',29); repelem([3 2 1 2 3]',29)];

subID = repmat(1:29,1,5*3)';
trianID = [reshape(res_o{1}.res_count,[],1);reshape(res_o{2}.res_count,[],1);reshape(res_o{3}.res_count,[],1)];

Pr_Subs = table(subID,Ners,Diffs,Sw,trianID);% make table
Pr_Subs.Ners(Pr_Subs.Ners==0)=nan;

Pr_Subs = rmmissing(Pr_Subs);
% Pr_Subs.Ners = categorical(Pr_Subs.Ners);
% Pr_Subs.Diffs = categorical(Pr_Subs.Diffs);

glme = fitglme(Pr_Subs, 'Sw ~ 1 + Diffs + Ners + (1 + Diffs + Ners|subID)', ...
  'Distribution','binomial','Weights',[Pr_Subs.trianID],'Link','logit','FitMethod','Laplace');

disp(glme)
[beta, betaCI, stats] = fixedEffects(glme);
exp(beta)





%% logistic regression (binary)
% % % clear; clc
% % % %load
% % % load('model_res_switch_fmri_m_4_29sub.mat');
% % % %
% % % beha_data=Input;
% % % for subs = 1:length(beha_data)
% % %     infer_data(subs).ID = repelem(subs, sum(beha_data(subs).StateExp==1))';
% % %     infer_data(subs).tDev = beha_data(subs).tDev(beha_data(subs).StateExp==1);
% % %         infer_data(subs).tDev(infer_data(subs).tDev==0.1)=3;
% % %         infer_data(subs).tDev(infer_data(subs).tDev==0.3)=2;
% % %         infer_data(subs).tDev(infer_data(subs).tDev==0.5)=1;
% % %         infer_data(subs).tDev(infer_data(subs).tDev==0.7)=2;
% % %         infer_data(subs).tDev(infer_data(subs).tDev==0.9)=3;
% % %     infer_data(subs).SW = beha_data(subs).SW(beha_data(subs).StateExp==1);
% % %     infer_data(subs).Nback = beha_data(subs).Nback_act(beha_data(subs).StateExp==1);
% % % end
% % % 
% % % pr_cat = [vertcat(infer_data.ID), vertcat(infer_data.Nback),...
% % %     vertcat(infer_data.tDev), vertcat(infer_data.SW)]; % 
% % % 
% % % % pr_cat(pr_cat(:,2)==0)=nan;
% % % pr_cat(pr_cat(:,2)>2)=nan;
% % % pr_cat = rmmissing(pr_cat);
% % % % pr_cat(:,2:3) = zscore(pr_cat(:,2:3));
% % % 
% % % Pr_Subs = array2table(pr_cat,'VariableNames', {'subID', 'Ners','Diffs', 'Sw'});% make table
% % % % Pr_Subs.Diffs = categorical(Pr_Subs.Diffs);
% % % % Pr_Subs.Ners = categorical(Pr_Subs.Ners);
% % % % Pr_Subs.subID = categorical(Pr_Subs.subID);
% % % 
% % % glme = fitglme(Pr_Subs, 'Sw ~ 1 + Ners + Diffs + (1 + Ners + Diffs|subID)', ...
% % %   'Distribution','binomial','Link','logit','FitMethod','Laplace');
% % % 
% % % disp(glme)
% % % [beta, betaCI, stats] = fixedEffects(glme);
% % % exp(beta)
% % % 


%% model comparasion
% compute log-likelihood, AIC, BIC
%2nd level
    % % % load model_res_switch_fmri_m_4_29sub.mat
    % % % for subs=1:length(Input)
    % % %     ntri_subs(subs) = sum(Input(subs).StateExp==1);
    % % %     k_total = 3;
    % % %     n_total = ntri_subs(subs);
    % % %     loglikelihood(subs) = MLE(subs);
    % % %     AIC_2nd(subs) = (2*k_total-2*loglikelihood(subs));
    % % %     BIC_2nd(subs) = (k_total*log(n_total)-2*loglikelihood(subs));
    % % % end
    % % % modcomp_m(1,1) = mean(loglikelihood)%CBM
    % % % modcomp_m(1,2) = mean(AIC_2nd)%CBM
    % % % modcomp_m(1,3) = mean(BIC_2nd)%CBM
    % % % 
    % % % %1st level
    % % % load model_res_psychometric_fmri_60sub_pAlphErr1_test3.mat
    % % % for subs=1:length(Input)
    % % %     ntri_subs(subs) = sum(Input(subs).StateExp==1);
    % % %     k_total = 6;
    % % %     n_total = ntri_subs(subs);
    % % %     loglikelihood(subs) = (model.psychometric.subj(subs).MLE);
    % % %     AIC_1nd(subs) = (2*k_total-2*loglikelihood(subs));
    % % %     BIC_1nd(subs) = (k_total*log(n_total)-2*loglikelihood(subs));
    % % % end
    % % % loglikelihood_m(1) = mean(loglikelihood)
    % % % AIC_1nd_m(1) = mean(AIC_1nd)
    % % % BIC_1nd_m(1) = mean(BIC_1nd)

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
        
        loglikelihood(subs)  = loglikelihood1(subs)+loglikelihood2(subs);

        AIC_2nd(subs) = (2*(k_total1+k_total2)-2*(loglikelihood1(subs)+loglikelihood2(subs)));
        BIC_2nd(subs) = (k_total1*log(n_total1) + k_total2*log(n_total2) -2*(loglikelihood1(subs)+loglikelihood2(subs)));

    end
    modcomp_m(1,1) = mean(loglikelihood)%CBM
    modcomp_m(1,2) = mean(AIC_2nd)%CBM
    modcomp_m(1,3) = mean(BIC_2nd)%CBM


    Sync_output = readtable('Beha_ModCom/Sync_data/Sync_output.csv','VariableNamingRule','preserve');
    modcomp_m(2,1) = mean(Sync_output.LogLik)%Sync
    modcomp_m(2,2) = mean(Sync_output.AIC)%Sync
    modcomp_m(2,3) = mean(Sync_output.BIC)%Sync

    Hybrid_output = readtable('Beha_ModCom/Hybrid_data/Hybrid_output.csv','VariableNamingRule','preserve');
    modcomp_m(3,1) = mean(Hybrid_output.LogLik)%Hybrid
    modcomp_m(3,2) = mean(Hybrid_output.AIC)%Hybrid
    modcomp_m(3,3) = mean(Hybrid_output.BIC)%Hybrid

    RW_output = readtable('Beha_ModCom/RW_data/RW_output.csv','VariableNamingRule','preserve');
    modcomp_m(4,1) = mean(RW_output.LogLik)%RW
    modcomp_m(4,2) = mean(RW_output.AIC)%RW
    modcomp_m(4,3) = mean(RW_output.BIC)%RW
    
% plot
addpath(genpath('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/matlab-special-heatmap'));
    figure
    set(gcf,'unit','pixels','position',[500 50 550 300]);
    set(gca,'Fontsize',14);

    % scatter(datasample((1:1:29)*0.005+0.7, 29), loglikelihood, 10, colorlist(50,:), 'filled');hold on
    % scatter(datasample((1:1:29)*0.005+0.95, 29), AIC_2nd, 10, colorlist(100,:), 'filled');hold on
    % scatter(datasample((1:1:29)*0.005+1.15, 29), BIC_2nd, 10, colorlist(150,:), 'filled');hold on
    % legend off
    % hold on

    hbar            = bar(modcomp_m);
    colorlist       = flipud(slanCM(21)); 
    hbar(1).FaceColor  = colorlist(50,:);  hbar(1).EdgeColor  = colorlist(50,:); 
    hbar(2).FaceColor  = colorlist(100,:); hbar(2).EdgeColor  = colorlist(100,:);
    hbar(3).FaceColor  = colorlist(150,:); hbar(3).EdgeColor  = colorlist(150,:);
    hbar(1).FaceAlpha=0.5; hbar(2).FaceAlpha=0.5; hbar(3).FaceAlpha=0.5;
    hbar(1).EdgeAlpha=0.5; hbar(2).EdgeAlpha=0.5; hbar(3).EdgeAlpha=0.5;
    hbar(1).LineWidth=2; hbar(2).LineWidth=2; hbar(3).LineWidth=2
    ylabel('Values'); xlabel('Models');
    xticklabels({'CBM','Sync','ALR','RW'})
    legend([hbar(1), hbar(2), hbar(3)],{'LL','AIC','BIC'},'Location','northwest')
    title('Model comparison','fontsize',18,'FontWeight','normal')
    box off;grid off

    print('HR_modelComparasions','-dsvg','-painters')




    
    %% Plot the behaviour for each alternative models (only for rule, not consider perceptual accuracy)
matdir = uigetdir('','select data folder');%select the location of time series 
dir1 = dir([matdir,filesep,'*Behavioral_data_subject*.csv']);
load model_res_switch_fmri_m_4_29sub.mat
for subs =1:length(dir1)
    subs_data = readtable([matdir,filesep,dir1(subs).name]);

    pr_of_swi{subs,1}  = double(Input(subs).tDev( Input(subs).StateExp==1) )';
    pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.9)=5; pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.1)=1;
    pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.7)=4; pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.3)=2;
    pr_of_swi{subs,1}(pr_of_swi{subs,1}==0.5)=3;

    pr_of_swi{subs,2}  = double(Input(subs).Nback_act( Input(subs).StateExp==1))';
    pr_of_swi{subs,3}  = 0;
    pr_of_swi{subs,4}  = double(subs_data.Response_likelihood)';
    pr_of_swi{subs,5}  = double(subs_data.Reward)';
    % accuracy
        pr_of_swi{subs,4}(pr_of_swi{subs,5}==0) = 1-pr_of_swi{subs,4}(pr_of_swi{subs,5}==0);
        %sync model
        Pr_Sw_inp_3_sync{1}(subs,1) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==1));
        Pr_Sw_inp_3_sync{1}(subs,2) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==2));
        Pr_Sw_inp_3_sync{1}(subs,3) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==3));
        Pr_Sw_inp_3_sync{1}(subs,4) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==4));
        Pr_Sw_inp_3_sync{1}(subs,5) = mean(pr_of_swi{subs,4}(pr_of_swi{subs,1}==5));
        %subs data
        Pr_Sw_inp_3_sync{2}(subs,1) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==1)/sum(pr_of_swi{subs,1}==1);
        Pr_Sw_inp_3_sync{2}(subs,2) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==2)/sum(pr_of_swi{subs,1}==2);
        Pr_Sw_inp_3_sync{2}(subs,3) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==3)/sum(pr_of_swi{subs,1}==3);
        Pr_Sw_inp_3_sync{2}(subs,4) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==4)/sum(pr_of_swi{subs,1}==4);
        Pr_Sw_inp_3_sync{2}(subs,5) = sum((Input(subs).TF( Input(subs).StateExp==1))' & pr_of_swi{subs,1}==5)/sum(pr_of_swi{subs,1}==5);
        %Bayesian simulation data
        Pr_Sw_inp_3_sync{3}(subs,1) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.1))/sum(MachineSimulation(subs).tDev==0.1 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.1));
        Pr_Sw_inp_3_sync{3}(subs,2) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.3))/sum(MachineSimulation(subs).tDev==0.3 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.3));
        Pr_Sw_inp_3_sync{3}(subs,3) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.5))/sum(MachineSimulation(subs).tDev==0.5 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.5));
        Pr_Sw_inp_3_sync{3}(subs,4) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.7))/sum(MachineSimulation(subs).tDev==0.7 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.7));
        Pr_Sw_inp_3_sync{3}(subs,5) = sum(MachineSimulation(subs).TF( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.9))/sum(MachineSimulation(subs).tDev==0.9 & MachineSimulation(subs).StateExp==1) * mean(MachineSimulation(subs).expectedAccuracy( MachineSimulation(subs).StateExp==1 & MachineSimulation(subs).tDev==0.9));

end

set(0, 'DefaultAxesFontName', 'Helvetica');   % Fonttype for axis
set(0, 'DefaultTextFontName', 'Helvetica');   % Fonttype for text
addpath(genpath('/Users/vincentxu/Desktop/MRI_behavioural/Hiearchica_Reasoning_model/depend/Functions'));
load('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/model_res_switch_fmri_m_4_29sub.mat')
%
rn=[];
figure
set(gcf,'unit','pixels','position',[200 50 450 400]);
subplot(1,1,1)
    set(gca,'FontSize',14,'tickdir','out')

fillsteplotgreen(Pr_Sw_inp_3_sync{1},4); hold on
fillsteplotblue(Pr_Sw_inp_3_sync{2},4); hold on
fillsteplotred(Pr_Sw_inp_3_sync{3},4); hold on

title('Model Comparison','fontsize',18,'FontWeight','Normal')

% print('Pr_switch_bayesian_rnn0', '-dsvg','-painters')
