clear;close all
c    = [[255/255 51/255 51/255];[051/255 153/255 255/255]]; % set colors
set(0, 'DefaultAxesFontName', 'Helvetica');   % Fonttype for axis
set(0, 'DefaultTextFontName', 'Helvetica');   % Fonttype for text

load('psychometric_fmri_exp_behvaior_29subj.mat');%******
load('model_res_psychometric_fmri_60sub_pAlphErr1_test3.mat')
thrs = thrs; %Difficulty levels
isSave=1; %save figures

addpath(genpath('Hiearchica_Reasoning_model/depend/Hierarchica_Reasoning_model_changeNB'));
addpath(genpath('Hiearchica_Reasoning_model/depend/FMINSEARCHBND'));
rng(66,'twister');

%% load data
    addpath(genpath('Hiearchica_Reasoning_model/depend/Functions'));
    load('model_res_switch_fmri_m_4_29sub.mat');
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



%% logistic regression (pr of only 1-er, 2-er); Interaction, model comparison
% logistic regression (binary)
%load
load('model_res_switch_fmri_m_4_29sub.mat');
beha_data_ori = readtable('data_29sub_raw.xlsx');
%
beta_errors = []; beta_diffic = []; infer_data = [];
beha_data=Input;
for subs = 1:length(beha_data)
    infer_data(subs).ID = repelem(subs, sum(beha_data(subs).StateExp==1))';

    % Keep original 5-level anti-probability for plotting / binning
    raw_tDev = beha_data(subs).tDev(beha_data(subs).StateExp==1);
    infer_data(subs).Dev5 = raw_tDev;

    % Convert to symmetric difficulty for logistic predictor
    infer_data(subs).Diffs = raw_tDev;
    infer_data(subs).Diffs(raw_tDev==0.1)=3;
    infer_data(subs).Diffs(raw_tDev==0.3)=2;
    infer_data(subs).Diffs(raw_tDev==0.5)=1;
    infer_data(subs).Diffs(raw_tDev==0.7)=2;
    infer_data(subs).Diffs(raw_tDev==0.9)=3;

    % infer_data(subs).SW = beha_data(subs).SW(beha_data(subs).StateExp==1);
    tri_idx = beha_data(subs).trialNum(beha_data(subs).StateExp==1);
    ruleSW = beha_data_ori.ruleSwitch(beha_data_ori.idIdx ==beha_data(subs).ID  & beha_data_ori.expIdx==1);    
    infer_data(subs).SW    = ruleSW(tri_idx);
    infer_data(subs).Nback = beha_data(subs).Nback_act(beha_data(subs).StateExp==1);
end

pr_cat = [vertcat(infer_data.ID), vertcat(infer_data.Nback), ...
          vertcat(infer_data.Diffs), vertcat(infer_data.Dev5), ...
          vertcat(infer_data.SW)];

pr_cat(pr_cat(:,2)==0,:) = nan;
pr_cat(pr_cat(:,2)>2,:)  = nan;
pr_cat = rmmissing(pr_cat);
% pr_cat(:,2:3) = zscore(pr_cat(:,2:3));

Pr_Subs = array2table(pr_cat, ...
    'VariableNames', {'subID', 'Ners', 'Diffs', 'Dev5', 'Sw'}); % make table
% Pr_Subs.Diffs = categorical(Pr_Subs.Diffs);
% Pr_Subs.Ners = categorical(Pr_Subs.Ners);
% Pr_Subs.subID = categorical(Pr_Subs.subID);

% Keep this full-sample GLME only for descriptive group-level effect display
% Cross-validation below is subject-wise, to match the CBM CV in m_7_4.
glme = fitglme(Pr_Subs, 'Sw ~ 1 + Ners + Diffs + (1 + Ners + Diffs|subID)', ...
  'Distribution','binomial','Link','logit','FitMethod','Laplace');

glme_int = fitglme(Pr_Subs, 'Sw ~ 1 + Ners + Diffs + Ners * Diffs + (1 + Ners + Diffs + Ners * Diffs |subID)', ...
  'Distribution','binomial','Link','logit','FitMethod','Laplace');


disp(glme);
disp(glme_int);
glme.ModelCriterion
glme_int.ModelCriterion
% [beta, betaCI, stats] = fixedEffects(glme); exp(beta);
[B, Bnames] = randomEffects(glme)