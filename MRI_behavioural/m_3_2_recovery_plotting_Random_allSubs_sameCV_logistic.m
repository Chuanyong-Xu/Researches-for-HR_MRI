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


%% Full model: out-of-sample testing, parameter recovery
% 2-fold CV, repeat 10 times

load data_29sub_raw.mat
% % % %[data,titleName]=xlsread('.\behavior_pilot_fmri_raw.xlsx','inferred1_instructed0');%******
idIdx         = 1;
expIdxIdx     = 3;
trialnumIdx   = 4;
rulelistIdx   = 5;
refAngIdx     = 6;
testAngIdx    = 7;
bias_levelIdx = 8;
ruleRespIdx   = 9;
pcptRespIdx   = 11;
ruleSwitchIdx = 13;
feedbackIdx   = 14

% % % % % % thrs          = [.01, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, .99];%********************************
thrs          = [0.1 0.3 0.5 0.7 0.9];%[0.01 1/4 1/2 3/4 0.99]; %[1/6 2/6 3/6 4/6 5/6]; %[0.0675, 5/16, 1/2, 11/16, .9325];%********************************
delMiss       = data;
delMiss(unique([find(isnan(data(:,ruleRespIdx)));find(isnan(data(:,pcptRespIdx)))]),:) = [];
id            = unique(delMiss(:,idIdx));

%% prepare data set
for iSubject = 1:length(id)
    InputAll(iSubject).ID = id(iSubject); % name of subject
    Index_of_subject = find(delMiss(:,idIdx)==id(iSubject));
    InputAll(iSubject).StateExp = delMiss(Index_of_subject,expIdxIdx); % defines the type of experiment (0: instructed task; 1: inferred task)
    InputAll(iSubject).Cue  = 1-delMiss(Index_of_subject, rulelistIdx); % ****** actual rule of environment at each trial(0:red; 1:blue)
    InputAll(iSubject).PrAn = 1-delMiss(Index_of_subject, pcptRespIdx); % ****** action (0: bottom, 1: up)
    InputAll(iSubject).tDev = thrs(delMiss(Index_of_subject, bias_levelIdx))'; % sample intervals
    InputAll(iSubject).RuleChoice = 1-delMiss(Index_of_subject, ruleRespIdx); % rule choice by subject
    InputAll(iSubject).TF = delMiss(Index_of_subject, feedbackIdx); % Feedback (accuracy of both action and rule)
    InputAll(iSubject).SW = delMiss(Index_of_subject, ruleSwitchIdx);
    InputAll(iSubject).DevValues = thrs; % sample interval ranges
    InputAll(iSubject).trialNum = delMiss(Index_of_subject, trialnumIdx); % for test

end


psycho = load('model_res_psychometric_fmri_60sub_pAlphErr1_test3.mat');
model  = psycho.model;

%% fit the belief update model (this is after we fitted the psychometric functions):
for iSubject = 1:length(id)
    % step #6: Sort the trials based on the coditions (1-B error,  2B error) in the inferred task
    temp_parameters_as_input(:) = model.psychometric.subj(iSubject).EstimatedParameters(2, :); % Subjective psychometric parameters. Index can be 1 (instructed task) or 2 (inferred task)
    temp_parameter_strings = model.psychometric.subj(iSubject).parameter_strings; % name of fitted parameters
    temp_SubjObj_flag = model.psychometric.subj(iSubject).SubjObj_flag; % The tag shows subjective or objective fit
    
    % prepare data format to feed the modelling sub-function
    % extract 1-B error trials (1-B Error and 2-B Reward)
    conditions_1Back_index = find( (InputAll(iSubject).StateExp(1:end-2) == 1) .* (InputAll(iSubject).TF(1:end-2)==1) .* (InputAll(iSubject).TF(2:end-1)==0) ) +1; % RW, ER  => {1B-Er}
    %conditions_2Back_index = find( (InputAll(iSubject).StateExp(1:end-3) == 1) .* (InputAll(iSubject).TF(1:end-3)==1) .* (InputAll(iSubject).TF(2:end-2)==0) .* (InputAll(iSubject).TF(3:end-1)==0) ) +2; % RW, ER, ER  => {2B-Er}
    conditions_2Back_index = find( (InputAll(iSubject).StateExp(1:end-3) == 1) .* (InputAll(iSubject).TF(1:end-3)==1) .* (InputAll(iSubject).TF(2:end-2)==0) .* (InputAll(iSubject).TF(3:end-1)==0) .* (InputAll(iSubject).RuleChoice(2:end-2)==InputAll(iSubject).RuleChoice(3:end-1)) ) +2; % RW, ER, ER  => {2B-Er} changed by wenshan 20220518
%     conditions_2Back_index = find( (InputAll(iSubject).StateExp(1:end-3) == 1) .* (InputAll(iSubject).TF(1:end-3)==1) .* (InputAll(iSubject).TF(2:end-2)==0) .* (InputAll(iSubject).TF(3:end-1)==0) ) +2; % RW, ER, ER  => {2B-Er} changed by wenshan 20220518


    clear switchInput
    
        ii_dev_2=0;
        ii_dev_1=0;
    % making a dataset for the modelling sub-function:
    for iTrial =1: length(conditions_1Back_index)
        switchInput(iTrial).T = 1; % number of 1Back errors
        switchInput(iTrial).tDev = InputAll(iSubject).tDev(conditions_1Back_index(iTrial)); % stimulus
        switchInput(iTrial).PrAn = InputAll(iSubject).PrAn(conditions_1Back_index(iTrial)); % Action
        switchInput(iTrial).RuleChoice = InputAll(iSubject).RuleChoice(conditions_1Back_index(iTrial)); % rule choice
        switchInput(iTrial).Cue = InputAll(iSubject).Cue(conditions_1Back_index(iTrial)); % actual rule
        switchInput(iTrial).TF = InputAll(iSubject).TF(conditions_1Back_index(iTrial)); % feedback
        switchInput(iTrial).SW = InputAll(iSubject).RuleChoice(conditions_1Back_index(iTrial)) ~= InputAll(iSubject).RuleChoice(conditions_1Back_index(iTrial)+1); % switch/notswitch
        switchInput(iTrial).DevValues = InputAll(iSubject).DevValues; % array of stimulus (samples)
    
        ii_dev_1=ii_dev_1+1;%
        input_dev_1back(ii_dev_1,1) = InputAll(iSubject).tDev(conditions_1Back_index(iTrial));


    end
    
    for iTrial = 1: length(conditions_2Back_index)
        switchInput(length(conditions_1Back_index)+iTrial).T = 2; % number of 2Back errors
        switchInput(length(conditions_1Back_index)+iTrial).tDev = [InputAll(iSubject).tDev(conditions_2Back_index(iTrial) -1), InputAll(iSubject).tDev(conditions_2Back_index(iTrial))];
        switchInput(length(conditions_1Back_index)+iTrial).PrAn = [InputAll(iSubject).PrAn(conditions_2Back_index(iTrial) -1), InputAll(iSubject).PrAn(conditions_2Back_index(iTrial))];
        switchInput(length(conditions_1Back_index)+iTrial).RuleChoice = [InputAll(iSubject).RuleChoice(conditions_2Back_index(iTrial) -1), InputAll(iSubject).RuleChoice(conditions_2Back_index(iTrial))];
        switchInput(length(conditions_1Back_index)+iTrial).Cue = [InputAll(iSubject).Cue(conditions_2Back_index(iTrial) -1), InputAll(iSubject).Cue(conditions_2Back_index(iTrial))];
        switchInput(length(conditions_1Back_index)+iTrial).TF = [InputAll(iSubject).TF(conditions_2Back_index(iTrial) -1), InputAll(iSubject).TF(conditions_2Back_index(iTrial))];
        switchInput(length(conditions_1Back_index)+iTrial).SW = InputAll(iSubject).RuleChoice(conditions_2Back_index(iTrial)) ~= InputAll(iSubject).RuleChoice(conditions_2Back_index(iTrial)+1);
        switchInput(length(conditions_1Back_index)+iTrial).DevValues = InputAll(iSubject).DevValues;
    
        ii_dev_2=ii_dev_2+1;%
        input_dev_2back(ii_dev_2,1:2) = [InputAll(iSubject).tDev(conditions_2Back_index(iTrial) -1), InputAll(iSubject).tDev(conditions_2Back_index(iTrial))];
    end

    input{iSubject} = switchInput;  

    % Fit the switch model
    % ****** by xu for model comparison, change the lamada (alpha_transition) in SwitchParameterOptimization.m and pr_switch_func.m;

    Set_param_variable = [1, 0.5, 0.5];%[1, 0.5, 0.5]; [4,.5,.5];%sigma,scale,pam3
     [output(iSubject), estimate_parameters(iSubject), MachineSimulation(iSubject),Input(iSubject),MLE(iSubject,1)] = SwitchParameterOptimization(switchInput, temp_parameters_as_input, temp_parameter_strings, temp_SubjObj_flag, InputAll(iSubject), Set_param_variable);
    % [~, ~, estimate_parameters, ~, ~, MachineSimulation] = SwitchParameterOptimization(switchInput, temp_parameters_as_input, temp_parameter_strings, temp_SubjObj_flag, InputAll(iSubject), Set_param_variable);       % pr_switch_model(iRuleChoice, itd, T_numOfBackError )
    % model parameters:
    % 1- estimate_parameters.sigma_switch_estimated (sigma parameter of the switch model) (Free parameter for optimization)
    % 2- estimate_parameters.pam3_estimated (alpha of the switch model) (Free parameter for optimization)
    % hazard parameter is fixed.
    % consecutive errors (for example. for T=2 => sigma2= sigma * sqrt(2), ...)
    

%% =========================================================
% Fig. S5B-like 2-fold cross-validation for switch model
% + subject-wise logistic benchmark using exactly the same train/test split
% Fit on random half of switchInput; test on the remaining half
%% =========================================================

nCV = 100;   % repeated 2-fold CV
nT  = 2;     % 1B-Er and 2B-Er
nDev = length(thrs);

if iSubject == 1
    CV_Pr_Sw_obs   = nan(length(id), nCV, nT, nDev);
    CV_Pr_Sw_model = nan(length(id), nCV, nT, nDev);
    CV_test_NLL    = nan(length(id), nCV);

    % Logistic benchmark, using the same switchInput and the same CV split
    LR_Pr_Sw_obs   = nan(length(id), nCV, nT, nDev);
    LR_Pr_Sw_model = nan(length(id), nCV, nT, nDev);
    LR_test_NLL    = nan(length(id), nCV);
end

% ---- expected accuracy benchmark, same as inside SwitchParameterOptimization.m ----
expectedAccuracy_Benchmark = [];
flag = 0;
for iCue = 0:1
    for tDev_i = 1:length(InputAll(iSubject).DevValues)
        flag = flag + 1;
        expectedAccuracy_BenchmarkInput.RuleChoice = iCue;
        expectedAccuracy_BenchmarkInput.tDev = InputAll(iSubject).DevValues(tDev_i);

        expectedAccuracy_Benchmark.RuleChoice(flag,1) = iCue;
        expectedAccuracy_Benchmark.tDev(flag,1) = InputAll(iSubject).DevValues(tDev_i);
        expectedAccuracy_Benchmark.expectedAccuracy(flag,1) = expected_Accuracy( ...
            expectedAccuracy_BenchmarkInput, ...
            temp_parameters_as_input, ...
            temp_parameter_strings, ...
            temp_SubjObj_flag);
    end
end

% ---- labels for stratified random half split ----
T_all = [switchInput.T]';   % 1 or 2
nTrial_switch = length(switchInput);

for iCV = 1:nCV

    trainMask = false(nTrial_switch,1);

    % Stratified 2-fold split: half of 1B-Er and half of 2B-Er go to training.
    % This single trainMask is shared by CBM and logistic.
    for iT = 1:2
        idxT = find(T_all == iT);
        idxT = idxT(randperm(length(idxT)));

        nTrainT = floor(length(idxT)/2);
        trainMask(idxT(1:nTrainT)) = true;
    end

    testMask = ~trainMask;

    switchInput_train = switchInput(trainMask);
    switchInput_test  = switchInput(testMask);

    % =====================================================
    % 1) CBM: fit model on training half
    % =====================================================
    [~, estimate_parameters_cv, ~, ~, ~] = SwitchParameterOptimization( ...
        switchInput_train, ...
        temp_parameters_as_input, ...
        temp_parameter_strings, ...
        temp_SubjObj_flag, ...
        InputAll(iSubject), ...
        Set_param_variable);

    % ---- convert fitted struct back to parameter vector for pr_switch_func ----
    switchPam_cv = [ ...
        estimate_parameters_cv.sigma_switch_estimated, ...
        estimate_parameters_cv.alpha_transition, ...
        estimate_parameters_cv.pam3_estimated];

    parameter_strings_cv = { ...
        'sigma_switch = switchPam(1);', ...
        'alpha_transition = switchPam(2);', ...
        'pam3 = switchPam(3);'};

    % ---- predict Pr(Switch) on held-out test half ----
    [Output_pr_of_switch_cv, Output_tDev_lastOne_cv, ~, Output_T_cv, Output_SW_cv, ~] = pr_switch_func( ...
        switchInput_test, ...
        switchPam_cv, ...
        parameter_strings_cv, ...
        temp_parameters_as_input, ...
        temp_parameter_strings, ...
        temp_SubjObj_flag, ...
        expectedAccuracy_Benchmark);

    % ---- CBM test-set negative log-likelihood ----
    p = Output_pr_of_switch_cv(:);
    y = Output_SW_cv(:);

    p = min(max(p, 1e-4), 1-1e-4);  % avoid log(0)
    CV_test_NLL(iSubject,iCV) = -mean(y .* log(p) + (1-y) .* log(1-p));

    % ---- group CBM held-out predictions by T and td ----
    for iT = 1:2
        for itd = 1:nDev

            idx = (Output_T_cv(:) == iT) & ...
                  (abs(Output_tDev_lastOne_cv(:) - thrs(itd)) < 1e-10);

            if any(idx)
                CV_Pr_Sw_obs(iSubject,iCV,iT,itd)   = mean(Output_SW_cv(idx));
                CV_Pr_Sw_model(iSubject,iCV,iT,itd) = mean(Output_pr_of_switch_cv(idx));
            end

        end
    end

    % =====================================================
    % 2) Logistic benchmark: same training/test trials
    %    Predictor: Sw ~ 1 + Ners + Diffs, no interaction
    % =====================================================

    % ---- convert switchInput_train to table ----
    nTrain = length(switchInput_train);
    Ners_train  = nan(nTrain,1);
    Dev5_train  = nan(nTrain,1);
    Diffs_train = nan(nTrain,1);
    Sw_train    = nan(nTrain,1);

    for k = 1:nTrain
        Ners_train(k) = switchInput_train(k).T;
        tmpDev = switchInput_train(k).tDev;
        Dev5_train(k) = tmpDev(end);       % for 2B-Er, use the second error's tDev
        Sw_train(k) = switchInput_train(k).SW;

        if abs(Dev5_train(k) - 0.1) < 1e-10
            Diffs_train(k) = 3;
        elseif abs(Dev5_train(k) - 0.3) < 1e-10
            Diffs_train(k) = 2;
        elseif abs(Dev5_train(k) - 0.5) < 1e-10
            Diffs_train(k) = 1;
        elseif abs(Dev5_train(k) - 0.7) < 1e-10
            Diffs_train(k) = 2;
        elseif abs(Dev5_train(k) - 0.9) < 1e-10
            Diffs_train(k) = 3;
        end
    end

    T_train_lr = table(Ners_train, Diffs_train, Dev5_train, Sw_train, ...
        'VariableNames', {'Ners','Diffs','Dev5','Sw'});

    % ---- convert switchInput_test to table ----
    nTest = length(switchInput_test);
    Ners_test  = nan(nTest,1);
    Dev5_test  = nan(nTest,1);
    Diffs_test = nan(nTest,1);
    Sw_test    = nan(nTest,1);

    for k = 1:nTest
        Ners_test(k) = switchInput_test(k).T;
        tmpDev = switchInput_test(k).tDev;
        Dev5_test(k) = tmpDev(end);        % for 2B-Er, use the second error's tDev
        Sw_test(k) = switchInput_test(k).SW;

        if abs(Dev5_test(k) - 0.1) < 1e-10
            Diffs_test(k) = 3;
        elseif abs(Dev5_test(k) - 0.3) < 1e-10
            Diffs_test(k) = 2;
        elseif abs(Dev5_test(k) - 0.5) < 1e-10
            Diffs_test(k) = 1;
        elseif abs(Dev5_test(k) - 0.7) < 1e-10
            Diffs_test(k) = 2;
        elseif abs(Dev5_test(k) - 0.9) < 1e-10
            Diffs_test(k) = 3;
        end
    end

    T_test_lr = table(Ners_test, Diffs_test, Dev5_test, Sw_test, ...
        'VariableNames', {'Ners','Diffs','Dev5','Sw'});

    % Save the held-out observed switch rate before fitting logistic.
    % These observed points are based on exactly the same test trials as CBM.
    for iT = 1:2
        for itd = 1:nDev
            idx = T_test_lr.Ners == iT & ...
                  abs(T_test_lr.Dev5 - thrs(itd)) < 1e-10;

            if any(idx)
                LR_Pr_Sw_obs(iSubject,iCV,iT,itd) = mean(T_test_lr.Sw(idx));
            end
        end
    end

    % Logistic needs both response classes in training.
    if ~isempty(T_train_lr) && ~isempty(T_test_lr) && numel(unique(T_train_lr.Sw)) >= 2
        try
            glm_lr_cv = fitglm(T_train_lr, ...
                'Sw ~ 1 + Ners + Diffs', ...
                'Distribution','binomial', ...
                'Link','logit', ...
                'Options', statset('MaxIter', 5000));

            p_pred_lr = predict(glm_lr_cv, T_test_lr);
            p_pred_lr = min(max(p_pred_lr, 1e-4), 1-1e-4);

            y_obs_lr = T_test_lr.Sw;
            LR_test_NLL(iSubject,iCV) = -mean(y_obs_lr .* log(p_pred_lr) + ...
                                              (1-y_obs_lr) .* log(1-p_pred_lr));

            for iT = 1:2
                for itd = 1:nDev

                    idx = T_test_lr.Ners == iT & ...
                          abs(T_test_lr.Dev5 - thrs(itd)) < 1e-10;

                    if any(idx)
                        LR_Pr_Sw_model(iSubject,iCV,iT,itd) = mean(p_pred_lr(idx));
                    end

                end
            end

        catch ME
            fprintf('Logistic CV failed: subject %d, CV %d, %s\n', iSubject, iCV, ME.message);
        end
    end

end



    % measuring switch probability of subject and also the machineSimulation
    % for subject:
    input_m(iSubject).T = [];
    input_m(iSubject).tDev = [];
    input_m(iSubject).pr_of_switch = [];
    for itd = 1: length(InputAll(iSubject).DevValues) % for each sample interval:
        % for rw:
        %# Finding the index of RuleChoice has switched; 
        %# SW_inp is from the model fitting (Input from SwitchParameterOptimization.m);
        %# SW_v0 is from the original dataset (InputAll);
        %# Pr_Sw calculate the original switch rate under the trial whose preceding trial was rewarded (InputAll)
        %# Inp_Sw calculate the switch rate under the trial whose preceding trial was rewarded from model fitting data (Input)
        %# Inp_Pr_Sw is the switch rate from model fitting data ()
        conditions_0Back_index = find( (InputAll(iSubject).StateExp(1:end-1) == 1) .* (InputAll(iSubject).TF(1:end-1)==1).* (InputAll(iSubject).tDev(1:end-1)==InputAll(iSubject).DevValues(itd)) ) ; % RW, ER  => {1B-Er}
        SW_v0 = InputAll(iSubject).RuleChoice(conditions_0Back_index) ~= InputAll(iSubject).RuleChoice(conditions_0Back_index+1);
        SW_inp = Input(iSubject).RuleChoice(conditions_0Back_index)  ~= Input(iSubject).RuleChoice(conditions_0Back_index+1);
        iT = 3;
        Pr_Sw{iT}(iSubject, itd) = sum(SW_v0) / length(SW_v0); res_o{iT}.res_switch(iSubject, itd) = sum(SW_v0);  res_o{iT}.res_count(iSubject, itd) = length(SW_v0);
        Inp_Sw{iT}(iSubject, itd) = sum(SW_inp) / length(SW_v0); res_inp{iT}.res_switch(iSubject, itd) = sum(SW_inp);  res_inp{iT}.res_count(iSubject, itd) = length(SW_inp);
% % % % % %         Inp_Pr_Sw{iT}(iSubject, itd) = mean(Input(iSubject).pr_of_switch(conditions_1Back_index)); % by xu 20231012？？？？？？？？
        Inp_Pr_Sw{iT}(iSubject, itd) = mean(Input(iSubject).pr_of_switch(conditions_0Back_index));

        % for 1-Back:
        conditions_1Back_index = find( (InputAll(iSubject).StateExp(1:end-2) == 1) .* (InputAll(iSubject).TF(1:end-2)==1) .* (InputAll(iSubject).TF(2:end-1)==0) .* (InputAll(iSubject).tDev(2:end-1)==InputAll(iSubject).DevValues(itd)) ) +1; % RW, ER  => {1B-Er}
        SW_v1 = InputAll(iSubject).RuleChoice(conditions_1Back_index) ~= InputAll(iSubject).RuleChoice(conditions_1Back_index+1);
        SW_inp = Input(iSubject).RuleChoice(conditions_1Back_index)  ~= Input(iSubject).RuleChoice(conditions_1Back_index+1);
        iT = 1; 
        Pr_Sw{iT}(iSubject, itd) = sum(SW_v1) / length(SW_v1); res_o{iT}.res_switch(iSubject, itd) = sum(SW_v1);  res_o{iT}.res_count(iSubject, itd) = length(SW_v1);
        Inp_Sw{iT}(iSubject, itd) = sum(SW_inp) / length(SW_v1); res_inp{iT}.res_switch(iSubject, itd) = sum(SW_inp);  res_inp{iT}.res_count(iSubject, itd) = length(SW_inp);
        Inp_Pr_Sw{iT}(iSubject, itd) = mean(Input(iSubject).pr_of_switch(conditions_1Back_index));
        
        conditions_1Back_index2{iSubject, itd}=conditions_1Back_index; %****************        
        % for 2-Back Er (cosecutive errors with no 1-B switch):
        conditions_2Back_index = find( (InputAll(iSubject).StateExp(1:end-3) == 1) .* (InputAll(iSubject).TF(1:end-3)==1) .* (InputAll(iSubject).TF(2:end-2)==0) .* (InputAll(iSubject).TF(3:end-1)==0) .* (InputAll(iSubject).RuleChoice(2:end-2) == InputAll(iSubject).RuleChoice(3:end-1)) .* (InputAll(iSubject).tDev(3:end-1)==InputAll(iSubject).DevValues(itd)) ) +2; % RW, ER, ER  => {2B-Er}
%        conditions_2Back_index = find( (InputAll(iSubject).StateExp(1:end-3) == 1) .* (InputAll(iSubject).TF(1:end-3)==1) .* (InputAll(iSubject).TF(2:end-2)==0) .* (InputAll(iSubject).TF(3:end-1)==0).* (InputAll(iSubject).tDev(3:end-1)==InputAll(iSubject).DevValues(itd)) ) +2; % RW, ER, ER  => {2B-Er} changed by wenshan 20220518

        SW_v2 = InputAll(iSubject).RuleChoice(conditions_2Back_index) ~= InputAll(iSubject).RuleChoice(conditions_2Back_index+1);
        SW_inp = Input(iSubject).RuleChoice(conditions_2Back_index)  ~= Input(iSubject).RuleChoice(conditions_2Back_index+1);
        iT = 2; 
        Pr_Sw{iT}(iSubject, itd) = sum(SW_v2) / length(SW_v2); res_o{iT}.res_switch(iSubject, itd) = sum(SW_v2);  res_o{iT}.res_count(iSubject, itd) = length(SW_v2);
        Inp_Sw{iT}(iSubject, itd) = sum(SW_inp) / length(SW_v2); res_inp{iT}.res_switch(iSubject, itd) = sum(SW_inp);  res_inp{iT}.res_count(iSubject, itd) = length(SW_inp);% length(SW_v2);%changed by wenshan 20220810
        Inp_Pr_Sw{iT}(iSubject, itd) = mean(Input(iSubject).pr_of_switch(conditions_2Back_index));
        
        conditions_2Back_index2{iSubject, itd}=conditions_2Back_index; %****************
        
        accuracy_infer{iSubject, itd} = sum((InputAll(iSubject).StateExp(1:end) == 1) .* (InputAll(iSubject).TF(1:end)==1).* (InputAll(iSubject).tDev(1:end)==InputAll(iSubject).DevValues(itd)) ) / sum((InputAll(iSubject).StateExp(1:end) == 1) .* (InputAll(iSubject).tDev(1:end)==InputAll(iSubject).DevValues(itd)) ); % RW, ER  => {1B-Er}
        accuracy_instr{iSubject, itd} = sum((InputAll(iSubject).StateExp(1:end) == 0) .* (InputAll(iSubject).TF(1:end)==1).* (InputAll(iSubject).tDev(1:end)==InputAll(iSubject).DevValues(itd)) ) / sum((InputAll(iSubject).StateExp(1:end) == 0) .* (InputAll(iSubject).tDev(1:end)==InputAll(iSubject).DevValues(itd)) ); % RW, ER  => {1B-Er}

        confidence_mu_2Back{iSubject, itd}=mean(Input(iSubject).mu_switch_estimated(conditions_2Back_index));
        confidence_mu_1Back{iSubject, itd}=mean(Input(iSubject).mu_switch_estimated(conditions_1Back_index));
        
        Pr_An_2Back{iSubject, itd}=mean(Input(iSubject).PrAn(conditions_2Back_index));
        Pr_An_1Back{iSubject, itd}=mean(Input(iSubject).PrAn(conditions_1Back_index));


    end
    
    % for machine simulation:
    for itd = 1: length(InputAll(iSubject).DevValues) % for each sample interval:
        % for rw:
        conditions_0Back_index = find( (MachineSimulation(iSubject).StateExp(1:end-1) == 1) .* (MachineSimulation(iSubject).TF(1:end-1)==1) .* (MachineSimulation(iSubject).tDev(1:end-1)==InputAll(iSubject).DevValues(itd)) ) ; % RW, ER  => {1B-Er}
        SW_v0 = MachineSimulation(iSubject).RuleChoice(conditions_0Back_index) ~= MachineSimulation(iSubject).RuleChoice(conditions_0Back_index+1);
        iT = 3; 
        M_Sw{iT}(iSubject, itd) = sum(SW_v0) / length(SW_v0); res_m{iT}.res_switch(iSubject, itd) = sum(SW_v0);  res_m{iT}.res_count(iSubject, itd) = length(SW_v0);
        M_Pr_Sw{iT}(iSubject, itd) = mean(MachineSimulation(iSubject).pr_of_switch(conditions_0Back_index));
        
        % for 1-Back:
        conditions_1Back_index = find( (MachineSimulation(iSubject).StateExp(1:end-2) == 1) .* (MachineSimulation(iSubject).TF(1:end-2)==1) .* (MachineSimulation(iSubject).TF(2:end-1)==0) .* (MachineSimulation(iSubject).tDev(2:end-1)==InputAll(iSubject).DevValues(itd)) ) +1; % RW, ER  => {1B-Er}
        SW_v1 = MachineSimulation(iSubject).RuleChoice(conditions_1Back_index) ~= MachineSimulation(iSubject).RuleChoice(conditions_1Back_index+1);
        iT = 1; 
        M_Sw{iT}(iSubject, itd) = sum(SW_v1) / length(SW_v1); res_m{iT}.res_switch(iSubject, itd) = sum(SW_v1);  res_m{iT}.res_count(iSubject, itd) = length(SW_v1);
        M_Pr_Sw{iT}(iSubject, itd) = mean(MachineSimulation(iSubject).pr_of_switch(conditions_1Back_index));
        
        % for 2-Back Er (cosecutive errors with no 1-B switch):
        conditions_2Back_index = find( (MachineSimulation(iSubject).StateExp(1:end-3) == 1) .* (MachineSimulation(iSubject).TF(1:end-3)==1) .* (MachineSimulation(iSubject).TF(2:end-2)==0) .* (MachineSimulation(iSubject).TF(3:end-1)==0) .* (MachineSimulation(iSubject).RuleChoice(2:end-2) == MachineSimulation(iSubject).RuleChoice(3:end-1)) .* (MachineSimulation(iSubject).tDev(3:end-1)==InputAll(iSubject).DevValues(itd)) ) +2; % RW, ER, ER  => {2B-Er}
%         conditions_2Back_index = find( (MachineSimulation(iSubject).StateExp(1:end-3) == 1) .* (MachineSimulation(iSubject).TF(1:end-3)==1) .* (MachineSimulation(iSubject).TF(2:end-2)==0) .* (MachineSimulation(iSubject).TF(3:end-1)==0) .* (MachineSimulation(iSubject).tDev(3:end-1)==InputAll(iSubject).DevValues(itd)) ) +2; % RW, ER, ER  => {2B-Er}
        
        SW_v2 = MachineSimulation(iSubject).RuleChoice(conditions_2Back_index) ~= MachineSimulation(iSubject).RuleChoice(conditions_2Back_index+1);
        iT = 2; 
        M_Sw{iT}(iSubject, itd) = sum(SW_v2) / length(SW_v2); res_m{iT}.res_switch(iSubject, itd) = sum(SW_v2);  res_m{iT}.res_count(iSubject, itd) = length(SW_v2);
        M_Pr_Sw{iT}(iSubject, itd) = mean(MachineSimulation(iSubject).pr_of_switch(conditions_2Back_index));
    end
    
end

save('model_res_switch_CV.mat', ...
     'CV_Pr_Sw_obs', 'CV_Pr_Sw_model', 'CV_test_NLL', ...
     'LR_Pr_Sw_obs', 'LR_Pr_Sw_model', 'LR_test_NLL', ...
     'thrs', 'id');


%% =========================================================
% Fig. S5B-like: group-level 2-fold cross-validated prediction
% CBM and logistic use the same held-out trials in each CV split
%% =========================================================

load('model_res_switch_CV.mat');
% CV_Pr_Sw_obs   : subj x CV x T x difficulty
% CV_Pr_Sw_model : subj x CV x T x difficulty
% LR_Pr_Sw_obs   : subj x CV x T x difficulty
% LR_Pr_Sw_model : subj x CV x T x difficulty

% This figure is kept for the original layout:
% parameter recovery occupies subplots 1-3 later; CBM CV occupies subplot 4.
fig_cvrec = figure('Color','w','Position',[100 100 1000 250]);
subplot(1,4,4); hold on;

% average across CV repeats first, then across subjects
CV_obs_subj   = squeeze(nanmean(CV_Pr_Sw_obs, 2));     % subj x T x difficulty
CV_model_subj = squeeze(nanmean(CV_Pr_Sw_model, 2));   % subj x T x difficulty

% colors: use the same color rule as your previous plot
c_cv = {[0.2784 0.4471 0.451], [0.7412 0.3569 0.0235]};
% c_cv{1}: 1B-Er
% c_cv{2}: 2B-Er

xvalue = 1:5;

for iT = 1:2

    y_obs = squeeze(CV_obs_subj(:,iT,:));      % subj x difficulty
    y_mod = squeeze(CV_model_subj(:,iT,:));    % subj x difficulty

    obs_mean = nanmean(y_obs, 1);
    obs_sem  = nanstd(y_obs, [], 1) ./ sqrt(sum(~isnan(y_obs), 1));

    mod_mean = nanmean(y_mod, 1);
    mod_sem  = nanstd(y_mod, [], 1) ./ sqrt(sum(~isnan(y_mod), 1));

    % observed held-out data: points + error bars
    errorbar(xvalue, obs_mean, obs_sem, 'ko', ...
        'LineWidth', 2, ...
        'MarkerFaceColor', c_cv{iT}, ...
        'MarkerSize', 10);
    hold on;

    % model held-out prediction: line + shaded SEM
    xx = [xvalue, fliplr(xvalue)];
    yy = [mod_mean + mod_sem, fliplr(mod_mean - mod_sem)];

    h_fill(iT) = fill(xx, yy, c_cv{iT}, ...
        'FaceAlpha', 0.25, ...
        'EdgeColor', 'none');
    hold on;

    h_line(iT) = plot(xvalue, mod_mean, '-', ...
        'Color', c_cv{iT}, ...
        'LineWidth', 2.5);
    hold on;

end

set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.1','0.3','0.5','0.7','0.9'});
axis([0.5 5.5 0 1]);
yticks(0:0.25:1);
set(gca,'tickdir','out');
set(gca,'box','off', 'TitleHorizontalAlignment','center');
set(gca,'FontSize',14);

xlabel('Difficulties (Anti probability)','fontsize',12);
ylabel('Pr(Switch)','fontsize',12);
title('2-fold CV: CBM prediction','fontsize',18,'FontWeight','normal');

legend([h_line(1), h_line(2)], {'1B-Er CV model','2B-Er CV model'}, ...
    'Location','best');
legend boxoff;


%% =========================================================
% Logistic CV plot and CBM-vs-logistic held-out NLL comparison
% =========================================================

fig_lr = figure('Color','w','Position',[150 150 900 300]);
subplot(1,2,1); hold on;

LR_obs_subj   = squeeze(nanmean(LR_Pr_Sw_obs, 2));      % subj x T x difficulty
LR_model_subj = squeeze(nanmean(LR_Pr_Sw_model, 2));    % subj x T x difficulty

clear h_line_lr h_fill_lr
for iT = 1:2

    y_obs = squeeze(LR_obs_subj(:,iT,:));
    y_mod = squeeze(LR_model_subj(:,iT,:));

    obs_mean = nanmean(y_obs, 1);
    obs_sem  = nanstd(y_obs, [], 1) ./ sqrt(sum(~isnan(y_obs), 1));

    mod_mean = nanmean(y_mod, 1);
    mod_sem  = nanstd(y_mod, [], 1) ./ sqrt(sum(~isnan(y_mod), 1));

    % observed held-out data; these should match the CBM observed points
    % because both models use the same switchInput and the same train/test split.
    errorbar(xvalue, obs_mean, obs_sem, 'ko', ...
        'LineWidth', 2, ...
        'MarkerFaceColor', c_cv{iT}, ...
        'MarkerSize', 9);
    hold on;

    % logistic held-out prediction
    xx = [xvalue, fliplr(xvalue)];
    yy = [mod_mean + mod_sem, fliplr(mod_mean - mod_sem)];

    h_fill_lr(iT) = fill(xx, yy, c_cv{iT}, ...
        'FaceAlpha', 0.25, ...
        'EdgeColor', 'none');
    hold on;

    h_line_lr(iT) = plot(xvalue, mod_mean, '-', ...
        'Color', c_cv{iT}, ...
        'LineWidth', 2.5);
    hold on;

end

set(gca,'XTick',1:5);
set(gca,'XTickLabel',{'0.1','0.3','0.5','0.7','0.9'});
axis([0.5 5.5 0 1]);
yticks(0:0.25:1);
set(gca,'tickdir','out');
set(gca,'box','off');
set(gca,'FontSize',14);

xlabel('Difficulties (Anti probability)','fontsize',12);
ylabel('Pr(Switch)','fontsize',12);
text(-0.16, 1.08, 'A','Units','normalized', ...
    'FontSize',18, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');
title('2-fold CV: subject-wise logistic','fontsize',14,'FontWeight','normal');

legend([h_line_lr(1), h_line_lr(2)], ...
    {'E1 logistic','E2 logistic'}, ...
    'Location','best');
legend boxoff;


% ---- Compare CBM vs logistic benchmark using held-out NLL ----
CBM_NLL_subj = nanmean(CV_test_NLL, 2);   % one value per subject
LR_NLL_subj  = nanmean(LR_test_NLL, 2);   % one value per subject

valid = ~isnan(CBM_NLL_subj) & ~isnan(LR_NLL_subj);

CBM_NLL_subj = CBM_NLL_subj(valid);
LR_NLL_subj  = LR_NLL_subj(valid);

% paired t-test
[~, p_t, ~, stats_t] = ttest(CBM_NLL_subj, LR_NLL_subj);

% non-parametric paired test, safer if non-normal
[p_w, ~, stats_w] = signrank(CBM_NLL_subj, LR_NLL_subj);

fprintf('\n=== Held-out NLL comparison: CBM vs logistic, same CV split ===\n');
fprintf('CBM mean NLL      = %.4f +/- %.4f SEM\n', ...
    mean(CBM_NLL_subj), std(CBM_NLL_subj)/sqrt(numel(CBM_NLL_subj)));
fprintf('Logistic mean NLL = %.4f +/- %.4f SEM\n', ...
    mean(LR_NLL_subj), std(LR_NLL_subj)/sqrt(numel(LR_NLL_subj)));

fprintf('Paired t-test: t(%d) = %.3f, p = %.4g\n', ...
    stats_t.df, stats_t.tstat, p_t);

fprintf('Wilcoxon signed-rank: z = %.3f, p = %.4g\n', ...
    stats_w.zval, p_w);

% difference: positive means CBM better
Delta_NLL = LR_NLL_subj - CBM_NLL_subj;

fprintf('Mean Delta_NLL = Logistic - CBM = %.4f\n', mean(Delta_NLL));
fprintf('Positive Delta_NLL means CBM predicts held-out switching better.\n');

subplot(1,2,2); hold on;
x1 = ones(size(CBM_NLL_subj));
x2 = 2 * ones(size(LR_NLL_subj));

for i = 1:numel(CBM_NLL_subj)
    plot([1 2], [CBM_NLL_subj(i), LR_NLL_subj(i)], ...
        '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1);
    hold on;
end

scatter(x1, CBM_NLL_subj, 55, 'filled');
scatter(x2, LR_NLL_subj, 55, 'filled');

errorbar(1, mean(CBM_NLL_subj), ...
    std(CBM_NLL_subj)/sqrt(numel(CBM_NLL_subj)), ...
    'k', 'LineWidth', 2);

errorbar(2, mean(LR_NLL_subj), ...
    std(LR_NLL_subj)/sqrt(numel(LR_NLL_subj)), ...
    'k', 'LineWidth', 2);

set(gca,'XTick',[1 2], 'XTickLabel',{'CBM','Logistic'});
xlim([0.5 2.5]);
ylabel('Held-out negative log likelihood');
text(-0.16, 1.08, 'B','Units','normalized', ...
    'FontSize',18, ...
    'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'VerticalAlignment','top');
title(sprintf('Held-out prediction, p = %.3g', p_w), ...
    'fontsize',14,'FontWeight','normal');
box off;
set(gca,'TickDir','out');
hold off;

if isSave
    figure(fig_lr);
    print(gcf, 'FigS_Logistic_2foldCV_sameSplit', '-dsvg', '-painters');
end

% Return to the original figure so parameter recovery plots are added to subplots 1-3.
figure(fig_cvrec);



%% =========================================================
% Fig. S5C-like parameter recovery
% Random true parameters -> simulate all subjects -> average recovered params
% =========================================================

nSim = 500;   % number of random true parameter sets

true_sigma = nan(nSim,1);
true_alpha = nan(nSim,1);
true_pam3  = nan(nSim,1);

rec_sigma  = nan(nSim,1);
rec_alpha  = nan(nSim,1);
rec_pam3   = nan(nSim,1);

rec_sigma_subj = nan(nSim,length(id));
rec_alpha_subj = nan(nSim,length(id));
rec_pam3_subj  = nan(nSim,length(id));

parameter_strings_sw = { ...
    'sigma_switch = switchPam(1);', ...
    'alpha_transition = switchPam(2);', ...
    'pam3 = switchPam(3);'};

for iSim = 1:nSim

    fprintf('Parameter recovery simulation %d / %d\n', iSim, nSim);

    % 1. Generate one true parameter set
    truePam = [ ...
        0.1  + (10   - 0.1)  * rand, ...   % sigma_switch
        0.01 + (0.99 - 0.01) * rand, ...   % alpha_transition
        0.01 + (0.99 - 0.01) * rand];      % pam3

    true_sigma(iSim) = truePam(1);
    true_alpha(iSim) = truePam(2);
    true_pam3(iSim)  = truePam(3);

    % 2. Use this same truePam to simulate and recover all subjects
    for iSubject = 1:length(id)

        switchInput0 = input{iSubject};

        temp_parameters_as_input = model.psychometric.subj(iSubject).EstimatedParameters(2,:);
        temp_parameter_strings   = model.psychometric.subj(iSubject).parameter_strings;
        temp_SubjObj_flag        = model.psychometric.subj(iSubject).SubjObj_flag;

        % expected accuracy benchmark
        clear expectedAccuracy_Benchmark
        flag = 0;

        for iCue = 0:1
            for tDev_i = 1:length(InputAll(iSubject).DevValues)

                flag = flag + 1;

                tmp.RuleChoice = iCue;
                tmp.tDev = InputAll(iSubject).DevValues(tDev_i);

                expectedAccuracy_Benchmark.RuleChoice(flag,1) = iCue;
                expectedAccuracy_Benchmark.tDev(flag,1) = InputAll(iSubject).DevValues(tDev_i);
                expectedAccuracy_Benchmark.expectedAccuracy(flag,1) = expected_Accuracy( ...
                    tmp, ...
                    temp_parameters_as_input, ...
                    temp_parameter_strings, ...
                    temp_SubjObj_flag);
            end
        end

        % generate Pr(Switch)
        [p_sw_true, ~, ~, ~, ~, ~] = pr_switch_func( ...
            switchInput0, ...
            truePam, ...
            parameter_strings_sw, ...
            temp_parameters_as_input, ...
            temp_parameter_strings, ...
            temp_SubjObj_flag, ...
            expectedAccuracy_Benchmark);

        p_sw_true = min(max(p_sw_true(:), 1e-4), 1-1e-4);

        % simulate switch choices
        simInput = switchInput0;

        for k = 1:length(simInput)
            simInput(k).SW = rand < p_sw_true(k);
        end

        % recover parameters for this subject
        Set_param_variable = [1, 0.5, 0.5];

        [~, est_rec, ~, ~, ~] = SwitchParameterOptimization( ...
            simInput, ...
            temp_parameters_as_input, ...
            temp_parameter_strings, ...
            temp_SubjObj_flag, ...
            InputAll(iSubject), ...
            Set_param_variable);

        rec_sigma_subj(iSim,iSubject) = est_rec.sigma_switch_estimated;
        rec_alpha_subj(iSim,iSubject) = est_rec.alpha_transition;
        rec_pam3_subj(iSim,iSubject)  = est_rec.pam3_estimated;

    end

    % 3. Average recovered parameters across subjects
    rec_sigma(iSim) = mean(rec_sigma_subj(iSim,:), 'omitnan');
    rec_alpha(iSim) = mean(rec_alpha_subj(iSim,:), 'omitnan');
    rec_pam3(iSim)  = mean(rec_pam3_subj(iSim,:),  'omitnan');

end

save('model_res_switch_parameter_recovery_allsubj.mat', ...
    'true_sigma', 'true_alpha', 'true_pam3', ...
    'rec_sigma', 'rec_alpha', 'rec_pam3', ...
    'rec_sigma_subj', 'rec_alpha_subj', 'rec_pam3_subj');



%% =========================================================
% Plot Fig. S5C-like parameter recovery
% =========================================================
load model_res_switch_parameter_recovery_allsubj.mat

subplot(1,4,1); hold on;
x = true_sigma(:);
y = rec_sigma(:);
valid = ~isnan(x) & ~isnan(y);

[r_sigma, p_sigma] = corr(x(valid), y(valid), 'Type', 'Pearson');

scatter(x(valid), y(valid), 28, 'filled', 'MarkerFaceAlpha', 0.55);
plot([0.1 10], [0.1 10], 'k--', 'LineWidth', 1.5);
xlabel('True \sigma');
ylabel('Recovered \sigma');
title(sprintf('\\sigma, r = %.2f, p = %.2g', r_sigma, p_sigma));
axis square; box off;
xlim([0 10]); ylim([0 10]);


subplot(1,4,2); hold on;
x = true_alpha(:);
y = rec_alpha(:);
valid = ~isnan(x) & ~isnan(y);

[r_alpha, p_alpha] = corr(x(valid), y(valid), 'Type', 'Pearson');

scatter(x(valid), y(valid), 28, 'filled', 'MarkerFaceAlpha', 0.55);
plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
xlabel('True \alpha');
ylabel('Recovered \alpha');
title(sprintf('\\alpha, r = %.2f, p = %.2g', r_alpha, p_alpha));
axis square; box off;
xlim([0 1]); ylim([0 1]);


subplot(1,4,3); hold on;
x = true_pam3(:);
y = rec_pam3(:);
valid = ~isnan(x) & ~isnan(y);

[r_pam3, p_pam3] = corr(x(valid), y(valid), 'Type', 'Pearson');

scatter(x(valid), y(valid), 28, 'filled', 'MarkerFaceAlpha', 0.55);
plot([0 1], [0 1], 'k--', 'LineWidth', 1.5);
xlabel('True pam3');
ylabel('Recovered pam3');
title(sprintf('pam3, r = %.2f, p = %.2g', r_pam3, p_pam3));
axis square; box off;
xlim([0 1]); ylim([0 1]);


% subplot(1,4,4); hold on;
% axis off;
% 
% text(0, 0.75, sprintf('\\sigma: r = %.2f, p = %.2g', r_sigma, p_sigma), 'FontSize', 12);
% text(0, 0.55, sprintf('\\alpha: r = %.2f, p = %.2g', r_alpha, p_alpha), 'FontSize', 12);
% text(0, 0.35, sprintf('pam3: r = %.2f, p = %.2g', r_pam3, p_pam3), 'FontSize', 12);
% 
% sgtitle('Parameter recovery');

if isSave
    print(gcf, 'FigS_parameter_recovery_random', '-dsvg', '-painters');
end