function [output,estimate_parameters, MachineSimulation,Input,MLE] = SwitchParameterOptimization(switchInput, psych_parameters, psych_parameter_strings, psych_SubjObj_flag, Input, Set_param_variable);       % pr_switch_model(iRuleChoice, itd, T_numOfBackError )

% Setting the initialized values of optimizer
options = optimset('fminsearch');
options.Display = 'iter'; % 'off'
options.Iter = 100000;
options.TolFun = 1e-6;
options.TolX = 1e-6;
options.MaxFunEvals = 100000;

% % % % % % thrs          = [.01, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, .99];%********************************
thrs          = [0.1 0.3 0.5 0.7 0.9];%[0.01 1/4 1/2 3/4 0.99];  %[1/6 2/6 3/6 4/6 5/6];%[0.0675, 5/16, 1/2, 11/16, .9325];%********************************
devlevel = (.5-(abs(thrs-.5)))/0.1;%( (0.5-abs(thrs-.5)));%[1:1:5,4:-1:1];%0.5-abs(thrs-.5);%abs(thrs-.5);%[5:-1:1,2:1:5];%(0.5-abs(thrs-.5))/.01;%
% % % devlevel = [1 3.42 7.28 3.42 1]; %sigmoid function, k=0.35, x0=5
mOfT = 3;

alpha_transition = 0.5;%0.14; %21/150
% ******#2 by xu 20231005, for defining lamada as free parameter;

% switch_bound = 1;

switchPam.sigma_switch_estimated = Set_param_variable(1);
%switchPam.switch_bound = Set_param_variable(2);
switchPam.scale = Set_param_variable(2);
switchPam.pam3  = Set_param_variable(3);%1; % Set_param_variable(3); %？？？？？？****** alpha??#############

switchPam.alpha_transition = alpha_transition;
% ****** by xu 20231005, for defining lamada as free parameter;

switchPam.beta1 = 1; %unused latter******
switchPam.beta2 = -10;%unused latter******

StartPointInitializedValues = [switchPam.sigma_switch_estimated, switchPam.alpha_transition,switchPam.pam3];% #############,switchPam.scale , switchPam.scale   ,switchPam.beta1, switchPam.beta2
parameter_strings = { 'sigma_switch = switchPam(1)','alpha_transition = switchPam(2)','pam3=switchPam(3)'};% #############, 'scalet = switchPam(3);' , 'scalet = switchPam(3);'  ,'beta1 = switchPam(3);', 'beta2 = switchPam(4);'
lb = [ 0.1 0.01 0.01]; %0.1 -20#############lb = [ 0.1 0.01 0.01];
ub = [ 10 0.99 0.99]; % 2 -0.1#############ub = [ 10 0.99 0.99];
%******#1 by xu 20231005, and change the 'alpha_transition' in pr_switch_func.m


% % % % % % StartPointInitializedValues = [  switchPam.sigma_switch_estimated ];% switchPam.alpha_transition,,switchPam.scale , switchPam.scale   ,switchPam.beta1, switchPam.beta2
% % % % % % parameter_strings = { 'sigma_switch = switchPam(1);' };%  'alpha_transition = switchPam(1);',, 'scalet = switchPam(3);' , 'scalet = switchPam(3);'  ,'beta1 = switchPam(3);', 'beta2 = switchPam(4);'
% % % % % % 
% % % % % % lb = [ 0 ]; %0.1 -20
% % % % % % ub = [10]; % 2 -0.1
switch_bound_str = '1/( 1 + exp(-1* log10(T*devlevel(find(dev==thrs))/.1) *switch_bound))';
%% expectedAccuracy_Benchmark
flag = 0;
for iCue = 0: 1
    for tDev = 1:5 %9 %********************************
        flag = flag+1;
        expectedAccuracy_BenchmarkInput.RuleChoice = iCue;
        expectedAccuracy_BenchmarkInput.tDev = Input.DevValues(tDev);
        expectedAccuracy_Benchmark.RuleChoice(flag,1) = iCue;
        expectedAccuracy_Benchmark.tDev(flag,1) = Input.DevValues(tDev);
        expectedAccuracy_Benchmark.expectedAccuracy(flag,1) = expected_Accuracy(expectedAccuracy_BenchmarkInput, psych_parameters, psych_parameter_strings, psych_SubjObj_flag);
    end
end

% fit the psychometric model: %******
[logLikeSwitch, FVAL] = fminsearchbnd(@(switchPam) logLike_of_pr_switch(switchInput, switchPam, parameter_strings, psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark, mOfT), StartPointInitializedValues, lb,ub,options);
%[logLikeSwitch, FVAL] = fminsearch(@(switchPam) logLike_of_pr_switch(switchInput, alpha_transition, switchPam, parameter_strings, psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark, mOfT), StartPointInitializedValues, options);
%[logLikeSwitch, FVAL] = fminsearch(@(switchPam) logLike_of_pr_switch(switchInput, alpha_transition, switchPam(1), switchPam(2), psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark, mOfT), Set_param_variable, options);

MLE = -FVAL;

% % % % % % estimate_parameters.sigma_switch_estimated = logLikeSwitch(1);%****** by xu 20231020
estimate_parameters.sigma_switch_estimated = logLikeSwitch(1); %******#refer to: parameter_strings = { 'sigma_switch = switchPam(1)','alpha_transition = switchPam(2)' };
estimate_parameters.pam3_estimated = logLikeSwitch(3);%1;%logLikeSwitch(3);#############
estimate_parameters.switch_bound = 1;
switch_bound = estimate_parameters.switch_bound;
estimate_parameters.scale = 1; %%%%%% 1;%logLikeSwitch(3);
scalet = estimate_parameters.scale;

if logLikeSwitch(2)<1 % ******#1 by xu 20231005, for defining lamada as free parameter;Inp_Pr_Sw
    alpha_transition = logLikeSwitch(2);% ******#1 by xu 20231005, for defining lamada as free parameter;
else% ******#1 by xu 20231005, for defining lamada as free parameter;
    alpha_transition = 1/logLikeSwitch(2);% ******#1 by xu 20231005, for defining lamada as free parameter;
end% ******#1 by xu 20231005, for defining lamada as free parameter;



estimate_parameters.alpha_transition = alpha_transition;

% beta1 = logLikeSwitch(3);
% beta2 = logLikeSwitch(4);
% estimate_parameters.beta1 = beta1;
% estimate_parameters.beta2 = beta2;

% ***************************#？？？？？？
[Output_pr_of_switch, Output_tDev_lastOne, Output_RuleChoice_lastOne, Output_T, Output_SW, mu_switch_estimated] = pr_switch_func(switchInput, logLikeSwitch, parameter_strings, psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark);
%[Output_pr_of_switch, Output_tDev_lastOne, Output_RuleChoice_lastOne, Output_T, Output_SW, mu_switch_estimated] = pr_switch_func(switchInput,  0.2, logLikeSwitch(1), logLikeSwitch(2), psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark);

%% fitting
unitsq_sgm = @(x,zeta)( x.^zeta ./ (x.^zeta+(1-x).^zeta) ) ;% un-used
% A: sigma = sigma * T
% % % % % % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(T) ) ) )  )))  ) ;% adpated by wenshan 20220511
% ******#1 by xu 20231006

% B: sigma = sigma * T^2
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( (T) ) ) )  )))  ) ;
% C: sigma = sigma * (scale*devlevel+T)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt( scalet * devlevel(find(dev==thrs)) + T ) ) ) )  )))  *  (scalet>0) *  (scalet<=5) * (sigma_switch<=4) ) ;%* (sigma_switch>1)
% D: sigma = sigma * (scale*devlevel+T)^2
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(     ( scalet * devlevel(find(dev==thrs)) + T ) ) ) )  ))) *  (scalet>0) *  (scalet<=5) * (sigma_switch<=4) ) ;% adpated by wenshan 20220421  *  (scalet>0) *  (scalet<1)* (sigma_switch>0)
% E: sigma = sigma * (scale*devlevel*T)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt( scalet * devlevel(find(dev==thrs)) * T ) ) ) )  )))  ) ;
% F: sigma = sigma * (scale*devlevel*T)^2
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(     ( scalet * devlevel(find(dev==thrs)) * T ) ) ) )  ))) ) ;% adpated by wenshan 20220421  *  (scalet>0) *  (scalet<1)* (sigma_switch>0)
% G: sigma = sigma * (devlevel*log(T))
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt( scalet * exp(devlevel(find(dev==thrs))) * (T) ) ) ) )  )))  ) ;
% F: sigma = sigma * ((exp/T)^devlevel)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt( scalet * (exp(1)/T)^(devlevel(find(dev==thrs))) ) ) ) )  )))  * (sigma_switch<=4) ) ;
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  (scalet/T)^(devlevel(find(dev==thrs))) ) ) ) )  )))  * (sigma_switch<=4) ) ;
% I: sigma = sigma * (T/devlevel)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  scalet*log10(T / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  scalet*T + (T / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  (  scalet*T + (1 / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  scalet*((-1) * sqrt(T) * log(devlevel(find(dev==thrs))) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  ((-1) * (T^scalet) * log(devlevel(find(dev==thrs))) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  scalet* ( sqrt(T) * (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  (  scalet* ( sqrt(T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  scalet* ( (T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  (  scalet* ( (T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  scalet* ( log(T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  ( log(scalet*T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  scalet* ( (T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  ( scalet*log(T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  (  scalet^(T / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)

pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt( scalet*  (T / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)

% % % % % %  pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt( scalet*  (T / (2*devlevel(find(dev==thrs))) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% ******#2 by xu 20231006

% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt( scalet*  (T / (2*dev) ) ) ) ) )  )))   ) ;
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (1/( 1 + exp(-1* log(devlevel(find(dev==thrs))/.01) *x) )-mu_switch)/(sqrt(2)*(sigma_switch*(  ( scalet*  (T / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  log(scalet*T/(devlevel(find(dev==thrs))) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( (  log10(scalet*T/ dev ) ) ) ) )  )))   ) ;
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  scalet*T/ dev ) ) ) ) )  )))   ) ;
% pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( (  (beta1*T + beta2* dev ) ) ) ) )  )))   ) ;

% A: sigma = sigma * T
% % % % % % cor_str = 'cor = sqrt(T);';
% ******#1 by xu 20231006

% B: sigma = sigma * T^2
% cor =     (T);
% C: sigma = sigma * (scale*devlevel+T)
% cor = sqrt( scalet * devlevel(find(dev==thrs)) + T );
% D: sigma = sigma * (scale*devlevel+T)^2
% cor =     ( scalet * devlevel(find(dev==thrs)) + T );
% E: sigma = sigma * (scale*devlevel*T)
% cor = sqrt( scalet * devlevel(find(dev==thrs)) * T );
% F: sigma = sigma * (scale*devlevel*T)^2
% cor =     ( scalet * devlevel(find(dev==thrs)) * T );
% G: sigma = sigma * (devlevel*log(T))
% cor = sqrt( scalet * exp(devlevel(find(dev==thrs))) * (T) );
% F: sigma = sigma * ((exp/T)^devlevel)
% cor = sqrt( scalet * (exp(1)/T)^(devlevel(find(dev==thrs))) );
% cor = sqrt(  (scalet/T)^(devlevel(find(dev==thrs))) );
% I: sigma = sigma * (T/devlevel)
% cor =  sqrt( scalet*log10(T / devlevel(find(dev==thrs)) ) );
% cor = sqrt( scalet*T + (T / devlevel(find(dev==thrs)) ) );
% cor = ( scalet*T + (1 / devlevel(find(dev==thrs)) ) );
% cor =    sqrt(scalet*((-1) * sqrt(T) * log(devlevel(find(dev==thrs))) ) );
% cor =    sqrt( scalet*((-1) * (T^scalet) * log(devlevel(find(dev==thrs))) ) );
% cor =    sqrt( scalet* ( sqrt(T) - log(devlevel(find(dev==thrs))) ) );
% cor =    sqrt( scalet* ( sqrt(T) * ( - log(devlevel(find(dev==thrs))) ) ) );
% cor =    ( scalet* ( sqrt(T) + ( - log(devlevel(find(dev==thrs))) ) ) );
% cor =    ( scalet* ( (T) + ( - log(devlevel(find(dev==thrs))) ) ) );
% cor = sqrt( scalet* ( log(T) + ( - log(devlevel(find(dev==thrs))) ) ) );
% cor = sqrt( ( log(scalet* T) + ( - log(devlevel(find(dev==thrs))) ) ) );
% cor = sqrt( ( scalet* log(T) + ( - log(devlevel(find(dev==thrs))) ) ) );
% cor = ( scalet^(T / devlevel(find(dev==thrs)) ) );

cor_str = 'cor =  sqrt( scalet* (T / devlevel(find(dev==thrs)) ) );';

% % % % % %  cor_str = 'cor =  sqrt( scalet* (T / (2*devlevel(find(dev==thrs))) ) );';% ******
% ******#2 by xu 20231006

% cor_str = 'cor =  sqrt( scalet* (T / (2*dev) ) );';
% cor_str = 'cor = sqrt( log( scalet* (T /(devlevel(find(dev==thrs)))) ) );';
% cor_str = 'cor = ( log10( scalet* ( T / dev ) ) );';
% cor_str = 'cor = ( ( beta1*  T +beta2* dev  ) );';
%% recover
MachineSimulation.StateExp = Input.StateExp;
% MachineSimulation.TF = Input.TF; %changed by wenshan 20220518
MachineSimulation.tDev = Input.tDev;
%MachineSimulation.simRuleChoice = MachineSimulation.RuleChoice;
MachineSimulation.Cue = Input.Cue;

fstIdx = find(MachineSimulation.StateExp ==1); %151:300
for iTrial = fstIdx(1):fstIdx(1)+length(fstIdx)-1 %151:(151+150-1)
    if iTrial == fstIdx(1) %fstIdx(1)=151, the first trial
        MachineSimulation.RuleChoice(iTrial,1) = Input.RuleChoice(iTrial);%the first ruleChoice of MS is equal to Input,
    else  %otherwise is equal to：...SW and ruleChoice of (i-1)trial;xor(1,0)=1, xor(1,1)=0; xor(0,0)=0, xor(0,1)=1.
         %if sw(i-1)=1, ruleChoice of i will change; if sw(i-1)==0, ruleChoice of i will not change;
         %the most important is the SW of (i-1) trial
        MachineSimulation.RuleChoice(iTrial,1) = double(xor( MachineSimulation.SW(iTrial-1), MachineSimulation.RuleChoice(iTrial-1)));
    end
    MachineSimulation.TF(iTrial,1) = MachineSimulation.RuleChoice(iTrial) == Input.Cue(iTrial); % the feedback is determined by actual rule of Input
    
    
    
    
    if MachineSimulation.StateExp(iTrial) ==1 %select the trials in infer experiment
        if iTrial == fstIdx(1) %fstIdx(1)=151, the first trial
            Input.Nback(iTrial,1) = double(~Input.TF(iTrial));%the Nback of 151's trial = the reversal feedback of 151's；(~1)=0
            MachineSimulation.Nback(iTrial,1) = double(~MachineSimulation.TF(iTrial));%the Nback of 151's trial = the reversal feedback of 151's；(~1)=0
        else %otherwise(152-300's trial) is equal to：
            Input.Nback(iTrial,1) = (double(~Input.TF(iTrial)) * (Input.Nback(iTrial-1) + 1))*  (~Input.SW(iTrial-1))+  (~Input.TF(iTrial)) * Input.SW(iTrial-1);%changed by wenshan 20220518
            %****************************************************************
            MachineSimulation.Nback(iTrial,1) = double(~MachineSimulation.TF(iTrial)) * (MachineSimulation.Nback(iTrial-1) + 1) * (~MachineSimulation.SW(iTrial-1)) + (~MachineSimulation.TF(iTrial)) * MachineSimulation.SW(iTrial-1);%changed by wenshan 20220518
        end
        
        if  MachineSimulation.Nback(iTrial)==0
            iCT_index = find( (expectedAccuracy_Benchmark.RuleChoice == MachineSimulation.RuleChoice(iTrial) ) .* (expectedAccuracy_Benchmark.tDev == MachineSimulation.tDev(iTrial) )  );
            A_m = expectedAccuracy_Benchmark.expectedAccuracy(iCT_index);
            MachineSimulation.expectedAccuracy(iTrial,1) = A_m;
            MachineSimulation.mu_switch_estimated(iTrial,1) = 0;
            MachineSimulation.pr_of_switch(iTrial,1) = 0;
            MachineSimulation.SW(iTrial:iTrial+1,1) = [0 , 0];%********************************
        else
            clear A H
            for iT = 1:MachineSimulation.Nback(iTrial)
                idx = iTrial-MachineSimulation.Nback(iTrial)+iT;
                iCT_index = find( (expectedAccuracy_Benchmark.RuleChoice == MachineSimulation.RuleChoice(idx) ) .* (expectedAccuracy_Benchmark.tDev == MachineSimulation.tDev(idx) )  );
                A(iT) = expectedAccuracy_Benchmark.expectedAccuracy(iCT_index);
                H(iT) = alpha_transition;
            end
            T = MachineSimulation.Nback(iTrial);
            dev = MachineSimulation.tDev(iTrial);
            syntheticInput.RuleChoice = MachineSimulation.RuleChoice(iTrial);
            syntheticInput.Cue = MachineSimulation.Cue(iTrial);
            syntheticInput.tDev = MachineSimulation.tDev(iTrial);
            %syntheticInput.DevValues = MachineSimulation.DevValues(iTrial);%******
            %dev = 0.5-abs(pr_anti_c_td(syntheticInput, psych_parameters, psych_parameter_strings, psych_SubjObj_flag)-.5);%******
            eval(cor_str);
            switch_bound_t = 1;% %****** eval(switch_bound_str); %
            MachineSimulation.mu_switch_estimated(iTrial,1) = abs(switch_bound_t - sqrt(2) * estimate_parameters.sigma_switch_estimated * (cor) * erfinv( (1 - PO(A, H, MachineSimulation.Nback(iTrial))) ./ (1 + PO(A, H, MachineSimulation.Nback(iTrial)))  ) ) * estimate_parameters.pam3_estimated;
            % MachineSimulation.mu_switch_estimated(iTrial,1) =  abs(switch_bound_t - sqrt(2) * estimate_parameters.sigma_switch_estimated * (cor) * erfinv( (1 - PO(A, H, MachineSimulation.Nback(iTrial))) ./ (1 + PO(A, H, MachineSimulation.Nback(iTrial)))  ) ) * unitsq_sgm( estimate_parameters.pam3_estimated,log10(10*dev));
            % MachineSimulation.mu_switch_estimated(iTrial,1) =  abs(switch_bound_t - sqrt(2) * estimate_parameters.sigma_switch_estimated * (cor) * erfinv( (1 - PO(A, H, MachineSimulation.Nback(iTrial))) ./ (1 + PO(A, H, MachineSimulation.Nback(iTrial)))  ) ) * estimate_parameters.pam3_estimated * exp(-dev);
            % MachineSimulation.mu_switch_estimated(iTrial,1) =  abs(switch_bound_t - sqrt(2) * estimate_parameters.sigma_switch_estimated * (cor) * erfinv( (1 - PO(A, H, MachineSimulation.Nback(iTrial))) ./ (1 + PO(A, H, MachineSimulation.Nback(iTrial)))  ) ) * estimate_parameters.pam3_estimated * cor;
            MachineSimulation.pr_of_switch(iTrial,1) = pr_of_switch(switch_bound_t, MachineSimulation.mu_switch_estimated(iTrial) , estimate_parameters.sigma_switch_estimated,scalet, MachineSimulation.Nback(iTrial),dev);
            % MachineSimulation.pr_of_switch(iTrial,1) = unitsq_sgm(MachineSimulation.pr_of_switch(iTrial,1),log(1/dev));
            % MachineSimulation.pr_of_switch(iTrial,1) = unitsq_sgm(MachineSimulation.pr_of_switch(iTrial,1),scalet*log(T)-log(2*dev));
            % MachineSimulation.pr_of_switch(iTrial) = 1/(1+exp(-1/dev*MachineSimulation.pr_of_switch(iTrial)));% added by wenshan 20020516
            
            if iTrial+1 <= fstIdx(1)+length(fstIdx)-1
            MachineSimulation.SW(iTrial:iTrial+1,1) = [binornd(1,MachineSimulation.pr_of_switch(iTrial) ), 0];
            % MachineSimulation.RuleChoice(iTrial+1) = double(xor(MachineSimulation.pr_of_switch(iTrial)>=.5, MachineSimulation.RuleChoice(iTrial)));
            % MachineSimulation.RuleChoice(iTrial+1) = double(xor(MachineSimulation.mu_switch_estimated(iTrial)>=switch_bound, MachineSimulation.RuleChoice(iTrial)));
            % MachineSimulation.RuleChoice(iTrial+1) = (MachineSimulation.pr_of_switch(iTrial)/(1-MachineSimulation.pr_of_switch(iTrial)))>=1; %changed by wenshan 20220514
            end
            MachineSimulation.expectedAccuracy(iTrial,1) = A(end);
        end
        
        if  Input.Nback(iTrial)==0
            iCT_index = find( (expectedAccuracy_Benchmark.RuleChoice == Input.RuleChoice(iTrial) ) .* (expectedAccuracy_Benchmark.tDev == Input.tDev(iTrial) )  );
            A_m = expectedAccuracy_Benchmark.expectedAccuracy(iCT_index);
            Input.expectedAccuracy(iTrial,1) = A_m;
            Input.alpha_transition(iTrial,1) = alpha_transition;
            Input.mu_switch_estimated(iTrial,1) = 0;
            Input.pr_of_switch(iTrial,1) = 0;
        else
            clear A H
            for iT = 1:Input.Nback(iTrial,1)
                idx = iTrial-Input.Nback(iTrial,1)+iT;
                iCT_index = find( (expectedAccuracy_Benchmark.RuleChoice == Input.RuleChoice(idx) ) .* (expectedAccuracy_Benchmark.tDev == Input.tDev(idx) )  );
                A(iT) = expectedAccuracy_Benchmark.expectedAccuracy(iCT_index);
                H(iT) = alpha_transition;
            end
            T = Input.Nback(iTrial,1);
            dev = Input.tDev(iTrial,1);
            syntheticInput.RuleChoice = Input.RuleChoice(iTrial);
            syntheticInput.Cue = Input.Cue(iTrial);
            syntheticInput.tDev = Input.tDev(iTrial);
            %syntheticInput.DevValues = Input.DevValues(iTrial);%******
            %dev = 0.5-abs(pr_anti_c_td(syntheticInput, psych_parameters, psych_parameter_strings, psych_SubjObj_flag)-.5);%******
            eval(cor_str);
            switch_bound_t = 1;% %****** eval(switch_bound_str); %
            Input.mu_switch_estimated(iTrial,1) = abs(switch_bound_t - sqrt(2) * estimate_parameters.sigma_switch_estimated * (cor) * erfinv( (1 - PO(A, H, Input.Nback(iTrial))) ./ (1 + PO(A, H, Input.Nback(iTrial)))  ) ) * estimate_parameters.pam3_estimated;
            % Input.mu_switch_estimated(iTrial,1) = abs(switch_bound_t - sqrt(2) * estimate_parameters.sigma_switch_estimated * (cor) * erfinv( (1 - PO(A, H, Input.Nback(iTrial))) ./ (1 + PO(A, H, Input.Nback(iTrial)))  ) ) * unitsq_sgm(estimate_parameters.pam3_estimated, log10(10*dev));
            % Input.mu_switch_estimated(iTrial,1) = abs(switch_bound_t - sqrt(2) * estimate_parameters.sigma_switch_estimated * (cor) * erfinv( (1 - PO(A, H, Input.Nback(iTrial))) ./ (1 + PO(A, H, Input.Nback(iTrial)))  ) ) * estimate_parameters.pam3_estimated * exp(-dev);
            % Input.mu_switch_estimated(iTrial,1) = abs(switch_bound_t - sqrt(2) * estimate_parameters.sigma_switch_estimated * (cor) * erfinv( (1 - PO(A, H, Input.Nback(iTrial))) ./ (1 + PO(A, H, Input.Nback(iTrial)))  ) ) * estimate_parameters.pam3_estimated * cor;
            Input.pr_of_switch(iTrial,1) = pr_of_switch(switch_bound_t, Input.mu_switch_estimated(iTrial) , estimate_parameters.sigma_switch_estimated, scalet,Input.Nback(iTrial),dev);
            % Input.pr_of_switch(iTrial,1) = unitsq_sgm(Input.pr_of_switch(iTrial,1),log(1/dev));
            % Input.pr_of_switch(iTrial,1) = unitsq_sgm(Input.pr_of_switch(iTrial,1),log10((T)/(dev)));
            % Input.pr_of_switch(iTrial,1) = unitsq_sgm(Input.pr_of_switch(iTrial,1),scalet*log(T)-log(2*dev));
            if iTrial < fstIdx(1)+length(fstIdx)-1
                Input.SW(iTrial,1) = binornd(1,Input.pr_of_switch(iTrial) );
            else
                Input.SW(iTrial,1) = 0; % the laste one
            end
            Input.expectedAccuracy(iTrial,1) = A(end);
            Input.alpha_transition(iTrial,1) = H(end);
            
        end
        
    end
    
    
    
    
end

% conditions_1Back_index = find( (Input.StateExp(1:end-2) == 1) .* (Input.TF(1:end-2)==1) .* (Input.TF(2:end-1)==0) ) +1; % RW, ER  => {1B-Er}
% conditions_2Back_index = find( (Input.StateExp(1:end-3) == 1) .* (Input.TF(1:end-3)==1) .* (Input.TF(2:end-2)==0) .* (Input.TF(3:end-1)==0) ) +2; % RW, ER, ER  => {2B-Er}
% Index_of_trials_forFit = [conditions_1Back_index; conditions_2Back_index];
% MachineSimulation.RuleChoice(Index_of_trials_forFit+1) = double(xor(Output_SW, Output_RuleChoice_lastOne));%%%%  
% 
output.Output_pr_of_switch = Output_pr_of_switch;
output.Output_tDev_lastOne = Output_tDev_lastOne;
output.Output_RuleChoice_lastOne = Output_RuleChoice_lastOne;
output.Output_T = Output_T;
output.Output_SW = Output_SW;
output.mu_switch_estimated = mu_switch_estimated;
