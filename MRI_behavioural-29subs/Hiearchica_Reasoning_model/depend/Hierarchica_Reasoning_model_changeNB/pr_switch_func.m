%% modeling the switching probability.
function [Output_pr_of_switch, Output_tDev_lastOne, Output_RuleChoice_lastOne, Output_T, Output_SW, mu_switch_estimated] = pr_switch_func(Input,  switchPam, parameter_strings, psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark);
% ******# commented by xu 20231005, for defining lamada as free parameter;

% function [Output_pr_of_switch, Output_tDev_lastOne, Output_RuleChoice_lastOne, Output_T, Output_SW, mu_switch_estimated] = pr_switch_func(Input, alpha_transition, sigma_switch, pam3,  psych_parameters, psych_parameter_strings, psych_SubjObj_flag, expectedAccuracy_Benchmark);
    % Setting the initialized values of optimizer
    options = optimset('fminsearch');
    options.Display = 'off'; % 'off'
    options.Iter = 1000000;
    options.TolFun = 1e-10;
    options.TolX = 1e-10;
    
    for iParameters = 1:length(parameter_strings)
        eval(parameter_strings{iParameters});
    end
    

% % % % % %     alpha_transition = 0.17;% ******#1 commented by xu 20231005, for defining lamada as free parameter;
    

%     pam3 = 1;%%%%%%1;% ******#1 commented by xu 20231005, for defining lamada as free parameter;
    scalet = 1;%%%%%%1
    %sigma_switch = 1;% ******#1 commented by xu 20231005, for defining lamada as free parameter;
    switch_bound = 1;
% % % % % % thrs          = [.01, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, .99];%********************************
    thrs = [0.1 0.3 0.5 0.7 0.9];%[0.01 1/4 1/2 3/4 0.99];  %[1/6 2/6 3/6 4/6 5/6];% [0.0675, 5/16, 1/2, 11/16, .9325];%********************************

    devlevel = (.5-(abs(thrs-.5)))/0.1;%( (0.5-abs(thrs-.5)));%[1:1:5,4:-1:1];%0.5-abs(thrs-.5);%abs(thrs-.5);%[5:-1:1,2:1:5];%(0.5-abs(thrs-.5))/.01;%

% some reminders:
    % The array "Input" includes only the error trials (sorted 1B, 2B, ...)
    % variables of switch_model.alpha_transition, switch_model.sigma_switch
    
    % How to contruct the Posterior odds
        % for 1 error: H(T) / [(1-H(T))(1-A(T))]   % arrayOfInput => (T):Array(1)
        % for 2 error: ( H(T-1) + H(T)[1-H(T-1)][1-A(T-1)] ) / ( [1-H(T-1)][1-A(T-1)][1-H(T)][1-A(T)] )   % arrayOfInput => (T-1):Array(1)  ,   (T):Array(2)
        % for 3 errors: ( H(T-2) + H(T-1)[1-H(T-2)][1-A(T-2)] + H(T)[1-H(T-1)][1-A(T-1)][1-H(T-2)][1-A(T-2)] )   /   ( [1-H(T-2)][1-A(T-2)][1-H(T-1)][1-A(T-1)][1-H(T)][1-A(T)] )   % arrayOfInput => (T-2):Array(1)  ,   (T-1):Array(2)   ,   (T):Array(3)
        % for n errors: Q_IOM in supplementary



    
    % (this function p_of_cr is not used in the code):
    % just for reminder:  [p_anti] = pr_anti_c_td(Input, psych_parameters, psych_parameter_strings, psych_SubjObj_flag);
    % Here, I wrote the probability of correct trial
    p_of_cr = @(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)(     (Input.RuleChoice==0)*(Input.tDev<0)*(1-pr_anti_c_td(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)) +  ...
                                                                                        (Input.RuleChoice==0)*(Input.tDev>0)*(pr_anti_c_td(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)) +  ...
                                                                                        (Input.RuleChoice==0)*(Input.tDev==0)*0.5 +  ...
                                                                                        (Input.RuleChoice==1)*(Input.tDev<0)*(pr_anti_c_td(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)) +  ...
                                                                                        (Input.RuleChoice==1)*(Input.tDev>0)*(1-pr_anti_c_td(Input, psych_model, psych_model_strings, psych_model_SubjObj_flag)) +  ...
                                                                                        (Input.RuleChoice==1)*(Input.tDev==0)*0.5);
    
    unitsq_sgm = @(x,zeta)( x.^zeta ./ (x.^zeta+(1-x).^zeta) ); 
    
    % finding the probability of switch (cosidering gaussian distribution)
    % A: sigma = sigma * T
% % % % % %     pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(T) ) ) )  ))) ) ;% adpated by wenshan 20220511 
    % ******#1 by xu 20231006
    
    % B: sigma = sigma * T^2
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( (T) ) ) )  )))    * (sigma_switch<=4)  ) ;
    % C: sigma = sigma * (scale*devlevel+T)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt( scalet * devlevel(find(dev==thrs)) + T ) ) ) )  )))  *  (scalet>0)  ) ;%* (sigma_switch>1) *  (scalet<=5) * (sigma_switch<=4) 
    % D: sigma = sigma * (scale*devlevel+T)^2
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(     ( scalet * devlevel(find(dev==thrs)) + T ) ) ) )  ))) *  (scalet>0) *  (scalet<=5) * (sigma_switch<=4) ) ;% adpated by wenshan 20220421  *  (scalet>0) *  (scalet<1)* (sigma_switch>0)
    % E: sigma = sigma * (scale*devlevel*T)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt( scalet * devlevel(find(dev==thrs)) * T ) ) ) )  )))  ) ;
    % F: sigma = sigma * (scale*devlevel*T)^2
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(     ( scalet * devlevel(find(dev==thrs)) * T ) ) ) )  ))) ) ;% adpated by wenshan 20220421  *  (scalet>0) *  (scalet<1)* (sigma_switch>0)
    % G: sigma = sigma * (exp(devlevel)*(T))
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt( scalet * exp(devlevel(find(dev==thrs))) * (T) ) ) ) )  )))  * (sigma_switch<=4) ) ;
    % F: sigma = sigma * ((exp/T)^devlevel)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt( scalet * (exp(1)/T)^(devlevel(find(dev==thrs))) ) ) ) )  )))  * (sigma_switch<=4) ) ;
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  (scalet/T)^(devlevel(find(dev==thrs))) ) ) ) )  )))   *  (scalet>2)) ;%* (sigma_switch<=4)
    % I: sigma = sigma * (T/devlevel)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  scalet* log10(T / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  scalet*T + (T / devlevel(find(dev==thrs)) ) ) ) ) )  )))  ) ;%* (sigma_switch<=4) *  abs(scalet)<=1 
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  scalet*((-1) * sqrt(T) * log(devlevel(find(dev==thrs))) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt(  ((-1) * (T^scalet) * log(devlevel(find(dev==thrs))) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  ( sqrt(T) * (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  (  scalet* ( sqrt(T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  (log(scalet* T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))  *  (scalet>0) ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  (  scalet*T + (1 / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt(  ( scalet*log(T) + (- log(devlevel(find(dev==thrs))) ) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
      
    pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt( scalet*  (T / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)

% % % % % %         pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt( scalet*  (T / (2*devlevel(find(dev==thrs))) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
    % ******#2 by xu 20231006
    
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*(  sqrt( scalet*  (T / (2*dev) ) ) ) ) )  )))   ) ;
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (1/( 1 + exp(-1* log(devlevel(find(dev==thrs))/.01) *x) )-mu_switch)/(sqrt(2)*(sigma_switch*(  ( scalet*  (T / devlevel(find(dev==thrs)) ) ) ) ) )  )))   ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( sqrt( log( scalet*T /(devlevel(find(dev==thrs))) ) )) ) ) )  ))   ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( ( log10( scalet*T / dev ) )) ) ) )  ))   ) ;%* (sigma_switch<=4)
    % pr_of_switch = @(x, mu_switch, sigma_switch, scalet, T ,dev) (  (0.5 * (1-erf( (x-mu_switch)/(sqrt(2)*(sigma_switch*( ( ( beta1*T +beta2* dev ) )) ) ) )  ))   ) ;
    

    %pr_of_switch = @(x, mu_switch, sigma_switch, T ) (  (0.5 * (1-erf((x-mu_switch)/(sqrt(2)*(sigma_switch*sqrt(T)) )))) );
    %pr_of_switch = @(x, mu_switch, sigma_switch, T ) (  (0.5 * (1-erf((x-mu_switch)/(sqrt(2)*(sigma_switch ) )))) );
    
    % cost function for finding the mean of distribution
    loss_function_pr_of_switch = @(switch_bound, mu_switch, sigma_switch, posteriorRation,dev, T) ( (pr_of_switch(switch_bound, mu_switch, sigma_switch, scalet, T, dev) - (posteriorRation/(posteriorRation+1)) )^2 );
    
    
  for iTrial = 1: length(Input)
    T = Input(iTrial).T;
    for iT = 1: T % saving the history for each trial in a new variable
        syntheticInput(iT).RuleChoice = Input(iTrial).RuleChoice(iT);
        syntheticInput(iT).Cue  = Input(iTrial).Cue(iT);
        syntheticInput(iT).tDev = Input(iTrial).tDev(iT);
        syntheticInput(iT).Anti = Input(iTrial).PrAn(iT);
        % syntheticInput(iT).Anti = (-1*Input(iTrial).PrAn(iT) +1)/2; % PrAn(-1: antisaccade, +1: prosaccade) => converted to ( 0: prosaccade, and 1: antisaccade )
        syntheticInput(iT).TF = Input(iTrial).TF(iT);
        % syntheticInput(iT).tdMean = Input(iTrial).tdMean;
        syntheticInput(iT).DevValues =Input(iTrial).DevValues;
        
        iCT_index = find( (expectedAccuracy_Benchmark.RuleChoice == Input(iTrial).RuleChoice(iT) ) .* (expectedAccuracy_Benchmark.tDev == Input(iTrial).tDev(iT) )  );
        A(iT) = expectedAccuracy_Benchmark.expectedAccuracy(iCT_index);
        H(iT) = alpha_transition;
        
    end
    
    % StartPointInitializedValues = [0.5*T];
    dev = Input(iTrial).tDev(iT);
    % dev = 0.5-abs(pr_anti_c_td(syntheticInput(end), psych_parameters, psych_parameter_strings, psych_SubjObj_flag)-.5);
    %switch_bound = 1;
    %option #1:
    % [mu_switch_estimated(iTrial), FVAL(iTrial)] = fminsearch(@(mu)  loss_function_pr_of_switch(switch_bound, mu, sigma_switch, PO(A, H, T),Input(iTrial).tDev(iT), T), StartPointInitializedValues, options);
    %option #2:
    % A: sigma = sigma * T
% % % % % %     cor = sqrt(T); % ******#1 by xu 20231006
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
    % G: sigma = sigma * (exp(devlevel)*(T))
    % cor = sqrt( scalet * exp(devlevel(find(dev==thrs))) * (T) );
    % F: sigma = sigma * ((exp/T)^devlevel)
    % cor = sqrt( scalet * (exp(1)/T)^(devlevel(find(dev==thrs))) );
    % cor = sqrt(  (scalet/T)^(devlevel(find(dev==thrs))) );
    % I: sigma = sigma * (T/devlevel)
    % cor =  sqrt( scalet*log10(T / devlevel(find(dev==thrs)) ) );
    % cor = sqrt( scalet + (T / devlevel(find(dev==thrs)) ) );
    % cor = ( scalet*T + (1 / devlevel(find(dev==thrs)) ) );
    % cor =    sqrt( scalet*((-1) * sqrt(T) * log(devlevel(find(dev==thrs))) ) );
    % cor =  sqrt( scalet*((-1) * (T^scalet) * log(devlevel(find(dev==thrs))) ) );
    % cor =    sqrt( scalet* ( sqrt(T) * ( - log(devlevel(find(dev==thrs))) ) ) ); 
    % cor =    ( scalet* ( sqrt(T) + ( - log(devlevel(find(dev==thrs))) ) ) ); 
    % cor = sqrt( ( log(scalet* T) + ( - log(devlevel(find(dev==thrs))) ) ) ); 
    % cor = sqrt( ( scalet* log(T) + ( - log(devlevel(find(dev==thrs))) ) ) ); 
    
    cor = sqrt( scalet* (T / devlevel(find(dev==thrs)) ) );

% % % % % %       cor = sqrt( scalet* (T / (2*devlevel(find(dev==thrs))) ) );
      % ******#2 by xu 20231006
      
    % cor = sqrt( scalet* (T / (2*dev) ) );
    % cor = sqrt( log(scalet* (T /(devlevel(find(dev==thrs)))) ) );
    % cor = ( log10(scalet* (T /dev) ) );
    % cor = ( (beta1*T + beta2* dev)  );
    
    
      switch_bound_t = 1;%1/( 1 + exp(-1* log10(T*devlevel(find(dev==thrs))/.1) *switch_bound));
      mu_switch_estimated(iTrial) = abs(switch_bound_t - sqrt(2) * sigma_switch * cor * erfinv( (1 - PO(A, H, T)) ./ (1 + PO(A, H, T))  ) );

    
    % implementing "alpha * mu" in the simple model
     mu_switch_estimated(iTrial) = mu_switch_estimated(iTrial) * pam3;
    % mu_switch_estimated(iTrial) = mu_switch_estimated(iTrial) * unitsq_sgm( pam3, log10(10*dev));
    % mu_switch_estimated(iTrial) = mu_switch_estimated(iTrial) *  pam3*exp(-dev);
    % mu_switch_estimated(iTrial) = mu_switch_estimated(iTrial) *  pam3 * cor;

    % to have a bound on the value of parameters:
%     if sigma_switch < 6 %rewrited from 4-->6 ,by wenshan20220219
    Output_pr_of_switch(iTrial) = pr_of_switch(switch_bound_t, mu_switch_estimated(iTrial), sigma_switch, scalet, T, dev);
    % Output_pr_of_switch(iTrial) = unitsq_sgm(Output_pr_of_switch(iTrial),log(1/dev));
    % Output_pr_of_switch(iTrial) = unitsq_sgm(Output_pr_of_switch(iTrial),scalet*log(T)-log(2*dev));
%     else
%     Output_pr_of_switch(iTrial) = 0.001; % This would help so algorithm tries to not exceed larger values (>=4)
%     end
    Output_tDev_lastOne(iTrial) = Input(iTrial).tDev(T);
    Output_RuleChoice_lastOne(iTrial) = Input(iTrial).RuleChoice(T);
    Output_T(iTrial) = T;
%     if sigma_switch < 6 %rewrited from 4-->6 ,by wenshan20220219
        Output_SW(iTrial) = Input(iTrial).SW;
%     else
%         Output_SW(iTrial) = 0; % This would make log-likelihood so small so that it helps the algorithm tries to not exceed larger values (<=4)
%     end

    clear A; clear H; clear syntheticInput;
  end

 
end

