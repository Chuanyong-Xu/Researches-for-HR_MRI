%% psychometric function Pr(Anti | rule, sample interval)
function [p_anti] = pr_anti_c_td(Input, parameters, parameter_strings, SubjObj_flag);
% Gaussian pdf
pdf_x = @(x, mu, sigma) ((1/sqrt(2*pi*(sigma^2))) .* exp(-((x - mu).^2) / (2*(sigma^2))) );
Pr_I_Given_C_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution) ((C==0)*sum(x_resolution*pdf_x(pBound_PrAn:x_resolution:x_max, pDrift*td, pWm*abs(td-.5)+pWm_offset)) + (C==1)*sum(x_resolution*pdf_x(x_min:x_resolution:pBound_AnPr, pDrift*td, pWm*abs(td-.5)+pWm_offset)) );
%Pr_I_Given_NC_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution) ((C==0)*sum(x_resolution*pdf_x(pBound_PrAn:x_resolution:x_max, pDrift*td, pWm*abs(td-.5)+pWm_offset)) + (C==1)*sum(x_resolution*pdf_x(x_min:x_resolution:pBound_AnPr, pDrift*td, pWm*abs(td-.5)+pWm_offset)) );
%Pr_Anti_Given_C_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, pAlphaErr1, lapseRate, x_min, x_max, x_resolution) ( ((1-(pAlphaErr1)) * ((1-lapseRate)*Pr_I_Given_C_td(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution)   + 0.5*lapseRate))  * (lapseRate<=1 && lapseRate>=0) * (pWm>0) * (pWm_offset>0) * ( abs(pAlphaErr1)<1 ) );%_pAlphErr1_test2   && pAlphaErr1>=0
%Pr_Anti_Given_C_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, pAlphaErr1, pAlphaErr2, lapseRate, x_min, x_max, x_resolution) ( (pAlphaErr2 + (1-abs(pAlphaErr1)) * ((1-lapseRate)*Pr_I_Given_C_td(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution)   + 0.5*lapseRate))  * (lapseRate<=1 && lapseRate>=0) * (pWm>0)  * (pAlphaErr1<1 && pAlphaErr1>=0) * (pAlphaErr2<1 && pAlphaErr2>=0));%_pAlphErr1_2 % * (pWm_offset>0)
%Pr_Anti_Given_C_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, pAlphaErr1, pAlphaErr2, lapseRate, x_min, x_max, x_resolution) ( (pAlphaErr2 * (((1-lapseRate)*Pr_I_Given_NC_td(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution)   + 0.5*lapseRate)) + (1-(pAlphaErr1)) * ((1-lapseRate)*Pr_I_Given_C_td(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution)   + 0.5*lapseRate))  * (lapseRate<=1 && lapseRate>=0) * (pWm>0)  * (pAlphaErr1<1 && pAlphaErr1>=0) * (pAlphaErr2<1 && pAlphaErr2>=0));%_pAlphErr1_2_test3 
Pr_Anti_Given_C_td = @(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, pAlphaErr1, pAlphaErr2, lapseRate, x_min, x_max, x_resolution) ( ((1-(pAlphaErr1)) * ((1-lapseRate)*Pr_I_Given_C_td(C, td, pDrift, pWm, pWm_offset, pBound_PrAn, pBound_AnPr, x_min, x_max, x_resolution)   + 0.5*lapseRate) +0.5*pAlphaErr1)  * (lapseRate<=1 && lapseRate>=0) * (pWm>0) * (pWm_offset>0) * ((pAlphaErr1)>0 && (pAlphaErr1)<1) );%

    for iParameters = 1: length(parameter_strings)
       eval(parameter_strings{iParameters}); 
    end
    local_parameters.pDrift = 1; % fixed
    local_parameters.pAlphaErr2 = 0;% added by wenshan 20220412

    for iTrial = 1: length(Input.tDev)
        if strcmp(SubjObj_flag, 'Subj')
        p_anti(iTrial) = Pr_Anti_Given_C_td(Input.RuleChoice(iTrial), (Input.tDev(iTrial)), local_parameters.pDrift, local_parameters.pWm, local_parameters.pWm_offset, local_parameters.pBound_PrAn, local_parameters.pBound_AnPr, local_parameters.pAlphaErr1, local_parameters.pAlphaErr2, local_parameters.lapseRate, -5, 5, .01);%%%%%%??????why the [-5 5]is selected??????
        % scalable model (uncommment for scalable model)
        % p_anti(iTrial) = Pr_Anti_Given_C_td(Input.RuleChoice(iTrial), Input.tDev(iTrial)*.001, local_parameters.pDrift, local_parameters.pWm, 0, local_parameters.pBound_PrAn, local_parameters.pBound_AnPr, local_parameters.pAlphaErr1, local_parameters.pAlphaErr2, local_parameters.lapseRate, -5, 5, 0.01);
        elseif strcmp(SubjObj_flag, 'Obj')
        p_anti(iTrial) = Pr_Anti_Given_C_td(Input.Cue(iTrial), (Input.tDev(iTrial)), local_parameters.pDrift, local_parameters.pWm, local_parameters.pWm_offset, local_parameters.pBound_PrAn, local_parameters.pBound_AnPr, local_parameters.pAlphaErr1, local_parameters.pAlphaErr2,  local_parameters.lapseRate,  -5, 5, .01);%%%%%%??????why the [-5 5]is selected??????
        % scalable model (uncommment for scalable model)
        % p_anti(iTrial) = Pr_Anti_Given_C_td(Input.Cue(iTrial), Input.tDev(iTrial), local_parameters.pDrift, local_parameters.pWm, 0, local_parameters.pBound_PrAn, local_parameters.pBound_AnPr, local_parameters.pAlphaErr1, local_parameters.pAlphaErr2, local_parameters.lapseRate, -5, 5, 0.01);
        end
    end
end
