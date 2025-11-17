function [thresh thresh_WH] = calThresh(subjectID)
% function [thresh] = calThresh(subjectID)
% thrs = [.01, 1/16, 1/8, 1/4, 1/2, 3/4, 7/8, 15/16, .99];
% thrs = [.01, 1/8, 1/4, 3/8, 1/2, 5/8, 3/4, 7/8, .99]; % changed by wenshan 20211228  %****** changed by xu 20240305
% % thrs = [1/6, 2/6, 3/6, 4/6, 5/6];  %****** changed by xu 20240305
    thrs= [0.1, 0.3, 0.5, 0.7, 0.9];%****** changed by xu 20240305
    
file = dir(['./result/subj',num2str(subjectID),'_threshGet_exp_*']);

if length(file) ==1
    load(fullfile(file.folder,file.name))
    data = [testAng'; respVector; countVector]';
    [co, curve, thresh] = ...
        fitPsyche.fitPsycheCurveLogit(data(:, 1), data(:, 2) ./ data(:, 3),'targets',thrs);
    thresh = round(thresh,4); % threshold calculated by GLM
    % ******plot the figure

%****************************************************************************% threshold calculated by WH2001
    % UpperLimits:
    UL = [0.01, 0.01,18,18]; % Limit upper bound of g and l to 5%
%     UL = [0.01, 0.01,16.5,16.5]; % Limit upper bound of g and l to 5%
    % StartPoints:
    SP = [0, 0, 0.5, 0.5];
    % LowerLimits:
    LL = [0.01, 0.01, -18,-18];
%     LL = [0.01, 0.01, -16.5,-16.5];
    
    ffit1=fitPsyche(data(:, 1), data(:, 2) ./ data(:, 3),'WH',[UL;SP;LL]);
%     ffit1=fitPsyche(data(:, 1), data(:, 2) ./ data(:, 3),'WH');
        display(ffit1.model)
    for thrs_i=1:length(thrs)
        F = @(g,l,u,v,x) g+(1-g-l)*0.5*(1+erf((x-u)/sqrt(2*v^2))); % from fitPsyche.fitPsycheCurveWH.m 
        target_F = thrs(thrs_i);
        equation = @(x) F(ffit1.model.g,ffit1.model.l,ffit1.model.u,ffit1.model.v,x) - target_F;
        initial_guess = 0.5;
        x_solution = fzero(equation, initial_guess);
%         disp(['对应的x值为: ', num2str(x_solution)]);
        thresh_WH(thrs_i)=x_solution; % threshold calculated by WH2001
    end
%****************************************************************************

else
    print = 'No file or multi files were found!'
end
end