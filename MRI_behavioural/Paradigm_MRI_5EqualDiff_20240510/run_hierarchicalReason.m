% % Get threshold; First day
% % Subject Info
% clc;clear; close all
% %******************************
% ID         =1101;
% sub_name   = 'LDW';
% gender     = 1; % man:1; woman:
% %******************************
% isExp = 1; % isExp: 0: practice (need >90% accuracy);1 experiment
% [thresh] = threshGet(sub_name,ID,gender, isExp)


%% Reasoning Task; Second day1
% Subject Info
clc;clear; close all
%******************************
ID         = 1101;
sub_name   = 'HRC';
gender     = 1; % man:1; woman:2
%******************************
devIdx = 0; % dev practice:1; no dev practice:0;
expIdx = 1; % 0: instructed; 1: infered1
%******************************
isExp  = 1; % isExp: 0: practice;1 experiment
%******************************

[thresh_glm thresh_WH] = calThresh(ID); %if use the 'thresh_WH', remember to change the UL & LL in calThresh.m & threshGet.m
thresh=thresh_glm;
% % % thresh=thresh_WH;


pct    = perceptureSet(thresh,expIdx,ID);
pct.thresh = thresh;
if devIdx==1
hierarchicalReason_dev(sub_name,ID,gender,expIdx,pct,isExp);
else
hierarchicalReason(sub_name,ID,gender,expIdx,pct,isExp);
end
clear