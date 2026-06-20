clear; clc; close all
restoredefaultpath
%addpath(genpath('.\depend\Hierarchica_Reasoning_model'));
addpath(genpath('Hiearchica_Reasoning_model/depend/Hierarchica_Reasoning_model_changeNB'));
addpath(genpath('Hiearchica_Reasoning_model/depend/FMINSEARCHBND'));
test_plot_enable = 1;
rng(66,'twister');
% rng(44,'twister');

%%prepare the data
matdir = uigetdir('','select data folder');%select the location of time series 
dir1 = dir([matdir,filesep,'*inferred_exp*.mat']);
%======  head motion: 2(), 10, 32 37 44
%======  no switch in 2-ers: 6
%======  sleep: 35
%======  poor task strategy: 9 16 21 22 23 34 39 40; or need permutation
%test in fMRI analysis (FSL randomise tfce)
exclud1=[]%[2 3 6 10 32 33 35 37 44     9 16 21 22 23 34 39 40];%[7 14 16 19];%%%%%%[20 21 25 27 33 35];
dir1(exclud1,:)=[];

dir0 = dir([matdir,filesep,'*instruct_exp*.mat']);
exclud2=[]%[2 3 6 10 32 33 35 37 44     9 16 21 22 23 34 39 40];%[7 14 16 19];%%%%%%[10 20];% excluding the subjects with low accuracy or bad quality

dir0(exclud2,:)=[];

for or_subj=1:length(dir0) %number of subjects
mat1=load([matdir,filesep,dir1(or_subj).name]);
% copyfile([matdir,filesep,dir1(or_subj).name], ['/Users/vincentxu/Documents/Chen/fMRI-wenshan/wenshan']);
mat0=load([matdir,filesep,dir0(or_subj).name]);

data(((or_subj-1)*300+1):(or_subj*300),1) = repelem(str2double(cell2mat(extractBetween(dir1(or_subj).name,'subj','_'))),300);
data(((or_subj-1)*300+1):(or_subj*300),2) = repelem(mat1.gender,300);
data(((or_subj-1)*300+1):(or_subj*300),3) = repelem([1,0],150);
data(((or_subj-1)*300+1):(or_subj*300),4) = [1:150 1:150];
data(((or_subj-1)*300+1):(or_subj*300),5) = [mat1.exp_setting.rulelist; mat0.exp_setting.rulelist];
data(((or_subj-1)*300+1):(or_subj*300),6) = [mat1.exp_setting.refAngle; mat0.exp_setting.refAngle];
data(((or_subj-1)*300+1):(or_subj*300),7) = [mat1.exp_setting.testAngle; mat0.exp_setting.testAngle];
data(((or_subj-1)*300+1):(or_subj*300),8) = data(((or_subj-1)*300+1):(or_subj*300),7)-data(((or_subj-1)*300+1):(or_subj*300),6);

% % % % % % data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(1))),9)=1001;%0.01
% % % % % % data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(2))),9)=1002;%1/8
% % % % % % data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(3))),9)=1003;%1/4
% % % % % % data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(4))),9)=1004;%3/8
% % % % % % data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(5))),9)=1005;%1/2
% % % % % % data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(6))),9)=1006;%5/8
% % % % % % data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(7))),9)=1007;%3/4
% % % % % % data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(8))),9)=1008;%7/8
% % % % % % data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(9))),9)=1009;%0.99

data((((or_subj-1)*300)+find(round(data(((or_subj-1)*300+1):(or_subj*300),8),4)==round(mat1.exp_setting.pct.thresh(1),4))),9)=1001;%0.01
data((((or_subj-1)*300)+find(round(data(((or_subj-1)*300+1):(or_subj*300),8),4)==round(mat1.exp_setting.pct.thresh(2),4))),9)=1002;%1/8
data((((or_subj-1)*300)+find(round(data(((or_subj-1)*300+1):(or_subj*300),8),4)==round(mat1.exp_setting.pct.thresh(3),4))),9)=1003;%1/4
data((((or_subj-1)*300)+find(round(data(((or_subj-1)*300+1):(or_subj*300),8),4)==round(mat1.exp_setting.pct.thresh(4),4))),9)=1004;%3/8
data((((or_subj-1)*300)+find(round(data(((or_subj-1)*300+1):(or_subj*300),8),4)==round(mat1.exp_setting.pct.thresh(5),4))),9)=1005;%1/2
% data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(6))),9)=1004;%5/8
% data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(7))),9)=1004;%3/4
% data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(8))),9)=1005;%7/8
% data((((or_subj-1)*300)+find(data(((or_subj-1)*300+1):(or_subj*300),8)==mat1.exp_setting.pct.thresh(9))),9)=1005;%0.99

data(((or_subj-1)*300+1):(or_subj*300),8) = mod(data(((or_subj-1)*300+1):(or_subj*300),9),10);%mod(vector, 10), return the last number of a vector

data(((or_subj-1)*300+1):(or_subj*300),9) =[mat1.ruleResp; mat0.ruleResp];
data(((or_subj-1)*300+1):(or_subj*300),10) = [mat1.ruleRT; mat0.ruleRT];
data(((or_subj-1)*300+1):(or_subj*300),11) = [mat1.pcptResp; mat0.pcptResp];
data(((or_subj-1)*300+1):(or_subj*300),12) = [mat1.pcptRT; mat0.pcptRT];
data(((or_subj-1)*300+1):(or_subj*300),14) = [double(mat1.ruleResp == mat1.exp_setting.rulelist & mat1.pcptResp==mat1.exp_setting.answ); double(mat0.ruleResp == mat0.exp_setting.rulelist & mat0.pcptResp==mat0.exp_setting.answ)];
data(((or_subj-1)*300+1):(or_subj*300),15) = [mat1.exp_setting.answ; mat0.exp_setting.answ];
%
    % onsetStart1 = repmat(mat1.onset_start',50,1);
    % onsetFB1  = [mat1.onset_fb - onsetStart1(:) - 10];
    % onsetRule1= [mat1.onset_rule - onsetStart1(:) - 10];
    % 
    % onsetStart0 = repmat(mat0.onset_start',50,1);
    % onsetFB0  = [mat0.onset_fb - onsetStart0(:) - 10];
    % onsetRule0= [mat0.onset_rule - onsetStart0(:) - 10];

    onsetStart1 = repmat(mat1.onset_start',50,1);
    onsetFB1  = [mat1.onset_fb - onsetStart1(:)];
    onsetRule1 = [mat1.onset_rule - onsetStart1(:)];
    onsetRef1 = [mat1.onset_refer - onsetStart1(:)]; 
    dura_grating1 = [mat1.onset_test - mat1.onset_rule]+0.21; 
    
    onsetStart0 = repmat(mat0.onset_start',50,1);
    onsetFB0  = [mat0.onset_fb - onsetStart0(:)];
    onsetRule0= [mat0.onset_rule - onsetStart0(:)];
    onsetRef0 = [mat0.onset_refer - onsetStart0(:)];
    dura_grating0 = [mat0.onset_test - mat0.onset_rule]+0.21; 

    data(((or_subj-1)*300+1):(or_subj*300),16)=[onsetFB1; onsetFB0];
    data(((or_subj-1)*300+1):(or_subj*300),17)=[onsetRule1; onsetRule0];
    data(((or_subj-1)*300+1):(or_subj*300),18)=[onsetRef1; onsetRef0];
    data(((or_subj-1)*300+1):(or_subj*300),19)=[dura_grating1; dura_grating0];
%
    feedback_acc(or_subj,1)=sum(double(mat0.ruleResp == mat0.exp_setting.rulelist & mat0.pcptResp==mat0.exp_setting.answ))/150;
    feedback_acc(or_subj,2)=sum(double(mat1.ruleResp == mat1.exp_setting.rulelist & mat1.pcptResp==mat1.exp_setting.answ))/150;

    thres_sub(or_subj,:)=mat0.exp_setting.pct.thresh;

    % #for m_2_plot* script #by 20231002
for td_i=1:5 % 9 %********************************
    % data(:,3),expIdxIdx; data(:,5),rulelistIdx; data(:,11),pcptRespIdx; data(:,8),bias_levelIdx; data(:,9),ruleRespIdx
    % instr-obj；
    respVector{1,1}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==0 & data(((or_subj-1)*300+1):(or_subj*300),5)==0 & data(((or_subj-1)*300+1):(or_subj*300),11)==1 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    respVector{1,2}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==0 & data(((or_subj-1)*300+1):(or_subj*300),5)==1 & data(((or_subj-1)*300+1):(or_subj*300),11)==1 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    % infer-subj；
    respVector{2,1}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==1 & data(((or_subj-1)*300+1):(or_subj*300),9)==0 & data(((or_subj-1)*300+1):(or_subj*300),11)==1 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    respVector{2,2}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==1 & data(((or_subj-1)*300+1):(or_subj*300),9)==1 & data(((or_subj-1)*300+1):(or_subj*300),11)==1 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    % infer-obj；
    respVector{3,1}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==1 & data(((or_subj-1)*300+1):(or_subj*300),5)==0 & data(((or_subj-1)*300+1):(or_subj*300),11)==1 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    respVector{3,2}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==1 & data(((or_subj-1)*300+1):(or_subj*300),5)==1 & data(((or_subj-1)*300+1):(or_subj*300),11)==1 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    
    % countvector
    % instr-obj；
    countVector{1,1}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==0 & data(((or_subj-1)*300+1):(or_subj*300),5)==0 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    countVector{1,2}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==0 & data(((or_subj-1)*300+1):(or_subj*300),5)==1 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    % infer-subj；
    countVector{2,1}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==1 & data(((or_subj-1)*300+1):(or_subj*300),9)==0 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    countVector{2,2}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==1 & data(((or_subj-1)*300+1):(or_subj*300),9)==1 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    % infer-obj；
    countVector{3,1}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==1 & data(((or_subj-1)*300+1):(or_subj*300),5)==0 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
    countVector{3,2}(or_subj,td_i)=sum(data(((or_subj-1)*300+1):(or_subj*300),3)==1 & data(((or_subj-1)*300+1):(or_subj*300),5)==1 & data(((or_subj-1)*300+1):(or_subj*300),8)==td_i);
end
    
end

% % for or_subj=1:(2*length(dir0)) % 2*(number of subjects)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% need consideration!!
% % for i_swi = 1:149
% %     if data(((or_subj-1)*150+1)+i_swi-1,9) ~= data(((or_subj-1)*150+1)+i_swi,9);
% %         data(((or_subj-1)*150+1)+i_swi,13) = 1;
% %     else
% %         data(((or_subj-1)*150+1)+i_swi,13) = 0;
% %     end
% % end
% % data(( (find(isnan(data(:,9))==1))-1),13)=0;  % exclude the effect of miss value
% % data(( (find(isnan(data(:,9))==1))),13)=0;    % exclude the effect of miss value
% % end

%%% (1) each run is independent.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for or_subj=1:(6*length(dir0)) % 6*(number of subjects)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% need consideration!!
% data(((or_subj-1)*50+1):((or_subj)*50),13)=[~(data(((or_subj-1)*50+1):((or_subj)*50-1),9)==data(((or_subj-1)*50+2):((or_subj)*50),9));0]; % if isswitch next trial
% data(( (find(isnan(data(:,9))==1))-1),13)=0;  % exclude the effect of miss value
% data(( (find(isnan(data(:,9))==1))),13)=0;    % exclude the effect of miss value
% end
%%% (2) each run is dependent.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%？？？？？need further checking
%%% (1) each run is independent.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for or_subj=1:(2*length(dir0)) % 6*(number of subjects)%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% need consideration!!
data(((or_subj-1)*150+1):((or_subj)*150),13)=[~(data(((or_subj-1)*150+1):((or_subj)*150-1),9)==data(((or_subj-1)*150+2):((or_subj)*150),9));0]; % if isswitch next trial
data(( (find(isnan(data(:,9))==1))-1),13)=0;  % exclude the effect of miss value
data(( (find(isnan(data(:,9))==1))),13)=0;    % exclude the effect of miss value
end



data=sortrows(data); % ******


% ******
iSubject=unique(data(:,1))
for subi=1:length(iSubject)
for itd=1:5
Dev_acc_infer(subi,itd) = sum((data(:,1) == iSubject(subi)) .* (data(1:end,3)==1).* (data(1:end,11)==data(1:end,15)).* (data(1:end,8)==itd))/sum((data(:,1) == iSubject(subi)) .* (data(1:end,3)==1).* (data(1:end,8)==itd));
Dev_acc_instr(subi,itd) = sum((data(:,1) == iSubject(subi)) .* (data(1:end,3)==0).* (data(1:end,11)==data(1:end,15)).* (data(1:end,8)==itd))/sum((data(:,1) == iSubject(subi)) .* (data(1:end,3)==0).* (data(1:end,8)==itd));
Dev_acc(subi,itd)=sum((data(:,1) == iSubject(subi)) .* (data(1:end,11)==data(1:end,15)).* (data(1:end,8)==itd))/sum((data(:,1) == iSubject(subi)) .* (data(1:end,8)==itd));

feedback_acc_infer(subi,itd)=sum((data(:,1) == iSubject(subi)) .* (data(1:end,3)==1).* (data(1:end,14)==1).* (data(1:end,8)==itd))/sum((data(:,1) == iSubject(subi)) .* (data(1:end,3)==1).* (data(1:end,8)==itd));
feedback_acc_instr(subi,itd)=sum((data(:,1) == iSubject(subi)) .* (data(1:end,3)==0).* (data(1:end,14)==1).* (data(1:end,8)==itd))/sum((data(:,1) == iSubject(subi)) .* (data(1:end,3)==0).* (data(1:end,8)==itd));

Dev_rt_rule_infer(subi,itd) = nanmean((data(:,1) == iSubject(subi)) .* (data(1:end,3)==1).* data(1:end,10).* (data(1:end,8)==itd));
Dev_rt__rule_instr(subi,itd) = nanmean((data(:,1) == iSubject(subi)) .* (data(1:end,3)==0).* data(1:end,10).* (data(1:end,8)==itd));

Dev_rt_pcpt_infer(subi,itd) = nanmean((data(:,1) == iSubject(subi)) .* (data(1:end,3)==1).* data(1:end,12).* (data(1:end,8)==itd));
Dev_rt__pcpt_instr(subi,itd) = nanmean((data(:,1) == iSubject(subi)) .* (data(1:end,3)==0).* data(1:end,12).* (data(1:end,8)==itd));
end
end

%--for MRI machine trouble, MRI scanning stopped before task end-----------
data((data(:,1)==1147 & data(:,4)>=98 & data(:,4)<=100 & data(:,3)==0),9)=NaN;
data((data(:,1)==1147 & data(:,4)>=98 & data(:,4)<=100 & data(:,3)==0),11)=NaN;

data((data(:,1)==1159 & data(:,4)>=97 & data(:,4)<=100 & data(:,3)==0),9)=NaN;
data((data(:,1)==1159 & data(:,4)>=97 & data(:,4)<=100 & data(:,3)==0),11)=NaN;
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
save data_29sub_raw.mat data

idIdx=data(:,1);gender=data(:,2);expIdx=data(:,3);trialNum=data(:,4);rulelist=data(:,5);
refAngle=data(:,6);testAngle=data(:,7);bias_level=data(:,8);ruleResp=data(:,9);
ruleRT=data(:,10);pcptResp=data(:,11);pcptRT=data(:,12);ruleSwitch=data(:,13);feedback=data(:,14);
answer_angle=data(:,15);onsetFB=data(:,16);onsetRule=data(:,17);onsetRef=data(:,18);dura_grating=data(:,19);

T=table(idIdx,gender,expIdx,trialNum,rulelist,...
    refAngle,testAngle,bias_level,ruleResp,...
    ruleRT,pcptResp,pcptRT,ruleSwitch,feedback,...
    answer_angle,onsetFB,onsetRule,onsetRef,dura_grating);
writetable(T,'data_29sub_raw.xlsx','Sheet','inferred1_instructed0');

%%
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
% % normalSubject = [1 3 5 6 8 11 13 21 24 27 28 30 31 33 202 203 207 210 305 310 401 402 407 408 410 411 412 413 415 417];%******
% % exclude2      = [2;20];%;%******
% % id            = setdiff(id,exclude2);% all 60 subjects %******
% % exclude1      = id(find(id>100));%setdiff(id,normalSubject);%******
% % [lia,excl_loc]= ismember(exclude1,id);%******
% % id            = setdiff(id,exclude1);%normalSubject';%%******


save psychometric_fmri_exp_behvaior_29subj respVector countVector id thrs
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



%% _psychometric
for iSubject = 1:length(id)    
% Step #1: find trials that are in (instructed rule task) and make an array of those trials
Index_of_trials_forFit = find( (InputAll(iSubject).StateExp==0) );
SelectedInput.tDev = InputAll(iSubject).tDev(Index_of_trials_forFit);
SelectedInput.Cue = InputAll(iSubject).Cue(Index_of_trials_forFit);
SelectedInput.RuleChoice = InputAll(iSubject).RuleChoice(Index_of_trials_forFit);
SelectedInput.TF = InputAll(iSubject).TF(Index_of_trials_forFit);
SelectedInput.PrAn = InputAll(iSubject).PrAn(Index_of_trials_forFit);
SelectedInput.DevValues = InputAll(iSubject).DevValues;

% Step #2: fitting the psychometric model (Subjective rule)
SubjObj_flag = 'Subj'; % this fits based on the subject rule tags (see paper for more info)
[EstimatedParameters, syntheticInput, parameter_strings] = ModelParameterOptimization(SelectedInput, SubjObj_flag);
model.psychometric.subj(iSubject).p_anti(1, :, :) = syntheticInput.p_anti; % fitted curve
model.psychometric.subj(iSubject).x_axis = syntheticInput.x_axis; % x axis of fit for psychometric function
model.psychometric.subj(iSubject).EstimatedParameters(1, :) = EstimatedParameters; % fitted parameters
model.psychometric.subj(iSubject).parameter_strings = parameter_strings; % name of fitted parameters
model.psychometric.subj(iSubject).SubjObj_flag = SubjObj_flag; % a tags shows subjective or objective fit
EstParameters.instrSubj(iSubject,:) = EstimatedParameters;

if test_plot_enable == 1 % plot the subjective psychometric function for the instructed task
    figure;
    temp(:,:) = model.psychometric.subj(iSubject).p_anti(1, :, :);
   x_axis = model.psychometric.subj(iSubject).x_axis;
   subplot(2,2,1); plot(x_axis, temp'); lH = legend('RuleA', 'RuleB'); set(lH, 'Box', 'off'); xlabel('td'); ylabel('Pr(up | rule, td)'); title('Subjective psychometric  in instructed rule task');   
%    axis([-.5,.5,0,1])
end


% Step #3: fitting the psychometric model (Objective rule) and make an array of those trials
SubjObj_flag = 'Obj'; % this fits based on the objective rule tags (see paper for more info)
[EstimatedParameters, syntheticInput, parameter_strings] = ModelParameterOptimization(SelectedInput, SubjObj_flag);
model.psychometric.obj(iSubject).p_anti(1, :, :) = syntheticInput.p_anti; % fitted curve
model.psychometric.obj(iSubject).x_axis = syntheticInput.x_axis;
model.psychometric.obj(iSubject).EstimatedParameters(1, :) = EstimatedParameters;
model.psychometric.obj(iSubject).parameter_strings = parameter_strings;
model.psychometric.obj(iSubject).SubjObj_flag = SubjObj_flag;
EstParameters.instrObj(iSubject,:) = EstimatedParameters;

if test_plot_enable == 1 % plot the objective psychometric function for the instructed rule task
   temp(:,:) = model.psychometric.obj(iSubject).p_anti(1, :, :);
   x_axis = model.psychometric.obj(iSubject).x_axis;
   subplot(2,2,2); plot(x_axis, temp'); lH = legend('RuleA', 'RuleB'); set(lH, 'Box', 'off'); xlabel('td'); ylabel('Pr(Anti | rule, td)'); title('Objective psychometric  in instructed rule task');
%    axis([-.5,.5,0,1])
end


% Step #4: find trials that are in (inferred rule task)
Index_of_trials_forFit = find( (InputAll(iSubject).StateExp==1) );
% make an array of those trials
SelectedInput.tDev = InputAll(iSubject).tDev(Index_of_trials_forFit);
SelectedInput.Cue = InputAll(iSubject).Cue(Index_of_trials_forFit);
SelectedInput.RuleChoice = InputAll(iSubject).RuleChoice(Index_of_trials_forFit);
SelectedInput.TF = InputAll(iSubject).TF(Index_of_trials_forFit);
SelectedInput.PrAn = InputAll(iSubject).PrAn(Index_of_trials_forFit);
% SelectedInput.tdMean = InputAll(iSubject).tdMean;
SelectedInput.DevValues = InputAll(iSubject).DevValues;

% Step #5: fitting the psychometric model (Subjective rule)
SubjObj_flag = 'Subj'; % this fits based on the subject rule tags (see paper for more info)
[EstimatedParameters, syntheticInput, parameter_strings, FVAL] = ModelParameterOptimization(SelectedInput, SubjObj_flag);
model.psychometric.subj(iSubject).p_anti(2, :, :) = syntheticInput.p_anti; % fitted curve
model.psychometric.subj(iSubject).x_axis = syntheticInput.x_axis; % x axis of fit for psychometric function
model.psychometric.subj(iSubject).EstimatedParameters(2, :) = EstimatedParameters; % fitted parameters
model.psychometric.subj(iSubject).parameter_strings = parameter_strings; % name of fitted parameters
model.psychometric.subj(iSubject).SubjObj_flag = SubjObj_flag; % a tags shows subjective or objective fit
model.psychometric.subj(iSubject).MLE = -FVAL; % log-likelihood

EstParameters.inferSubj(iSubject,:) = EstimatedParameters;


if test_plot_enable == 1 % plot the subjective psychometric function for the inferred task
   temp(:,:) = model.psychometric.subj(iSubject).p_anti(2, :, :);
   x_axis = model.psychometric.subj(iSubject).x_axis;
   subplot(2,2,3); plot(x_axis, temp'); lH = legend('RuleA', 'RuleB'); set(lH, 'Box', 'off'); xlabel('td'); ylabel('Pr(Anti | rule, td)'); title('Subjective psychometric  in inferred rule task');
%    axis([-.5,.5,0,1])
end

% Step #6: fitting the psychometric model (Objective rule) and make an array of those trials
SubjObj_flag = 'Obj'; % this fits based on the objective rule tags (see paper for more info)
[EstimatedParameters, syntheticInput, parameter_strings] = ModelParameterOptimization(SelectedInput, SubjObj_flag);
model.psychometric.obj(iSubject).p_anti(2, :, :) = syntheticInput.p_anti; % fitted curve
model.psychometric.obj(iSubject).x_axis = syntheticInput.x_axis;
model.psychometric.obj(iSubject).EstimatedParameters(2, :) = EstimatedParameters;
model.psychometric.obj(iSubject).parameter_strings = parameter_strings;
model.psychometric.obj(iSubject).SubjObj_flag = SubjObj_flag;
EstParameters.inferObj(iSubject,:) = EstimatedParameters;

if test_plot_enable == 1 % plot the objective psychometric function for the instructed rule task
   temp(:,:) = model.psychometric.obj(iSubject).p_anti(2, :, :);
   x_axis = model.psychometric.obj(iSubject).x_axis;
   subplot(2,2,4); plot(x_axis, temp'); lH = legend('RuleA', 'RuleB'); set(lH, 'Box', 'off'); xlabel('td'); ylabel('Pr(Anti | rule, td)'); title('Objective psychometric  in inferred rule task');
   %axis([0,1,0,1])
end
print(['Psychometric_' extractBefore(dir1(iSubject).name,"_exp_")], '-dtiff','-r300'); %******

end
sprintf('end!!!!!')   

save model_res_psychometric_fmri_60sub_pAlphErr1_test3 model id EstParameters
%note: psychometric model has two options for scalar and non-scalar model.
%and it can be selected in the psychometric and optimizer code.
% However, it is not suitable when td=0! (wenshan 20220103)





%%
psycho = load('model_res_psychometric_fmri_60sub_pAlphErr1_test3.mat');
model  = psycho.model;
% % % % model.psychometric.subj(excl_loc) = [];%******
% ii_dev_1=0;
% ii_dev_2=0;
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
    
    plotFigure = 1;
    if plotFigure
    figure; 
    xvalue = thrs;
    subplot(2,2,1); hold on;
    iT = 1; plot(xvalue, Pr_Sw{iT}(iSubject, :), 'o', 'Color', [1,0,0], 'LineWidth', 1.5); % subject , 1B error
    iT = 2; plot(xvalue, Pr_Sw{iT}(iSubject, :), 'o', 'Color', [0.5,0,0], 'LineWidth', 1.5); % subject , 2B error
    
    iT = 1; plot(xvalue, Inp_Sw{iT}(iSubject, :), '-', 'Color', [1,0,0],  'LineWidth', 1.5); % machine simulation (model) , 1B error
    iT = 2; plot(xvalue, Inp_Sw{iT}(iSubject, :), '-', 'Color', [0.5,0,0],  'LineWidth', 1.5); % machine simulation (model) , 2B error
    axis([0 1 0 1]);
    xticks(xvalue)
    ylabel('Pr(SW)'); %xlabel('td');
    %legend('O\_1B-Er', 'O\_2B-Er', 'M\_1B-Er', 'M\_2B-Er');
    title('fitting\_sw')
    
    subplot(2,2,2); hold on;
    iT = 1; plot(xvalue, Pr_Sw{iT}(iSubject, :), 'o', 'Color', [1,0,0], 'LineWidth', 1.5); % subject , 1B error
    iT = 2; plot(xvalue, Pr_Sw{iT}(iSubject, :), 'o', 'Color', [0.5,0,0], 'LineWidth', 1.5); % subject , 2B error
    
    iT = 1; plot(xvalue, Inp_Pr_Sw{iT}(iSubject, :), '-', 'Color', [1,0,0],  'LineWidth', 1.5); % machine simulation (model) , 1B error
    iT = 2; plot(xvalue, Inp_Pr_Sw{iT}(iSubject, :), '-', 'Color', [0.5,0,0],  'LineWidth', 1.5); % machine simulation (model) , 2B error
    axis([0 1 0 1]);
    xticks(xvalue)
    ylabel('Pr(SW)'); %xlabel('td');
   % legend('O\_1B-Er', 'O\_2B-Er', 'M\_1B-Er', 'M\_2B-Er');
    title('fitting\_pr\_sw')
    
    subplot(2,2,3); hold on;
    iT = 1; plot(xvalue, Pr_Sw{iT}(iSubject, :), 'o', 'Color', [1,0,0], 'LineWidth', 1.5); % subject , 1B error
    iT = 2; plot(xvalue, Pr_Sw{iT}(iSubject, :), 'o', 'Color', [0.5,0,0], 'LineWidth', 1.5); % subject , 2B error
    
    iT = 1; plot(xvalue, M_Sw{iT}(iSubject, :), '-', 'Color', [1,0,0],  'LineWidth', 1.5); % machine simulation (model) , 1B error
    iT = 2; plot(xvalue, M_Sw{iT}(iSubject, :), '-', 'Color', [0.5,0,0],  'LineWidth', 1.5); % machine simulation (model) , 2B error
    axis([0 1 0 1]);
    xticks(xvalue)
    ylabel('Pr(SW)'); %xlabel('td');
   % legend('O\_1B-Er', 'O\_2B-Er', 'M\_1B-Er', 'M\_2B-Er');
    title('m\_sw')
    
    subplot(2,2,4); hold on;
    iT = 1; plot(xvalue, Pr_Sw{iT}(iSubject, :), 'o', 'Color', [1,0,0], 'LineWidth', 1.5); % subject , 1B error
    iT = 2; plot(xvalue, Pr_Sw{iT}(iSubject, :), 'o', 'Color', [0.5,0,0], 'LineWidth', 1.5); % subject , 2B error
    try
    iT = 1; plot(xvalue, M_Pr_Sw{iT}(iSubject, :), '-', 'Color', [1,0,0],  'LineWidth', 1.5); % machine simulation (model) , 1B error
    iT = 2; plot(xvalue, M_Pr_Sw{iT}(iSubject, :), '-', 'Color', [0.5,0,0],  'LineWidth', 1.5); % machine simulation (model) , 2B error
    end
    axis([0 1 0 1]);
    xticks(xvalue)
    ylabel('Pr(SW)'); %xlabel('td');
  %  legend('O\_1B-Er', 'O\_2B-Er', 'M\_1B-Er', 'M\_2B-Er');
    title('m\_pr\_sw')
    end

% figure
% scatterhist(input_dev_2back(:,1),input_dev_2back(:,2),'Kernel','on','color','b')
% title('Dev in 1st and 2nd rrror',FontSize=20)
% xlabel('1st error',FontSize=15)
% ylabel('2nd error',FontSize=15)
% xticks([0.1 0.3 0.5 0.7 0.9])
% yticks([0.1 0.3 0.5 0.7 0.9])
% % Dev in 1st rrror
% figure
% histogram(input_dev_1back)
% xticks([0.1 0.3 0.5 0.7 0.9])
% histfit(input_dev_1back)
% title('Dev in 1st error only',FontSize=20)
% xlabel('1st error',FontSize=15)
% ylabel('Times of error',FontSize=15)


end

%% *****************************
for iSubject=1:length(id)

fstIdx = find(InputAll(iSubject).StateExp ==1);

    for iTrial = fstIdx(1):fstIdx(1)+length(fstIdx)-1
        if iTrial == fstIdx(1) %fstIdx(1)=151, the first trial
            InputAll(iSubject).Nback_act(iTrial,1) = double(~InputAll(iSubject).TF(iTrial));%the Nback of 151's trial = the reversal feedback of 151's；(~1)=0
        else %otherwise(152-300's trial) is equal to：
            InputAll(iSubject).Nback_act(iTrial,1) = (double(~InputAll(iSubject).TF(iTrial)) * (InputAll(iSubject).Nback_act(iTrial-1) + 1))*  (~InputAll(iSubject).SW(iTrial-1))+  (~InputAll(iSubject).TF(iTrial)) * InputAll(iSubject).SW(iTrial-1);%changed by wenshan 20220518
        end
    end
Input(iSubject).Nback_act = InputAll(iSubject).Nback_act;
end

fprintf('end!!!!!\n')
save model_res_switch_fmri_m_4_29sub estimate_parameters MachineSimulation res_o res_inp res_m id thrs Pr_Sw Inp_Pr_Sw M_Sw M_Pr_Sw Inp_Sw Input MLE

    % % Dev in 1st and 2nd rrror
    % figure
    % scatterhist(input_dev_2back(:,1),input_dev_2back(:,2),'Kernel','on','color','b')
    % title('Dev in 1st and 2nd rrror',FontSize=20)
    % xlabel('1st error',FontSize=15)
    % ylabel('2nd error',FontSize=15)
    % xticks([0.1 0.3 0.5 0.7 0.9])
    % yticks([0.1 0.3 0.5 0.7 0.9])
    % print('dev_1st&2nd_error','-dtiff','-r300')
    % close
    % % Dev in 1st rrror
    % figure
    % histogram(input_dev_1back)
    % xticks([0.1 0.3 0.5 0.7 0.9])
    % histfit(input_dev_1back)
    % title('Dev in 1st error only',FontSize=20)
    % xlabel('1st error',FontSize=15)
    % ylabel('Times of error',FontSize=15)
    % print('dev_1st_error','-dtiff','-r300')
    % close



    %% compute log-likelihood, AIC, BIC
    %1st level
    load model_res_psychometric_fmri_60sub_pAlphErr1_test3.mat

    %2nd level
    load model_res_switch_fmri_m_4_29sub.mat
    loglikelihood_2nd_m = mean(MLE);    %log-likelihood

    for subs=1:length(Input)
        ntri_subs(subs) = sum(Input(subs).StateExp==1);
        k_total = 3;
        n_total = ntri_subs(subs);
        loglikelihood_subs = MLE(subs);
        AIC_2nd(subs) = (2*k_total-2*loglikelihood_subs);
        BIC_2nd(subs) = (k_total*log(n_total)-2*loglikelihood_subs);
    end
    AIC_2nd_m = mean(AIC_2nd)
    BIC_2nd_m = mean(BIC_2nd)


    



