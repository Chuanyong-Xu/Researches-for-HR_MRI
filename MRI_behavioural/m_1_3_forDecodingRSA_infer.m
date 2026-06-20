clear; clc; close all
restoredefaultpath
%addpath(genpath('.\depend\Hierarchica_Reasoning_model'));
addpath(genpath('Hiearchica_Reasoning_model/depend/Hierarchica_Reasoning_model_changeNB'));
addpath(genpath('Hiearchica_Reasoning_model/depend/FMINSEARCHBND'));

%%prepare the data
matdir = uigetdir('','select data folder');%select the location of time series 
dir1 = dir([matdir,filesep,'*inferred_exp*.mat']);
exclud1=[]%[2 3 6 10 32 33 35 37 44     9 16 21 22 23 34 39 40];%[7 14 16 19];%%%%%%[20 21 25 27 33 35];
dir1(exclud1,:)=[];

dir0 = dir([matdir,filesep,'*instruct_exp*.mat']);
exclud2=[]%[2 3 6 10 32 33 35 37 44     9 16 21 22 23 34 39 40];%[7 14 16 19];%%%%%%[10 20];% excluding the subjects with low accuracy or bad quality
dir0(exclud2,:)=[];

load data_29sub_raw.mat;%***********

for or_subj=1:length(dir0) %number of subjects
mat1=load([matdir,filesep,dir1(or_subj).name]);
% copyfile([matdir,filesep,dir1(or_subj).name], ['/Users/vincentxu/Documents/Chen/fMRI-wenshan/wenshan']);
mat0=load([matdir,filesep,dir0(or_subj).name]);
%
    dura_Rule1 = [mat1.onset_refer - mat1.onset_rule];%
    dura_dev1 = [mat1.onset_fb - mat1.onset_refer];%
    dura_feedback1 = [mat1.onset_iti - mat1.onset_fb]+mat1.exp_setting.iti;%

    dura_Rule0 = [mat0.onset_refer - mat0.onset_rule];%
    dura_dev0 = [mat0.onset_fb - mat0.onset_refer];%
    dura_feedback0 = [mat0.onset_iti - mat0.onset_fb]+mat0.exp_setting.iti;%


%    data(((or_subj-1)*300+1):(or_subj*300),16)=[onsetFB
%    data(((or_subj-1)*300+1):(or_subj*300),17)=[onsetRule
%    data(((or_subj-1)*300+1):(or_subj*300),18)=[onsetRef
%    data(((or_subj-1)*300+1):(or_subj*300),19)=[dura_grating
    data(((or_subj-1)*300+1):(or_subj*300),21)=[dura_Rule0; dura_Rule1];%
    data(((or_subj-1)*300+1):(or_subj*300),22)=[dura_dev0; dura_dev1];%
    data(((or_subj-1)*300+1):(or_subj*300),20)=[dura_feedback0; dura_feedback1];%

end


%% for decoding; only exp=1 (inference)
load model_res_switch_fmri_m_4_29sub.mat;
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
feedbackIdx   = 14;
confidenceIdx = 23;
ruleRTIdx = 10;

delMiss       = data;
delMiss(unique([find(isnan(data(:,ruleRespIdx)));find(isnan(data(:,pcptRespIdx)))]),:) = [];
id = unique(delMiss(:,idIdx));


for subj=1:length(id)
sub_data=delMiss((delMiss(:,1)==id(subj)),:);%select subjects
sub_data(:,confidenceIdx)=Input(subj).mu_switch_estimated; %confidence
sub_data(:,26)=Input(subj).Nback_act; %error type
sub_data=sub_data(sub_data(:,3)==1,:); %only exp=1 (inference)

% for ruleRespIdx
% ---------------------sub_data_ruleresp---------------------
sub_data(:,24)=ceil((sub_data(:,17)+0.85)/0.85)+4; %onsetFB--feedback duration
% % % sub_data(:,25)=floor((sub_data(:,17)+sub_data(:,21))/0.85);
sub_data(:,25)=sub_data(:,24)+4;

for ii=1:length(sub_data)
    if 1<=sub_data(ii,4) && sub_data(ii,4)<=50
        time_run=0; chunks(ii,1)=1;
    elseif 51<=sub_data(ii,4) && sub_data(ii,4)<=100
        time_run=1148; chunks(ii,1)=2;
    elseif 101<=sub_data(ii,4) && sub_data(ii,4)<=150
        time_run=2296; chunks(ii,1)=3;
    end
time=([sub_data(ii,24)]:1:[sub_data(ii,25)])+time_run;
time_ind=num2str(time);
time_indices{ii,:}=strjoin(strsplit(time_ind),'+');

time_coo=num2str(time*0.85);
time_coords{ii,:}=strjoin(strsplit(time_coo),'+');
end
targets=sub_data(:,3);trials=sub_data(:,4)-1; context=sub_data(:,expIdxIdx);
rulelist=sub_data(:,rulelistIdx); rule_respons=sub_data(:,ruleRespIdx); %%%%%%%%
T=table(chunks,context,time_coords,time_indices,trials,rulelist,rule_respons);
writetable(T,['decoding_rsa_behaviors/' 'data_raw_decodingRSA_ruleresponse_' num2str(id(subj)) '.csv']);
clear chunks context time_coords time_indices trials rule_respons T

%---------------------for bias_levelIdx---------------------
sub_data(:,24)=ceil((sub_data(:,18)+0.85)/0.85)+4; %onsetFB--feedback duration
% % % sub_data(:,25)=floor((sub_data(:,18)+sub_data(:,22))/0.85);
sub_data(:,25)=sub_data(:,24)+4;

for ii=1:length(sub_data)
    if 1<=sub_data(ii,4) && sub_data(ii,4)<=50
        time_run=0; chunks(ii,1)=1;
    elseif 51<=sub_data(ii,4) && sub_data(ii,4)<=100
        time_run=1148; chunks(ii,1)=2;
    elseif 101<=sub_data(ii,4) && sub_data(ii,4)<=150
        time_run=2296; chunks(ii,1)=3;
    end
time=([sub_data(ii,24)]:1:[sub_data(ii,25)])+time_run;
time_ind=num2str(time);
time_indices{ii,:}=strjoin(strsplit(time_ind),'+');

time_coo=num2str(time*0.85);
time_coords{ii,:}=strjoin(strsplit(time_coo),'+');
end

targets=sub_data(:,3);trials=sub_data(:,4)-1;context=sub_data(:,expIdxIdx); 
Dev_lev=sub_data(:,bias_levelIdx);
Dev_lev(Dev_lev==1 | Dev_lev==5)=5; Dev_lev(Dev_lev==2 | Dev_lev==4)=4; %%%%%%%%
Dev_lev=Dev_lev - 2;
pcptResp=sub_data(:,pcptRespIdx); 
T=table(chunks,context,targets,time_coords,time_indices,trials,Dev_lev,pcptResp);
writetable(T,['decoding_rsa_behaviors/' 'data_raw_decodingRSA_Dev_lev_' num2str(id(subj)) '.csv']);
clear chunks context time_coords time_indices trials rule_respons T

%% ---------------------for confidence---------------------
% ***======================================================================
sub_data(:,24)=ceil((sub_data(:,16)+0.85)/0.85)+4; %onsetFB--feedback duration
% % % sub_data(:,25)=floor((sub_data(:,16)+sub_data(:,20))/0.85);
sub_data(:,25)=sub_data(:,24)+4;

for ii=1:length(sub_data)
    if 1<=sub_data(ii,4) && sub_data(ii,4)<=50
        time_run=0;chunks(ii,1)=1;
    elseif 51<=sub_data(ii,4) && sub_data(ii,4)<=100
        time_run=1148;chunks(ii,1)=2;
    elseif 101<=sub_data(ii,4) && sub_data(ii,4)<=150
        time_run=2296;chunks(ii,1)=3;
    end
time=([sub_data(ii,24)]:1:[sub_data(ii,25)])+time_run;
time_ind=num2str(time);
time_indices{ii,:}=strjoin(strsplit(time_ind),'+');

time_coo=num2str(time*0.85);
time_coords{ii,:}=strjoin(strsplit(time_coo),'+');
end

targets=sub_data(:,3);trials=sub_data(:,4)-1; context=sub_data(:,expIdxIdx);
confidence_type=sub_data(:,confidenceIdx); confidence_type(confidence_type>=1)=2;confidence_type(confidence_type<1 & confidence_type>0)=1; %%%%%%%%
%feedback=sub_data(:,feedbackIdx); 
type_dev_error(:,1)=sub_data(:,26); type_dev_error(:,2)=Dev_lev;type_dev_error(:,3)=sub_data(:,confidenceIdx);
type_dev_error(type_dev_error(:,1)==0 & type_dev_error(:,2)==3,4)=3; type_dev_error(type_dev_error(:,1)==0 & type_dev_error(:,2)==2,4)=2;
type_dev_error(type_dev_error(:,1)==0 & type_dev_error(:,2)==1,4)=1; type_dev_error(type_dev_error(:,1)==1 & type_dev_error(:,2)==3,4)=6;
type_dev_error(type_dev_error(:,1)==1 & type_dev_error(:,2)==2,4)=5; type_dev_error(type_dev_error(:,1)==1 & type_dev_error(:,2)==1,4)=4;
type_dev_error(type_dev_error(:,1)>=2 & type_dev_error(:,2)==3,4)=9; type_dev_error(type_dev_error(:,1)>=2 & type_dev_error(:,2)==2,4)=8;
type_dev_error(type_dev_error(:,1)>=2 & type_dev_error(:,2)==1,4)=7;

confidence=sub_data(:,confidenceIdx); confidence(confidence == Inf)=max(confidence(confidence < Inf));

% ---------------------for bias_level(3) & ruleResp(2) & pcptResp(2): RT of pcptResp 
ruleRT_next=[];
ruleRT_next = sub_data([2:end],ruleRTIdx); %extract the RT in next trial ***
ruleRT_next(end+1) = nan; % *** ***
% ---------------------for the consistence between ruleswitch and confidence
ruleSwitch=sub_data(:,ruleSwitchIdx); %%%%%%%% ruleswitch RSA - confidence RSA

T=table(chunks,context,targets,time_coords,time_indices,trials,ruleSwitch,...
    ruleRT_next,confidence_type,confidence,type_dev_error);
writetable(T,['decoding_rsa_behaviors/' 'data_raw_decodingRSA_confidence_' num2str(id(subj)) '.csv']);
clear chunks context time_coords time_indices trials rule_respons T type_dev_error

% ***======================================================================

%% ---------------------for ruleSwitchIdx---------------------
sub_data(:,24)=ceil((sub_data(:,16)+0.85)/0.85)+4; %onsetFB--feedback duration
% % % sub_data(:,25)=floor((sub_data(:,16)+sub_data(:,20))/0.85);
sub_data(:,25)=sub_data(:,24)+4;

for ii=1:length(sub_data)
    if 1<=sub_data(ii,4) && sub_data(ii,4)<=50
        time_run=0;chunks(ii,1)=1;
    elseif 51<=sub_data(ii,4) && sub_data(ii,4)<=100
        time_run=1148;chunks(ii,1)=2;
    elseif 101<=sub_data(ii,4) && sub_data(ii,4)<=150
        time_run=2296;chunks(ii,1)=3;
    end
time=([sub_data(ii,24)]:1:[sub_data(ii,25)])+time_run;
time_ind=num2str(time);
time_indices{ii,:}=strjoin(strsplit(time_ind),'+');

time_coo=num2str(time*0.85);
time_coords{ii,:}=strjoin(strsplit(time_coo),'+');
end

targets=sub_data(:,3);trials=sub_data(:,4)-1; context=sub_data(:,expIdxIdx);
%feedback=sub_data(:,feedbackIdx); 
confidence_type=sub_data(:,confidenceIdx); confidence_type(confidence_type>=1)=2;confidence_type(confidence_type<1 & confidence_type>0)=1; %%%%%%%%
ruleSwitch=sub_data(:,ruleSwitchIdx); %%%%%%%%
T=table(chunks,context,targets,time_coords,time_indices,trials,confidence_type,ruleSwitch);
writetable(T,['decoding_rsa_behaviors/' 'data_raw_decodingRSA_ruleSwitch_' num2str(id(subj)) '.csv']);
clear chunks context time_coords time_indices trials rule_respons T


%% ---------------------for bias_level(3) & ruleResp(2) & pcptResp(2): RT of pcptResp ----------------
% 
% 
% % A: 100x2矩阵, 第一列为类别(1~5), 第二列为数值
% 
% uniqueCats = unique(type_dev_error2(:,4));              % 找到所有类别
% numCats    = length(uniqueCats);
% 
% % 计算每个类别的均值
% catMeans   = arrayfun(@(c)nanmean(type_dev_error2(type_dev_error2(:,4)==c,5)), uniqueCats);
% 
% % 构造类别之间的距离矩阵 (一维数值 => |mean_i - mean_j|)
% distMatrix = abs(catMeans - catMeans');
% 
% distMatrix_all(:,:,subj) = distMatrix;
% % 绘制距离矩阵
% imagesc(distMatrix);
% % colormap('jet'); 
% colorbar; axis square;
% title('各类别之间距离矩阵');
% set(gca,'XTick',1:numCats,'XTickLabel',uniqueCats);
% set(gca,'YTick',1:numCats,'YTickLabel',uniqueCats);
% xlabel('类别'); ylabel('类别');
% 

end
