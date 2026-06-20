%% organize the data
clear; clc; close all

%for all the data in one file
beha_data=readtable('data_29sub_raw.xlsx');
load('model_res_switch_fmri_m_4_29sub.mat');
%for sub1147, missed 98-100th trials, corresponding to the 97-99th rows


% instr1: for each individuals' data in one single file
mkdir multiconditions_glm_cue_feedback01/
mkdir multiconditions_glm_cue_feedback02/
mkdir multiconditions_glm_cue_feedback03/
mkdir multiconditions_glm_cue_feedback04/
mkdir multiconditions_glm_cue_feedback05/
mkdir multiconditions_glm_cue_feedback06/

 matdir=dir('/mnt/HR_project_SZU/preprocess_by_spm12/fun_s1_4d04/'); %read the list of subjects
 subdir={matdir([matdir.isdir]).name};
 subdir=subdir(~ismember(subdir,{'.','..'}));
 
 imdir2='/mnt/HR_project_SZU/preprocess_by_spm12/'; %% ******
 
%% 
for ii=1:length(subdir)
    ii
    beha_data1=beha_data(((ii-1)*300+1):(ii*300),:);
    mu=isnan(beha_data1.ruleResp) | isnan(beha_data1.pcptResp);
    mu=double(mu);
    mu(mu==0)=2;
    mu=mu-1;
    Nback_act=mu;

    mu(mu==1)=Input(ii).mu_switch_estimated;
    beha_data1.mu=mu; %may need demean to avoid multicollinearity???

    Nback_act(Nback_act==1)=Input(ii).Nback_act;
    beha_data1.Nback_act=Nback_act;

    for i=4:6
    beha_data2=beha_data1(((i-1)*50+1):(i*50),:);
    beha_data2=rmmissing(beha_data2);%%******
    beha_data3=beha_data2;
    
% variables will be modulated----------------------------------------
    names=cell(1,4);
    names{1}='cue';
    names{2}='pcpt'
    names{3}='feedback';
    names{4}='isMissing';
    %
    onsets=cell(1,4);
    % cue
    onsets{1}=(beha_data2.onsetRule)';
    % pcpt
    onsets{2}=(beha_data2.onsetRef)';
    
    % feedback
    onsets{3}=(beha_data2.onsetFB)';

%*******************************************
    beha_data2=beha_data1(((i-1)*50+1):(i*50),:);%%******
% missing trials
    ruleres=beha_data2.ruleResp;
    ruleres(ruleres==0)=0; ruleres(ruleres==1)=0;ruleres(isnan(ruleres))=1;
    pcpres=beha_data2.pcptResp;
    pcpres(pcpres==0)=0; pcpres(pcpres==1)=0;pcpres(isnan(pcpres))=1;
    miss_t=ruleres + pcpres;
    miss_t(miss_t>=1)=1;

    onsets{4}=[(miss_t .* beha_data2.onsetRule);...
        (miss_t .* beha_data2.onsetRef); (miss_t .* beha_data2.onsetFB)]';
    
    onsets{4}(find(onsets{4}==0))=[];   
    if length(onsets{4})<1
        onsets{4}=990; % 711;
    else
    end

%durations
    durations=cell(1,4);
    durations{1}=repmat(1.2,1,length(onsets{1}));
    durations{2}=repmat(2.2,1,length(onsets{2}));
    durations{3}=repmat(1,1,length(onsets{3}));
    durations{4}=repmat(0,1,length(onsets{4}));

    
% pmod: variables to modulate----------------------------------------
    beha_data2=beha_data1(((i-1)*50+1):(i*50),:);
    beha_data2=rmmissing(beha_data2);%%******
    pmod=struct();
    pmod(1).name={'ruleresponse'};
    pmod(2).name={'pcptresponse'};
    pmod(3).name={'ner','Dif','conf'};
    
    
    pmod(1).poly={1};
    pmod(2).poly={1};
    pmod(3).poly={1,1,1};
    %rule response
    pmod(1).param{1}=zscore(beha_data2.ruleResp)';%demean
    %bias level
%     pmod(3).param{2}=(beha_data2.bias_level)';
%     pmod(3).param{2}(pmod(3).param{2}==5)=1;
%     pmod(3).param{2}(pmod(3).param{2}==4)=2;
%     pmod(3).param{2}=zscore(pmod(3).param{2});%demean
    %pcptresponse
    pmod(2).param{1}=zscore(beha_data2.pcptResp)';%demean
    
    
    %is error
    %pmod(3).param{2}=zscore((beha_data2.feedback==0)');%demean
    %is switch
    %------------------------------------------------------------------------------------
    beha_data2.bias_level(beha_data2.bias_level==1)=5;
    beha_data2.bias_level(beha_data2.bias_level==2)=4;
    beha_data2.bias_level=beha_data2.bias_level-2;
%     pmod(3).param{1}=(beha_data2.bias_level);
% 
%     pmod(3).param{2}=(beha_data2.Nback_act);

    pmod(3).param{1}=(beha_data2.mu);
    %------------------------------------------------------------------------------------
%     pmod_cor{ii,i}=corrcoef(pmod(3).param{1},pmod(3).param{2});
%********************************************************
    save(['multiconditions_glm_cue_feedback0' num2str(i) filesep subdir{ii} 'nerDifconf.mat'], 'names','onsets','durations','pmod');%notice:**********
%********************************************************
    end

end
%make the MRI scanning corresponding to the right behavioral data sequence

