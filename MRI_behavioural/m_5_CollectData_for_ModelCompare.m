clear; close all; clc
%load
load('model_res_switch_fmri_m_4_29sub.mat');
load('data_29sub_raw.mat');
delMiss       = data;
delMiss(unique([find(isnan(data(:,9)));find(isnan(data(:,11)))]),:) = [];
mkdir(pwd, '/Beha_ModCom/')
%
beha_data=Input;
for subs = 1:length(beha_data)
    ID = beha_data(subs).ID;
    StateExp = beha_data(subs).StateExp(beha_data(subs).StateExp==1);
    Rule = beha_data(subs).Cue(beha_data(subs).StateExp==1);
    % Stimulus = beha_data(subs).PrAn(beha_data(subs).StateExp==1); % is anti
    Stimulus = 1-delMiss((delMiss(:,3)==1 & delMiss(:,1)==Input(subs).ID),15); % is anti
    Response = beha_data(subs).PrAn(beha_data(subs).StateExp==1); % is key "up"/"2"/"middle finger"/"right"
    Response(Rule==1) = 1-Response(Rule==1); %% keep the original key press, different pcpt key in different rule

    % Response = delMiss((delMiss(:,3)==1 & delMiss(:,1)==Input(subs).ID),11);
    Reward = beha_data(subs).TF(beha_data(subs).StateExp==1);
    tDev = beha_data(subs).tDev(beha_data(subs).StateExp==1);
        tDev(tDev==0.1)=3;
        tDev(tDev==0.3)=2;
        tDev(tDev==0.5)=1;
        tDev(tDev==0.7)=2;
        tDev(tDev==0.9)=3;

    SW = beha_data(subs).SW(beha_data(subs).StateExp==1);
    Nback = beha_data(subs).Nback_act(beha_data(subs).StateExp==1);
    mu_switch_estimated = beha_data(subs).mu_switch_estimated(beha_data(subs).StateExp==1);
    pr_of_switch = beha_data(subs).pr_of_switch(beha_data(subs).StateExp==1);

T=table(StateExp,Rule,Stimulus,tDev,Response,Reward,SW,Nback,mu_switch_estimated,pr_of_switch);
writetable(T,[pwd, '/Beha_ModCom/', 'data_29sub_sub', num2str(ID), '.csv']);
end

