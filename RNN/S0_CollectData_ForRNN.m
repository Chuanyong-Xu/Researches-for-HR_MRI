%load
load('model_res_switch_fmri_m_4_29sub.mat');
%
beha_data=Input;
for subs = 1:length(beha_data)
    infer_data(subs).ID = beha_data(subs).ID
    infer_data(subs).StateExp = beha_data(subs).StateExp(beha_data(subs).StateExp==1);
    infer_data(subs).Cue = beha_data(subs).Cue(beha_data(subs).StateExp==1);
    infer_data(subs).PrAn = beha_data(subs).PrAn(beha_data(subs).StateExp==1);
    infer_data(subs).tDev = beha_data(subs).tDev(beha_data(subs).StateExp==1);
        infer_data(subs).tDev(infer_data(subs).tDev==0.1)=3;
        infer_data(subs).tDev(infer_data(subs).tDev==0.3)=2;
        infer_data(subs).tDev(infer_data(subs).tDev==0.5)=1;
        infer_data(subs).tDev(infer_data(subs).tDev==0.7)=2;
        infer_data(subs).tDev(infer_data(subs).tDev==0.9)=3;
    infer_data(subs).RuleChoice = beha_data(subs).RuleChoice(beha_data(subs).StateExp==1);
    infer_data(subs).TF = beha_data(subs).TF(beha_data(subs).StateExp==1);
    infer_data(subs).SW = beha_data(subs).SW(beha_data(subs).StateExp==1);
    infer_data(subs).Nback = beha_data(subs).Nback_act(beha_data(subs).StateExp==1);
    infer_data(subs).mu_switch_estimated = beha_data(subs).mu_switch_estimated(beha_data(subs).StateExp==1);
    infer_data(subs).pr_of_switch = beha_data(subs).pr_of_switch(beha_data(subs).StateExp==1);
end

save infer_data.mat infer_data