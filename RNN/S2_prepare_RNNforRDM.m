clear; clc
restoredefaultpath
load('/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/DeepRL_RNN/data/update_data_switch_weight_diff0.mat');
%collect data
for subs =1:length(infer_data)
    StateExp=(infer_data{subs}.StateExp)';

    Cue=(infer_data{subs}.Cue)';

    PrAn=(infer_data{subs}.PrAn)';

    tDev=(infer_data{subs}.tDev)';

    RuleChoice=(infer_data{subs}.RuleChoice)';

    TF=(infer_data{subs}.TF)';

    SW=(infer_data{subs}.SW)';

    Nback=(infer_data{subs}.Nback)';

    mu_switch_estimated=(infer_data{subs}.mu_switch_estimated)';

    pr_of_switch=(infer_data{subs}.pr_of_switch)';

    rnn_pre_prSwi=infer_data{subs}.rnn_pre;

    hidden_state=infer_data{subs}.hidden_state;
    
    conds_9 = zeros(length(Nback),1);
    conds_9(Nback==0 & tDev==1) = 1; conds_9(Nback==0 & tDev==2) = 2; conds_9(Nback==0 & tDev==3) = 3;
    conds_9(Nback==1 & tDev==1) = 4; conds_9(Nback==1 & tDev==2) = 5; conds_9(Nback==1 & tDev==3) = 6;
    conds_9(Nback>=2 & tDev==1) = 7; conds_9(Nback>=2 & tDev==2) = 8; conds_9(Nback>=2 & tDev==3) = 9;
    %
    T=table(StateExp,Cue,PrAn,tDev,RuleChoice,TF,SW,Nback,mu_switch_estimated,pr_of_switch,rnn_pre_prSwi,conds_9,hidden_state);
    writetable(T,['decoding_RDMforRNN/' 'data_RNN_prediction_' num2str(infer_data{subs}.ID) '.csv']);
end

