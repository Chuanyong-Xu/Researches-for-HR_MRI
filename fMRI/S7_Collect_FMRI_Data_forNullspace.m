clear; close; clc
%load
load('/mnt/HR_project_SZU/preprocess_by_spm12/model_res_switch_fmri_m_4_29sub.mat');
%% set dir
Niidatapath = '/mnt/HR_project_SZU/preprocess_by_spm12/Decoding_mvpa/data/RawData';
matdir=dir(Niidatapath); %read the list of subjects
subdir={matdir([matdir.isdir]).name};
subdir=subdir(~ismember(subdir,{'.','..'}));
%
beha_data=Input;
for subs = 1:length(beha_data)
    infer_data.ID = beha_data(subs).ID
    infer_data.StateExp = beha_data(subs).StateExp(beha_data(subs).StateExp==1);
    infer_data.Cue = beha_data(subs).Cue(beha_data(subs).StateExp==1);
    infer_data.PrAn = beha_data(subs).PrAn(beha_data(subs).StateExp==1);
    infer_data.tDev = beha_data(subs).tDev(beha_data(subs).StateExp==1);
        infer_data.tDev(infer_data.tDev==0.1)=3;
        infer_data.tDev(infer_data.tDev==0.3)=2;
        infer_data.tDev(infer_data.tDev==0.5)=1;
        infer_data.tDev(infer_data.tDev==0.7)=2;
        infer_data.tDev(infer_data.tDev==0.9)=3;
    infer_data.RuleChoice = beha_data(subs).RuleChoice(beha_data(subs).StateExp==1);
    infer_data.TF = beha_data(subs).TF(beha_data(subs).StateExp==1);
    infer_data.SW = beha_data(subs).SW(beha_data(subs).StateExp==1);
    infer_data.Nback = beha_data(subs).Nback_act(beha_data(subs).StateExp==1);
    infer_data.mu_switch_estimated = beha_data(subs).mu_switch_estimated(beha_data(subs).StateExp==1);
    infer_data.pr_of_switch = beha_data(subs).pr_of_switch(beha_data(subs).StateExp==1);

    
    %% ROI_name = arrayfun(@(i) sprintf('Schaefer%03d',i),1:400,'UniformOutput',false);
    ROI_list = dir('mask/mask_mvpa_ner_dif.nii'); %mu_p_FWEc227_mask mask_mvpa_ner_dif
    mask = spm_read_vols(spm_vol([ROI_list(1).folder filesep ROI_list(1).name]));
    mask(mask>=1)=1;
    
    raw_vols = spm_read_vols(spm_vol([Niidatapath filesep subdir{subs} filesep 'whole_brain_average.nii.gz']));%**********
    for conds = 1:size(raw_vols,4)
        raw_vols_conds = raw_vols(:,:,:,conds);
        mask_inds = find(mask>=1);
        mask_data(:,conds) = raw_vols_conds(mask_inds);
    end
    
    infer_data.BOLD = mask_data';
    infer_data_all{subs} = infer_data;
    clear mask_inds mask_data raw_vols infer_data

end

save infer_data_all.mat infer_data_all

