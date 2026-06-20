%-----------------------------------------------------------------------
% Job saved on 13-Dec-2023 19:29:34 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
clear; close all; clc

matdir='/sharehome/xuchuanyong/Documents/DATA_analysis/MRI_SZU_HR/preprocess_by_spm12';
dir_spm=dir([matdir, filesep, 'GLM2_spm_1st_infer_nerDifconf']);
folder_names={dir_spm([dir_spm.isdir]).name};
folder_names=folder_names(~ismember(folder_names,{'.','..'}));
k=0
for ii=1:length(folder_names); %[1:5 7:11 13:16 20:27 29:30 33:37];
    k=k+1
     spm1st_Dif{k,:}=([matdir filesep 'GLM2_spm_1st_infer_nerDifconf' filesep folder_names{1,ii} filesep 'con_0006.nii']);
     spm1st_ner{k,:}=([matdir filesep 'GLM2_spm_1st_infer_nerDifconf' filesep folder_names{1,ii} filesep 'con_0007.nii']);
     spm1st_conf{k,:}=([matdir filesep 'GLM2_spm_1st_infer_nerDifconf' filesep folder_names{1,ii} filesep 'con_0008.nii']);     
end


%% ner
spm_get_defaults;
global defaults;
spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {'/sharehome/xuchuanyong/Documents/DATA_analysis/MRI_SZU_HR/preprocess_by_spm12/GLM2_spm_2nd_nerDifconf_ner'};% mk direc to save
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = (spm1st_ner);

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;%************************
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/sharehome/xuchuanyong/Documents/Software/DPABI_V7.0_230110/Templates/BrainMask_05_91x109x91.img,1'}% explicit mask
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'ner_P';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'ner_N';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch); % ******
clear matlabbatch % ******


%% Dif
spm_get_defaults;
global defaults;
spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {'/sharehome/xuchuanyong/Documents/DATA_analysis/MRI_SZU_HR/preprocess_by_spm12/GLM2_spm_2nd_nerDifconf_Dif'};% mk direc to save
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = (spm1st_Dif);

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;%************************
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/sharehome/xuchuanyong/Documents/Software/DPABI_V7.0_230110/Templates/BrainMask_05_91x109x91.img,1'}% explicit mask
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'Dif_P';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Dif_N';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch); % ******
clear matlabbatch % ******


%% ner
spm_get_defaults;
global defaults;
spm_jobman('initcfg');

matlabbatch{1}.spm.stats.factorial_design.dir = {'/sharehome/xuchuanyong/Documents/DATA_analysis/MRI_SZU_HR/preprocess_by_spm12/GLM2_spm_2nd_nerDifconf_conf'};% mk direc to save
matlabbatch{1}.spm.stats.factorial_design.des.t1.scans = (spm1st_conf);

matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;%************************
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/sharehome/xuchuanyong/Documents/Software/DPABI_V7.0_230110/Templates/BrainMask_05_91x109x91.img,1'}% explicit mask
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;

matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'conf_P';
matlabbatch{3}.spm.stats.con.consess{1}.tcon.weights = 1;
matlabbatch{3}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'conf_N';
matlabbatch{3}.spm.stats.con.consess{2}.tcon.weights = -1;
matlabbatch{3}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch); % ******
clear matlabbatch % ******