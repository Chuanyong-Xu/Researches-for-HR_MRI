%-----------------------------------------------------------------------
% Job saved on 23-Nov-2023 09:47:18 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
clear all; close all; clc
addpath('/opt/spm12');

spm_get_defaults;
global defaults;
spm_jobman('initcfg');
%% ******

matdir='/mnt/HR_project_SZU/preprocess_by_spm12';
dir_spm=dir([matdir, filesep, 'GLM1_spm_1st_infer']);
folder_names={dir_spm([dir_spm.isdir]).name};
folder_names=folder_names(~ismember(folder_names,{'.','..'}));

for ii=1:length(folder_names)
    spm1st_all{(ii-1)*4+1,:}=([matdir filesep 'GLM1_spm_1st_infer' filesep folder_names{1,ii} filesep 'con_0001.nii']);
    spm1st_all{(ii-1)*4+2,:}=([matdir filesep 'GLM1_spm_1st_infer' filesep folder_names{1,ii} filesep 'con_0002.nii']);
    spm1st_all{(ii-1)*4+3,:}=([matdir filesep 'GLM1_spm_1st_instr' filesep folder_names{1,ii} filesep 'con_0001.nii']);
    spm1st_all{(ii-1)*4+4,:}=([matdir filesep 'GLM1_spm_1st_instr' filesep folder_names{1,ii} filesep 'con_0002.nii']);
end

imat=[ones(1,116)',repelem([1:29],4)',repmat([1;1;2;2],29,1), repmat([1;2],58,1)]; %4*num of sub


%% spm 2nd
matlabbatch{1}.spm.stats.factorial_design.dir = {'/mnt/HR_project_SZU/preprocess_by_spm12/GLM1_spm_2nd'};
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).name = 'subject';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).dept = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(1).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).name = 'task';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(2).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).name = 'iserr';
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).dept = 1;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).variance = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).gmsca = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fac(3).ancova = 0;
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.scans =(spm1st_all); % ******
%
matlabbatch{1}.spm.stats.factorial_design.des.fblock.fsuball.specall.imatrix = imat; % ******
%
matlabbatch{1}.spm.stats.factorial_design.des.fblock.maininters{1}.inter.fnums = [2
                                                                                  3];
matlabbatch{1}.spm.stats.factorial_design.cov = struct('c', {}, 'cname', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.multi_cov = struct('files', {}, 'iCFI', {}, 'iCC', {});
matlabbatch{1}.spm.stats.factorial_design.masking.tm.tm_none = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.im = 1;
matlabbatch{1}.spm.stats.factorial_design.masking.em = {'/opt/DPABI_V4.3_200401/Templates/BrainMask_05_91x109x91.img,1'}; 
matlabbatch{1}.spm.stats.factorial_design.globalc.g_omit = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.gmsca.gmsca_no = 1;
matlabbatch{1}.spm.stats.factorial_design.globalm.glonorm = 1;
%
matlabbatch{2}.spm.stats.fmri_est.spmmat(1) = cfg_dep('Factorial design specification: SPM.mat File', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{2}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;
%
matlabbatch{3}.spm.stats.con.spmmat(1) = cfg_dep('Model estimation: SPM.mat File', substruct('.','val', '{}',{2}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','spmmat'));
matlabbatch{3}.spm.stats.con.consess{1}.fcon.name = 'main_task';
matlabbatch{3}.spm.stats.con.consess{1}.fcon.weights = [1 1 -1 -1];
matlabbatch{3}.spm.stats.con.consess{1}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{2}.fcon.name = 'main_outcome';
matlabbatch{3}.spm.stats.con.consess{2}.fcon.weights = [-1 1 -1 1];
matlabbatch{3}.spm.stats.con.consess{2}.fcon.sessrep = 'none';
matlabbatch{3}.spm.stats.con.consess{3}.fcon.name = 'interact_task_outcome';
matlabbatch{3}.spm.stats.con.consess{3}.fcon.weights = [-1 1 1 -1];
matlabbatch{3}.spm.stats.con.consess{3}.fcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{4}.tcon.name = 'err_infer-err_instr';
matlabbatch{3}.spm.stats.con.consess{4}.tcon.weights = [0 1 0 -1];
matlabbatch{3}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{5}.tcon.name = 'rew_infer-rew_instr';
matlabbatch{3}.spm.stats.con.consess{5}.tcon.weights = [1 0 -1 0];
matlabbatch{3}.spm.stats.con.consess{5}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{6}.tcon.name = 'infer-instr';
matlabbatch{3}.spm.stats.con.consess{6}.tcon.weights = [1 1 -1 -1];
matlabbatch{3}.spm.stats.con.consess{6}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{7}.tcon.name = 'err-rew';
matlabbatch{3}.spm.stats.con.consess{7}.tcon.weights = [-1 1 -1 1];
matlabbatch{3}.spm.stats.con.consess{7}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{8}.tcon.name = 'rew-err';
matlabbatch{3}.spm.stats.con.consess{8}.tcon.weights = [1 -1 1 -1];
matlabbatch{3}.spm.stats.con.consess{8}.tcon.sessrep = 'none';

matlabbatch{3}.spm.stats.con.consess{9}.tcon.name = 'main_task2';
matlabbatch{3}.spm.stats.con.consess{9}.tcon.weights = [-1 -1 1 1];
matlabbatch{3}.spm.stats.con.consess{9}.tcon.sessrep = 'none';

% 
matlabbatch{3}.spm.stats.con.delete = 1;

spm_jobman('run', matlabbatch); % ******
clear matlabbatch % ******

