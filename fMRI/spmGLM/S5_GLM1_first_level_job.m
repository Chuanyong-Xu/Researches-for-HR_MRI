%-----------------------------------------------------------------------
% Job saved on 14-Nov-2023 19:10:07 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
clear; close all; clc
restoredefaultpath
addpath('/opt/spm12');

%% ******
 matdir=dir('/mnt/HR_project_SZU/preprocess_by_spm12/fun_s1_4d01'); %read the list of subjects
 subdir={matdir([matdir.isdir]).name};
 subdir=subdir(~ismember(subdir,{'.','..'}));

 imdir2='/mnt/HR_project_SZU/preprocess_by_spm12';


%% for first level model defining
mkdir GLM1_spm_1st_instr/;
mkdir GLM1_spm_1st_infer/;

for ii=1:length(subdir) % number of subjects ******
% reading the file structs and names
% instr files

 n=ii %% ******
 
for kk=1:3 %% ******instr
    nii_dir=struct(); nii_dir_sub={1};
    nii_dir = dir([imdir2,filesep, 'fun_s1_4d0', num2str(kk), filesep,  subdir{ii}, filesep, 'swrafun*.nii']); % find the nii
    for iii=1:length(spm_vol([nii_dir.folder filesep nii_dir.name]))% ******
        nii_dir_sub{iii,:}=([nii_dir(1).folder filesep nii_dir(1).name ',' num2str(iii)]);
    end

% for spm first level model specification: instr
% basic parameters
matlabbatch{n}.spm.stats.fmri_spec.dir = {[imdir2 filesep 'GLM1_spm_1st_instr' filesep subdir{ii}]}; % ******
matlabbatch{n}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{n}.spm.stats.fmri_spec.timing.RT = 0.85;
matlabbatch{n}.spm.stats.fmri_spec.timing.fmri_t = 66;
matlabbatch{n}.spm.stats.fmri_spec.timing.fmri_t0 = 33; % ******????????????
% instr session parameters
matlabbatch{n}.spm.stats.fmri_spec.sess(kk).scans = (nii_dir_sub); % instr******
matlabbatch{n}.spm.stats.fmri_spec.sess(kk).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{n}.spm.stats.fmri_spec.sess(kk).multi = {[imdir2 filesep 'multiconditions0' num2str(kk) filesep subdir{ii} '.mat']};% instr******
matlabbatch{n}.spm.stats.fmri_spec.sess(kk).regress = struct('name', {}, 'val', {});
matlabbatch{n}.spm.stats.fmri_spec.sess(kk).multi_reg = {[imdir2 filesep 'rp0' num2str(kk) filesep subdir{ii} '.txt']};% instr******
matlabbatch{n}.spm.stats.fmri_spec.sess(kk).hpf = 128; % 60; % ******
% basic parameters
matlabbatch{n}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{n}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];%*****************
% matlabbatch{n}.spm.stats.fmri_spec.bases.fir.length = 19;
% matlabbatch{n}.spm.stats.fmri_spec.bases.fir.order = 22;

matlabbatch{n}.spm.stats.fmri_spec.volt = 1;
matlabbatch{n}.spm.stats.fmri_spec.global = 'None';
matlabbatch{n}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{n}.spm.stats.fmri_spec.mask = {'/opt/DPABI_V4.3_200401/Templates/BrainMask_05_91x109x91.img,1'};
matlabbatch{n}.spm.stats.fmri_spec.cvi = 'AR(1)';
end


% for spm first level model specification: infer
% basic parameters

 n=ii+length(subdir) %% ******

for kk=4:6 %% ******infer
    nii_dir=struct(); nii_dir_sub2={1};
    nii_dir = dir([imdir2,filesep, 'fun_s1_4d0', num2str(kk), filesep,  subdir{ii}, filesep, 'swrafun*.nii']); % find the nii
    for iii=1:length(spm_vol([nii_dir.folder filesep nii_dir.name]))% ******
        nii_dir_sub2{iii,:}=([nii_dir(1).folder filesep nii_dir(1).name ',' num2str(iii)]);
    end
matlabbatch{n}.spm.stats.fmri_spec.dir = {[imdir2 filesep 'GLM1_spm_1st_infer' filesep subdir{ii}]}; % ******
matlabbatch{n}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{n}.spm.stats.fmri_spec.timing.RT = 0.85;
matlabbatch{n}.spm.stats.fmri_spec.timing.fmri_t = 66;
matlabbatch{n}.spm.stats.fmri_spec.timing.fmri_t0 = 33; % ******
% infer session parameters
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).scans = (nii_dir_sub2); % infer******
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).multi = {[imdir2 filesep 'multiconditions0' num2str(kk) filesep subdir{ii} '.mat']};% infer******
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).regress = struct('name', {}, 'val', {});
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).multi_reg = {[imdir2 filesep 'rp0' num2str(kk) filesep subdir{ii} '.txt']};% infer******
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).hpf = 128; % 60; % ******
% basic parameters
matlabbatch{n}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{n}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
% matlabbatch{n}.spm.stats.fmri_spec.bases.fir.length = 19;
% matlabbatch{n}.spm.stats.fmri_spec.bases.fir.order = 22;

matlabbatch{n}.spm.stats.fmri_spec.volt = 1;
matlabbatch{n}.spm.stats.fmri_spec.global = 'None';
matlabbatch{n}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{n}.spm.stats.fmri_spec.mask = {'/opt/DPABI_V4.3_200401/Templates/BrainMask_05_91x109x91.img,1'};
matlabbatch{n}.spm.stats.fmri_spec.cvi = 'AR(1)';

end
end
spm_get_defaults;
global defaults;
spm_jobman('initcfg');

parfor iii=1:length(matlabbatch)
    spm_jobman('run', matlabbatch(iii));%% ******
end
clear matlabbatch phase_nii mag_nii epi_nii%% ******



%% for first level model estimations

dir_spm=dir([imdir2, filesep, 'GLM1_spm_1st_instr']);
folder_names={dir_spm([dir_spm.isdir]).name};
folder_names=folder_names(~ismember(folder_names,{'.','..'}));

for iii=1:length(subdir)  % number of subjects ******
    n=iii %% ******
matlabbatch{n}.spm.stats.fmri_est.spmmat = {[imdir2 filesep 'GLM1_spm_1st_instr' filesep folder_names{1,iii} filesep 'SPM.mat']};
matlabbatch{n}.spm.stats.fmri_est.write_residuals = 0; %% ******
matlabbatch{n}.spm.stats.fmri_est.method.Classical = 1;

matlabbatch{n+length(subdir)}.spm.stats.fmri_est.spmmat = {[imdir2 filesep 'GLM1_spm_1st_infer' filesep folder_names{1,iii} filesep 'SPM.mat']};
matlabbatch{n+length(subdir)}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{n+length(subdir)}.spm.stats.fmri_est.method.Classical = 1;

end

spm_get_defaults;
global defaults;
spm_jobman('initcfg');

parfor iii=1:length(matlabbatch)
    spm_jobman('run', matlabbatch(iii));%% ******
end
clear matlabbatch phase_nii mag_nii epi_nii%% ******




%% for first level contrast
dir_spm=dir([imdir2, filesep, 'GLM1_spm_1st_infer']);
folder_names={dir_spm([dir_spm.isdir]).name};
folder_names=folder_names(~ismember(folder_names,{'.','..'}));

for iii=1:length(subdir)
matlabbatch{1}.spm.stats.con.spmmat = {[imdir2 filesep 'GLM1_spm_1st_infer' filesep folder_names{1,iii} filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'infer_rew';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [0 1/3 0 0 0 0 0 0 0 0 1/3 0 0 0 0 0 0 0 0 1/3 0 0 0 0 0 0 0 0 0 0]; %[isSwitch noSwitch miss 6rp parameters]
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none'; 
% if 'none', the weights must be divided by number of runs; if repl, should not be divided.
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'infer_err';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [1/3 0 0 0 0 0 0 0 0 1/3 0 0 0 0 0 0 0 0 1/3 0 0 0 0 0 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
%----------------
matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'infer_err-rew';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [1/3 -1/3 0 0 0 0 0 0 0 1/3 -1/3 0 0 0 0 0 0 0 1/3 -1/3 0 0 0 0 0 0 0]; %[isSwitch noSwitch miss 6rp parameters]
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
%----------------

matlabbatch{2}.spm.stats.con.spmmat = {[imdir2 filesep 'GLM1_spm_1st_instr' filesep folder_names{1,iii} filesep 'SPM.mat']};
matlabbatch{2}.spm.stats.con.consess{1}.tcon.name = 'instr_rew';
matlabbatch{2}.spm.stats.con.consess{1}.tcon.weights = [0 1/3 0 0 0 0 0 0 0 0 1/3 0 0 0 0 0 0 0 0 1/3 0 0 0 0 0 0 0 0 0 0];
matlabbatch{2}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.consess{2}.tcon.name = 'instr_err';
matlabbatch{2}.spm.stats.con.consess{2}.tcon.weights = [1/3 0 0 0 0 0 0 0 0 1/3 0 0 0 0 0 0 0 0 1/3 0 0 0 0 0 0 0 0 0 0 0];
matlabbatch{2}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
%----------------
matlabbatch{2}.spm.stats.con.consess{3}.tcon.name = 'instr_err-rew';
matlabbatch{2}.spm.stats.con.consess{3}.tcon.weights = [1/3 -1/3 0 0 0 0 0 0 0 1/3 -1/3 0 0 0 0 0 0 0 1/3 -1/3 0 0 0 0 0 0 0];
matlabbatch{2}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{2}.spm.stats.con.delete = 1;
%----------------
spm_jobman('run', matlabbatch); % ******
clear matlabbatch % ******
end

