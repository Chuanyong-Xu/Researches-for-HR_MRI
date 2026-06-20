%-----------------------------------------------------------------------
% Job saved on 14-Nov-2023 19:10:07 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
clear; close all; clc
addpath('/sharehome/xuchuanyong/Documents/Software/spm12');

%% ******
 matdir=dir('/sharehome/xuchuanyong/Documents/DATA_analysis/MRI_SZU_HR/preprocess_by_spm12/fun_s1_4d04'); %read the list of subjects
 subdir={matdir([matdir.isdir]).name};
 subdir=subdir(~ismember(subdir,{'.' '..'}));
 
 imdir2='/sharehome/xuchuanyong/Documents/DATA_analysis/MRI_SZU_HR/preprocess_by_spm12';

%% for first level model defining
mkdir GLM2_spm_1st_infer_nerDifconf/;

for ii=1:length(subdir)% number of subjects ******
% for spm first level model specification: infer
% basic parameters

 n=ii %% ******

for kk=4:6 %% ******infer
    nii_dir = dir([imdir2,filesep, 'fun_s1_4d0', num2str(kk), filesep,  subdir{ii}, filesep, 'swrafun*.nii']); % find the nii
    for iii=1:length(spm_vol([nii_dir.folder filesep nii_dir.name]))% ******
        nii_dir_sub2{iii,:}=([nii_dir(1).folder filesep nii_dir(1).name ',' num2str(iii)]);
    end
matlabbatch{n}.spm.stats.fmri_spec.dir = {[imdir2 filesep 'GLM2_spm_1st_infer_nerDifconf' filesep subdir{ii}]}; % ******
matlabbatch{n}.spm.stats.fmri_spec.timing.units = 'secs';
matlabbatch{n}.spm.stats.fmri_spec.timing.RT = 0.85;
matlabbatch{n}.spm.stats.fmri_spec.timing.fmri_t = 66;
matlabbatch{n}.spm.stats.fmri_spec.timing.fmri_t0 = 33; % ******
% infer session parameters
    load([imdir2 filesep 'multiconditions_glm_cue_feedback0' num2str(kk) filesep subdir{ii} 'nerDifconf.mat']);

matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).scans = (nii_dir_sub2); % infer******    
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(1).name = names{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(1).onset = onsets{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(1).duration = durations{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(1).tmod = 0;
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(1).pmod(1).name = pmod(1).name{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(1).pmod(1).param = pmod(1).param{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(1).pmod(1).poly = pmod(1).poly{1};

matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(2).name = names{2};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(2).onset = onsets{2};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(2).duration = durations{2};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(2).tmod = 0;
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(2).pmod(1).name = pmod(2).name{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(2).pmod(1).param = pmod(2).param{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(2).pmod(1).poly = pmod(2).poly{1};

matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).name = names{3};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).onset = onsets{3};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).duration = durations{3};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).tmod = 0;
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).pmod(1).name = pmod(3).name{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).pmod(1).param = pmod(3).param{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).pmod(1).poly = pmod(3).poly{1};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).pmod(2).name = pmod(3).name{2};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).pmod(2).param = pmod(3).param{2};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).pmod(2).poly = pmod(3).poly{2};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).pmod(3).name = pmod(3).name{3};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).pmod(3).param = pmod(3).param{3};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).pmod(3).poly = pmod(3).poly{3};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(3).orth = 0;

matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(4).name = names{4};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(4).onset = onsets{4};
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond(4).duration = durations{4};

% matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {}, 'orth', {});
% matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).multi = {[imdir2 filesep 'multiconditions_glm_cue_feedback0' num2str(kk) filesep subdir{ii} 'nerDifconf.mat']};% instr******
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).regress = struct('name', {}, 'val', {});
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).multi_reg = {[imdir2 filesep 'rp0' num2str(kk) filesep subdir{ii} '.txt']};% instr******
matlabbatch{n}.spm.stats.fmri_spec.sess(kk-3).hpf = 128; % 60
% basic parameters
matlabbatch{n}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{n}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];%*****************
matlabbatch{n}.spm.stats.fmri_spec.volt = 1;
matlabbatch{n}.spm.stats.fmri_spec.global = 'None';
matlabbatch{n}.spm.stats.fmri_spec.mthresh = 0.8;
matlabbatch{n}.spm.stats.fmri_spec.mask = {'/sharehome/xuchuanyong/Documents/Software/DPABI_V7.0_230110/Templates/BrainMask_05_91x109x91.img,1'};
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

dir_spm=dir([imdir2, filesep, 'GLM2_spm_1st_infer_nerDifconf']);
folder_names={dir_spm([dir_spm.isdir]).name};
folder_names=folder_names(~ismember(folder_names,{'.','..'}));

for iii=1:length(subdir)  % number of subjects ******
    n=iii %% ******

matlabbatch{n}.spm.stats.fmri_est.spmmat = {[imdir2 filesep 'GLM2_spm_1st_infer_nerDifconf' filesep folder_names{1,iii} filesep 'SPM.mat']};
matlabbatch{n}.spm.stats.fmri_est.write_residuals = 0;
matlabbatch{n}.spm.stats.fmri_est.method.Classical = 1;
end

spm_get_defaults;
% global defaults;
spm_jobman('initcfg');

parfor iii=1:length(matlabbatch)
    spm_jobman('run', matlabbatch(iii));%% ******
end
clear matlabbatch phase_nii mag_nii epi_nii%% ******




%% for first level contrast

dir_spm=dir([imdir2, filesep, 'GLM2_spm_1st_infer_nerDifconf']);
folder_names={dir_spm([dir_spm.isdir]).name};
folder_names=folder_names(~ismember(folder_names,{'.','..'}));

%#(muSwitch, isSwitch) to modulate feedback
for iii=1:length(subdir) %[1:6 8:12 14:18 21:28 30] 
matlabbatch{1}.spm.stats.con.spmmat = {[imdir2 filesep 'GLM2_spm_1st_infer_nerDifconf' filesep folder_names{1,iii} filesep 'SPM.mat']};
matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'cue';
matlabbatch{1}.spm.stats.con.consess{1}.tcon.weights = [1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Ruleresponse';
matlabbatch{1}.spm.stats.con.consess{2}.tcon.weights = [0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'pcpt';
matlabbatch{1}.spm.stats.con.consess{3}.tcon.weights = [0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.name = 'Pcptresponse';
matlabbatch{1}.spm.stats.con.consess{4}.tcon.weights = [0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{4}.tcon.sessrep = 'none';

matlabbatch{1}.spm.stats.con.consess{5}.tcon.name = 'Feedback';
matlabbatch{1}.spm.stats.con.consess{5}.tcon.weights = [0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{5}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.name = 'dif';
matlabbatch{1}.spm.stats.con.consess{6}.tcon.weights = [0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{6}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.name = 'ner';
matlabbatch{1}.spm.stats.con.consess{7}.tcon.weights = [0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{7}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.consess{8}.tcon.name = 'conf';
matlabbatch{1}.spm.stats.con.consess{8}.tcon.weights = [0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0];
matlabbatch{1}.spm.stats.con.consess{8}.tcon.sessrep = 'none';
matlabbatch{1}.spm.stats.con.delete = 1;
spm_jobman('run', matlabbatch); % ******
clear matlabbatch % ******
end

