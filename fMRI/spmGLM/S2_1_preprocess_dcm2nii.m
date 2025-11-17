%-----------------------------------------------------------------------
% Job saved on 24-Nov-2023 18:29:41 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------

%% ******
addpath('/opt/spm12');
clear; close all; clc

spm_get_defaults;
global defaults;
spm_jobman('initcfg');

%% 4D to 3D split
matdir=dir('/mnt/HR_project_SZU/Preprocess/S1_FunRaw'); %read the list of subjects
subdir={matdir([matdir.isdir]).name};
subdir=subdir(~ismember(subdir,{'.','..'})); %excluding the first two characters

subdir = {'sub1151'}; %*************

if ~exist('preprocess_by_spm122', 'dir')
    mkdir('preprocess_by_spm122');
end

imdir='/mnt/HR_project_SZU/Preprocess';
imdir2='/mnt/HR_project_SZU/preprocess_by_spm122'; %% ******
 
    n=0
for subs=1:length(subdir) % number of subjects ******
    for runs=1:6
        mkdir(['preprocess_by_spm122/fun_s1_4d0' num2str(runs) filesep subdir{subs}]) % creat the subjects' folder
        n=n+1
        dcm_dir=struct(); dcm_dir_sub={1};
        dcm_dir = dir([imdir,filesep, 'S', num2str(runs), '_FunRaw', filesep,  subdir{subs}, filesep, '*.IMA']); % find the DICOM
%     dcm_dir(1:2)=[];% %excluding the first two characters, *******************************
        for iii=1:length(dcm_dir)% number of dicom files******
            dcm_dir_sub{iii,:}=([dcm_dir(iii).folder filesep dcm_dir(iii).name]);
        end

        matlabbatch{n}.spm.util.import.dicom.data = (dcm_dir_sub);
        matlabbatch{n}.spm.util.import.dicom.root = 'flat';
        matlabbatch{n}.spm.util.import.dicom.outdir = {['preprocess_by_spm122/fun_s1_4d0' num2str(runs) filesep subdir{subs}]};%output dir
        matlabbatch{n}.spm.util.import.dicom.protfilter = '.*';
        matlabbatch{n}.spm.util.import.dicom.convopts.format = 'nii';
        matlabbatch{n}.spm.util.import.dicom.convopts.meta = 0;
        matlabbatch{n}.spm.util.import.dicom.convopts.icedims = 0;
    end
end
    spm_get_defaults;
    global defaults;
    spm_jobman('initcfg');

    parfor jobs=1:length(matlabbatch)
        spm_jobman('run', matlabbatch(jobs));%% ******
    end
    clear matlabbatch dcm_dir dcm_dir_sub%% ******

%%  
    n=0
for runs=1:6 %% ******
    for subs=1:length(subdir) %% ******
        nii_dir=struct(); nii_dir_sub={1};
        nii_dir = dir([imdir2,filesep, 'fun_s1_4d0', num2str(runs), filesep, subdir{subs}, filesep, 'fS*.nii']); % find the nii
        for iii=1:length(nii_dir)% number of dicom files******
            nii_dir_sub{iii,:}=([nii_dir(iii).folder filesep nii_dir(iii).name]);
        end
    n=n+1 

        matlabbatch{n}.spm.util.cat.vols = (nii_dir_sub);
% % % % % %     matlabbatch{2}.spm.util.cat.name = 'original_fun_4D.nii';
        matlabbatch{n}.spm.util.cat.name = 'fun_4D.nii';
        matlabbatch{n}.spm.util.cat.dtype = 4;
        matlabbatch{n}.spm.util.cat.RT = 0.85;
     
    end
end

    spm_get_defaults;
    global defaults;
    spm_jobman('initcfg');

    parfor jobs=1:length(matlabbatch)
        spm_jobman('run', matlabbatch(jobs));%% ******
    end
    clear matlabbatch nii_dir nii_dir_sub%% ******

%%
for runs=1:6 % ******
    for subs=1:length(subdir) %% ******
        delete(['preprocess_by_spm122/fun_s1_4d0' num2str(runs) filesep subdir{subs} filesep 'fS*.nii']);

%     delete(['preprocess_by_spm122/fun_s1_3d0' num2str(runs) filesep subdir{subs} filesep 'f20*-00001-000001-01.nii']) % delete the first 5 time points
%     delete(['preprocess_by_spm122/fun_s1_3d0' num2str(runs) filesep subdir{subs} filesep 'f20*-00002-000002-01.nii']) % delete the first 5 time points
%     delete(['preprocess_by_spm122/fun_s1_3d0' num2str(runs) filesep subdir{subs} filesep 'f20*-00003-000003-01.nii']) % delete the first 5 time points
%     delete(['preprocess_by_spm122/fun_s1_3d0' num2str(runs) filesep subdir{subs} filesep 'f20*-00004-000004-01.nii']) % delete the first 5 time points
%     delete(['preprocess_by_spm122/fun_s1_3d0' num2str(runs) filesep subdir{subs} filesep 'f20*-00005-000005-01.nii']) % delete the first 5 time points

    end
% anatomical images

end
