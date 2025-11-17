%% -----------------------------------------------------------------------
% Job saved on 27-Nov-2023 15:35:56 by cfg_util (rev $Rev: 7345 $)
% spm SPM - SPM12 (7771)
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
 clear; close all; clc

 matdir=dir('/mnt/HR_project_SZU/preprocess_by_spm122/fun_s1_4d01'); %read the list of subjects
 subdir={matdir([matdir.isdir]).name};
 subdir=subdir(~ismember(subdir,{'.','..'}));
 
 subdir = { 'sub1151'}; %*************

 imdir='/mnt/HR_project_SZU/Preprocess'; %% ******
 imdir2='/mnt/HR_project_SZU/preprocess_by_spm122'; %% ******


%% step1 parfor slice timing 
    n=0
    for runs=1:6 %% ******
    for subs=1:length(subdir) %% ******
        nii_dir=struct(); nii_dir_sub={1}; vols=[];
        nii_dir = dir([imdir2,filesep, 'fun_s1_4d0', num2str(runs), filesep, subdir{subs}, filesep, 'fun*.nii']); % find the nii
        for vols=1:length(spm_vol([nii_dir.folder filesep nii_dir.name])) %length(nii_dir)% ******
            nii_dir_sub{vols,:}=([nii_dir(1).folder filesep nii_dir(1).name ',' num2str(vols)]);
        end
    n=n+1

    matlabbatch{n}.spm.temporal.st.scans = {(nii_dir_sub)}';% *******note the form:{[]}' ******
    %%
    matlabbatch{n}.spm.temporal.st.nslices = 66;
    matlabbatch{n}.spm.temporal.st.tr = 0.85;
    matlabbatch{n}.spm.temporal.st.ta = 0; % 0, this will be calculated auto,if the slice order has been provided.
                                           % otherwise, must be inputed: ta=TR-(TR/nslices)
    matlabbatch{n}.spm.temporal.st.so = [0 0.4525 0.0750 0.5275 0.1500 0.6025 0.2250 0.6800 0.3025 0.7550 0.3775 0 0.4525 0.0750 0.5275 0.1500 0.6025 0.2250 0.6800 0.3025 0.7550 0.3775 0 0.4525 0.0750 0.5275 0.1500 0.6025 0.2250 0.6800 0.3025 0.7550 0.3775 0 0.4525 0.0750 0.5275 0.1500 0.6025 0.2250 0.6800 0.3025 0.7550 0.3775 0 0.4525 0.0750 0.5275 0.1500 0.6025 0.2250 0.6800 0.3025 0.7550 0.3775 0 0.4525 0.0750 0.5275 0.1500 0.6025 0.2250 0.6800 0.3025 0.7550 0.3775]*1000;
    matlabbatch{n}.spm.temporal.st.refslice = 0.3775*1000; % middle slice
    matlabbatch{n}.spm.temporal.st.prefix = 'a';
    end
    end

    spm_get_defaults;
    global defaults;
    spm_jobman('initcfg');

    parfor jobs=1:length(matlabbatch)
        spm_jobman('run', matlabbatch(jobs));%% ******
    end
    clear matlabbatch nii_dir nii_dir_sub jobs%% ******
% creat new files named 'af2022*.nii' in 'fun_s1_3d0*' after slice timing;


%% step2 parfor realign: Estimate & Reslice
n=0
    for runs=1:6 %% ******
    for subs=1:length(subdir) %% ******
        st_dir=struct(); st_dir_sub={1}; vols=[];
        st_dir = dir([imdir2,filesep, 'fun_s1_4d0', num2str(runs), filesep,  subdir{subs}, filesep, 'afun*.nii']); % find the nii
        for vols=1:length(spm_vol([st_dir.folder filesep st_dir.name]))% ******
            st_dir_sub{vols,:}=([st_dir(1).folder filesep st_dir(1).name ',' num2str(vols)]);
        end

    n=n+1

    matlabbatch{n}.spm.spatial.realign.estwrite.data = {(st_dir_sub)}';
    matlabbatch{n}.spm.spatial.realign.estwrite.eoptions.quality = 0.9;
    matlabbatch{n}.spm.spatial.realign.estwrite.eoptions.sep = 4;
    matlabbatch{n}.spm.spatial.realign.estwrite.eoptions.fwhm = 5;
    matlabbatch{n}.spm.spatial.realign.estwrite.eoptions.rtm = 1;
    matlabbatch{n}.spm.spatial.realign.estwrite.eoptions.interp = 2;
    matlabbatch{n}.spm.spatial.realign.estwrite.eoptions.wrap = [0 0 0];
    matlabbatch{n}.spm.spatial.realign.estwrite.eoptions.weight = '';
    matlabbatch{n}.spm.spatial.realign.estwrite.roptions.which = [2 1];
    matlabbatch{n}.spm.spatial.realign.estwrite.roptions.interp = 4;
    matlabbatch{n}.spm.spatial.realign.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{n}.spm.spatial.realign.estwrite.roptions.mask = 1;
    matlabbatch{n}.spm.spatial.realign.estwrite.roptions.prefix = 'r';

    end
    end
    spm_get_defaults;
    global defaults;
    spm_jobman('initcfg');

        parfor jobs=1:length(matlabbatch)
            spm_jobman('run', matlabbatch(jobs));%% ******
        end
        clear matlabbatch st_dir st_dir_sub jobs%% ******
% creat new files named 'af2022*_uw.mat' & 'craf2022*.nii' &
% 'meancraf2022*.nii' & 'rp_af2022*.txt'
% in 'fun_s1_3d0*' after slice timing;


%% step3 parfor: coregister: Estimate
n=0
for runs=1:6 %% ******
for subs=1:length(subdir) %% ******
    anat_nii=struct(); realigned_dir_mean=struct();
    %
    anat_nii=dir([imdir2 '/T1Img0' num2str(runs) filesep subdir{subs} filesep 'co*.nii']);
    
    realigned_dir_mean = dir([imdir2,filesep, 'fun_s1_4d0', num2str(runs), filesep,  subdir{subs}, filesep, 'meanafun*.nii']); % find the nii
    %
    vols=[]; realigned_dir_mean_sub={1};
    for vols=1:length(spm_vol([realigned_dir_mean.folder filesep realigned_dir_mean.name]))% ******
        realigned_dir_mean_sub{vols,:}=([realigned_dir_mean(1).folder filesep realigned_dir_mean(1).name ',' num2str(vols)]);
    end
    %
    vols=[]; realigned_dir=struct(); realigned_dir_sub={1};
    realigned_dir = dir([imdir2,filesep, 'fun_s1_4d0', num2str(runs), filesep,  subdir{subs}, filesep, 'rafun*.nii']); % find the nii
    for vols=1:length(spm_vol([realigned_dir.folder filesep realigned_dir.name]))% ******
        realigned_dir_sub{vols,:}=([realigned_dir(1).folder filesep realigned_dir(1).name ',' num2str(vols)]);
    end
    %
    n=n+1
    matlabbatch{n}.spm.spatial.coreg.estimate.ref = {[anat_nii.folder filesep anat_nii.name ',1']};
    matlabbatch{n}.spm.spatial.coreg.estimate.source = (realigned_dir_mean_sub);
    %
    matlabbatch{n}.spm.spatial.coreg.estimate.other = (realigned_dir_sub);
    %
    matlabbatch{n}.spm.spatial.coreg.estimate.eoptions.cost_fun = 'nmi';
    matlabbatch{n}.spm.spatial.coreg.estimate.eoptions.sep = [4 2];
    matlabbatch{n}.spm.spatial.coreg.estimate.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{n}.spm.spatial.coreg.estimate.eoptions.fwhm = [7 7];

end
end

spm_get_defaults;
global defaults;
spm_jobman('initcfg');

parfor jobs=1:length(matlabbatch)
    spm_jobman('run', matlabbatch(jobs));%% ******
end
clear matlabbatch realigned_dir realigned_dir_sub realigned_dir_mean_sub jobs%% ******
%creat 'c1..c6mean*.nii' segmented files,'mmeancraf2022*.nii'(bias corrected),
% 'biasField2022*.nii', 'meancraf2022_seg8.mat' in folders



%% step4 parfor normalize: Estimate & write
n=0
for runs=1:6 % ******
for subs=1:length(subdir) %% ******
    T1_df_nii=struct(); co_seg_dir=struct(); co_seg_dir_sub={1}; vols=[];
    T1_df_nii=dir([imdir2 '/T1Img0' num2str(runs) filesep subdir{subs} filesep 'co20*.nii']);
    
    co_seg_dir = dir([imdir2,filesep, 'fun_s1_4d0', num2str(runs), filesep,  subdir{subs}, filesep, 'rafun_4D.nii']); % find the nii
    for vols=1:length(spm_vol([co_seg_dir.folder filesep co_seg_dir.name]))% ******
        co_seg_dir_sub{vols,:}=([co_seg_dir(1).folder filesep co_seg_dir(1).name ',' num2str(vols)]);
    end
n=n+1

matlabbatch{n}.spm.spatial.normalise.estwrite.subj.vol = {[T1_df_nii.folder filesep T1_df_nii.name ',1']};
%
matlabbatch{n}.spm.spatial.normalise.estwrite.subj.resample = (co_seg_dir_sub);
%
matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.tpm = {'/opt/spm12/tpm/TPM.nii'};
matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
matlabbatch{n}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
matlabbatch{n}.spm.spatial.normalise.estwrite.woptions.bb = [-90 -126 -72
                                                             90 90 108] %[-78 -112 -70
                                                                       %78 76 85];
matlabbatch{n}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
matlabbatch{n}.spm.spatial.normalise.estwrite.woptions.interp = 4;
matlabbatch{n}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';

end
end
spm_get_defaults;
global defaults;
spm_jobman('initcfg');

parfor jobs=1:length(matlabbatch)
    spm_jobman('run', matlabbatch(jobs));%% ******
end
clear matlabbatch T1_df_nii co_seg_dir co_seg_dir_sub jobs%% ******

for runs=1:6 % ******
for subs=1:length(subdir) %% ******
delete(['preprocess_by_spm122/fun_s1_4d0' num2str(runs) filesep subdir{subs} filesep 'afun_*.nii']); %******
delete(['preprocess_by_spm122/fun_s1_4d0' num2str(runs) filesep subdir{subs} filesep 'rafun_*.nii']); %******
% % % % % % delete(['preprocess_by_spm122/fun_s1_4d0' num2str(runs) filesep subdir{subs} filesep 'fun_*.nii']); %******
end
end

%% smooth
n=0
for runs=1:6 % ******
for subs=1:length(subdir) %% ******
    norm_dir=struct(); norm_dir_sub={1};  vols=[];
    norm_dir = dir([imdir2, filesep, 'fun_s1_4d0', num2str(runs), filesep,  subdir{subs}, filesep, 'wrafun*.nii']); % find the nii
    for vols=1:length(spm_vol([norm_dir.folder filesep norm_dir.name]))% ******
        norm_dir_sub{vols,:}=([norm_dir(1).folder filesep norm_dir(1).name ',' num2str(vols)]);
    end
n=n+1

matlabbatch{n}.spm.spatial.smooth.data = (norm_dir_sub);%% ******
matlabbatch{n}.spm.spatial.smooth.fwhm = [6 6 6]; %% ******
matlabbatch{n}.spm.spatial.smooth.dtype = 0;
matlabbatch{n}.spm.spatial.smooth.im = 0;
matlabbatch{n}.spm.spatial.smooth.prefix = 's';
end
end
spm_get_defaults;
global defaults;
spm_jobman('initcfg');

parfor jobs=1:length(matlabbatch)
    spm_jobman('run', matlabbatch(jobs));%% ******
end
clear matlabbatch norm_dir norm_dir_sub jobs%% ******

