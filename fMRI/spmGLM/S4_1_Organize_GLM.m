%
clear; clc; close all

%for all the data in one file
beha_data=readtable('data_29sub_raw.xlsx');


% instr1: for each individuals' data in one single file
mkdir multiconditions01/
mkdir multiconditions02/
mkdir multiconditions03/
mkdir multiconditions04/
mkdir multiconditions05/
mkdir multiconditions06/

 matdir=dir('/mnt/HR_project_SZU/preprocess_by_spm12/fun_s1_4d01'); %read the list of subjects
 subdir={matdir([matdir.isdir]).name};
 subdir=subdir(~ismember(subdir,{'.','..'}));
 sub_index=subdir;
 imdir2='/mnt/HR_project_SZU/preprocess_by_spm12'; %% ******

 sub_index=table2array(unique(beha_data(:,1)));

for ii=1:length(sub_index)
    for i=1:6
    beha_data1=beha_data(((ii-1)*300+(i-1)*50+1):((ii-1)*300+i*50),:);
    beha_data1=rmmissing(beha_data1);%%******

    names=cell(1,3);
    names{1}='wrong';
    names{2}='right';
    names{3}='isMissing';
    %
    onsets=cell(1,3);
    % onsets{1}=beha_data1.expIdx .* beha_data1.onsetFB;
    % wrong
    feedback_wrong = zeros(length(beha_data1.feedback),1);
    feedback_wrong(find(beha_data1.feedback==0))=1;
    onsets{1}=feedback_wrong .* beha_data1.onsetFB;
    onsets{1}(find(onsets{1}==0))=[];
    % right
    onsets{2}=beha_data1.feedback .* beha_data1.onsetFB;
    onsets{2}(find(onsets{2}==0))=[];


%
    beha_data1=beha_data(((ii-1)*300+(i-1)*50+1):((ii-1)*300+i*50),:);%%******
    % missing trials
    ruleres=beha_data1.ruleResp;
    ruleres(ruleres==0)=0; ruleres(ruleres==1)=0;ruleres(isnan(ruleres))=1;
    pcpres=beha_data1.pcptResp;
    pcpres(pcpres==0)=0; pcpres(pcpres==1)=0;pcpres(isnan(pcpres))=1;
    miss_t=ruleres + pcpres;
    miss_t(miss_t>=1)=1;

    onsets{3}=miss_t .* beha_data1.onsetFB;
    onsets{3}(find(onsets{3}==0))=[];
    if length(onsets{3})<1
        onsets{3}=999; % 711 ************
    else
    end

    %
    durations=cell(1,3);
    durations{1}=repmat(1,length(onsets{1}),1);
    durations{2}=repmat(1,length(onsets{2}),1);
    durations{3}=repmat(0,length(onsets{3}),1);

%********************************************************
    
    save(['multiconditions0' num2str(i) filesep 'sub' num2str(sub_index(ii)) '.mat'], 'names','onsets','durations');%notice:**********

    
%********************************************************
    clear names onsets durations
    end

end
%make the MRI scanning corresponding to the right behavioral data sequence


%% ********************************notice*************************
% % % % % % movefile fun_s1_4d01/ fun_s1_4d01_prep/
% % % % % % movefile fun_s1_4d02/ fun_s1_4d02_prep/
% % % % % % movefile fun_s1_4d03/ fun_s1_4d03_prep/
% % % % % % movefile fun_s1_4d04/ fun_s1_4d04_prep/
% % % % % % movefile fun_s1_4d05/ fun_s1_4d05_prep/
% % % % % % movefile fun_s1_4d06/ fun_s1_4d06_prep/
% % % % % % 
% % % % % % for ii=1:length(subdir) % number of subjects ******
% % % % % %     for kk=1:6
% % % % % %     mkdir(['fun_s1_4d0' num2str(kk) filesep subdir{ii}]) % creat the subjects' folder
% % % % % %     end
% % % % % % end
% % % % % % 
% % % % % % sub_cha={'sub1104','sub1108','sub1111','sub1116',...
% % % % % %         'sub1117','sub1120','sub1123','sub1124','sub1126','sub1130',...
% % % % % %         'sub1142','sub1144','sub1147','sub1148','sub1155'} %*******************
% % % % % % 
% % % % % % parfor n=1:length(sub_cha)
% % % % % %     movefile(['fun_s1_4d01_prep/' sub_cha{n}], ['fun_s1_4d04/'])
% % % % % %     movefile(['fun_s1_4d02_prep/' sub_cha{n}], ['fun_s1_4d05/'])
% % % % % %     movefile(['fun_s1_4d03_prep/' sub_cha{n}], ['fun_s1_4d06/'])
% % % % % %     movefile(['fun_s1_4d04_prep/' sub_cha{n}], ['fun_s1_4d01/'])
% % % % % %     movefile(['fun_s1_4d05_prep/' sub_cha{n}], ['fun_s1_4d02/'])
% % % % % %     movefile(['fun_s1_4d06_prep/' sub_cha{n}], ['fun_s1_4d03/'])
% % % % % % end
% % % % % % 
% % % % % % sub_no_cha={'sub1101','sub1103', 'sub1105', 'sub1110','sub1114',...
% % % % % %             'sub1121','sub1122','sub1125','sub1128','sub1129','sub1132',...
% % % % % %             'sub1141','sub1146','sub1150','sub1153'} %*******************
% % % % % % % sub1121, ear hurt; 
% % % % % %         
% % % % % % parfor n=1:length(sub_no_cha)
% % % % % %     movefile(['fun_s1_4d01_prep/' sub_no_cha{n}], ['fun_s1_4d01/'])
% % % % % %     movefile(['fun_s1_4d02_prep/' sub_no_cha{n}], ['fun_s1_4d02/'])
% % % % % %     movefile(['fun_s1_4d03_prep/' sub_no_cha{n}], ['fun_s1_4d03/'])
% % % % % %     movefile(['fun_s1_4d04_prep/' sub_no_cha{n}], ['fun_s1_4d04/'])
% % % % % %     movefile(['fun_s1_4d05_prep/' sub_no_cha{n}], ['fun_s1_4d05/'])
% % % % % %     movefile(['fun_s1_4d06_prep/' sub_no_cha{n}], ['fun_s1_4d06/'])
% % % % % % end
% % % % % % clear sub_cha sub_no_cha n

% ********************************notice*************************
%% for head motion parameters
 mkdir rp01/
 mkdir rp02/
 mkdir rp03/
 mkdir rp04/
 mkdir rp05/
 mkdir rp06/
 mkdir rp_plot/
 mkdir rp_plot_FD/
 
for kk=1:6 %% ******
for ii=1:length(sub_index) %% ******
 epi_rp=struct();
 epi_rp=dir([imdir2 '/fun_s1_4d0' num2str(kk) filesep subdir{ii} filesep 'rp_afun*']);
 copyfile([epi_rp.folder filesep epi_rp.name], ['rp0' num2str(kk) filesep subdir{ii} '.txt']);
  
% plot head motions
%   headm = [];
%   headm=load([epi_rp.folder filesep epi_rp.name]);
%   figure
%   set(gcf,'unit','pixels','position',[200 50 600 250])
%   plot(headm)
%   print(['rp_plot/fun_s1_4d0' num2str(kk) '_' subdir{ii}], '-dtiff', '-r300') 
%   close
  
%   max_headm{kk}(ii,:)=max(abs(headm));
  
  % calculate the FD value
%   translations = headm(:,1:3); %x, y, z
%   rotations = headm(:,4:6); % pitch, roll, yaw rotation
%   
%   radius = 50; %mm
%   rotations_mm = rotations * (radius * pi /180);
%   
%   num_frames = size(headm,1);
%   FD_trans = zeros(num_frames-1, 1);
%   FD_rot = zeros(num_frames-1, 1);
%   for t = 2:num_frames
%       FD_trans(t-1) = sqrt((translations(t,1)-translations(t-1,1))^2 + ...
%           (translations(t,2)-translations(t-1,2))^2 + (translations(t,3)-translations(t-1,3))^2);
%       FD_rot(t-1) = sqrt((rotations(t,1)-rotations(t-1,1))^2 + ...
%           (rotations(t,2)-rotations(t-1,2))^2 + (rotations(t,3)-rotations(t-1,3))^2);
%   end
  
% % % % % %   for t = 2:num_frames
% % % % % %       delta_trans = abs(translations(t,:)-translations(t-1,:));
% % % % % %       delta_rot = abs(rotations_mm(t,:)-rotations_mm(t-1,:));
% % % % % %       FD(t-1) = sum(delta_trans)+sum(delta_rot);
% % % % % %       mean_FD(t-1) = mean([delta_trans delta_rot]);
% % % % % %       max_FD(t-1) = max([delta_trans delta_rot]);
% % % % % %   end
  
%   figure
%   set(gcf,'unit','pixels','position',[200 50 600 250])
%   plot(FD_trans)
%   hold on
%   plot(FD_rot)
%   print(['rp_plot_FD/fun_s1_4d0' num2str(kk) '_' subdir{ii}], '-dtiff', '-r300') 
%   close
  
  
end
end
