%% setup
clear; clc
restoredefaultpath
addpath('/opt/spm12')
%% set dir
matdir=dir('/mnt/HR_project_SZU/preprocess_by_spm12/tRSA_confidence_9cond_fan_searchlight_10mm_p/data'); %read the list of subjects
subdir={matdir([matdir.isdir]).name};
subdir=subdir(~ismember(subdir,{'.','..'}));

% ROI_name = arrayfun(@(i) sprintf('Schaefer%03d',i),1:400,'UniformOutput',false);
ROI_list = dir('mask/mask_mvpa_ner_dif.nii'); %mu_p_FWEc227_mask mask_mvpa_ner_dif
mask = spm_read_vols(spm_vol([ROI_list(1).folder filesep ROI_list(1).name]));
mask(mask>=1)=1;

save_dir = 'infer_geo';

%% IC analysis, data preparation
% load and save ROI data

for subs=1:length(subdir)
    subs
    raw_vols = spm_read_vols(spm_vol([pwd filesep 'data' filesep subdir{subs} filesep 'whole_brain_average.nii']));%**********
    for conds = 1:size(raw_vols,4)
        raw_vols_conds = raw_vols(:,:,:,conds);
        mask_inds = find(mask>=1);
        mask_data(:,conds) = raw_vols_conds(mask>=1);
    end
    save([save_dir filesep subdir{subs} '_dACC_conf.mat'], 'mask_inds', 'mask_data', '-v7.3');%notice:**********
    clear mask_inds mask_data raw_vols
end

% load the mask for selecting the voxels in mask only
% select the mask you need to do mvpa ******
save_dir = 'infer_geo';

%% for predicting MRI RDM with theroical behavioral RDM, on the group level
%construct the list for the combination of theta & scale, and theorical
%behavioral RDM
iters=0;

params = {[90*pi/180, 1], [90*pi/180, 0], [90*pi/180, 1e+2], [84*pi/180, .9]};
for idx = 1:length(params) %[-pi, -90*pi/180, -45*pi/180, 0, 45*pi/180 ,90*pi/180, pi] % 0:pi/180:pi
        iters                = iters+1;
        iters_all(iters,1:2) = params{idx};
        if idx==3
            beha_d_single = pdist([1 2 3 1 2 3 1 2 3]', 'euclidean');
        else
            beha_d_single = SolveCosFunction([iters_all(iters,1),iters_all(iters,2)]);
        end
        imagesc(squareform(beha_d_single));
        beha_d_all(:,iters)        = repmat((beha_d_single)',29,1);
        subj_c = (repelem(1:29,36)'); % 30 is the length(subs), 36 is length(RDM)
end



%% extract the neural data and construct the neural RDM
    addpath(genpath('/opt/matlab-special-heatmap'));
    colorlist = flipud(slanCM(21));
for subs = 1:length(subdir)
    subs
% ===================== for fMRI ==========================================
        load([save_dir filesep subdir{subs} '_dACC_conf.mat']);
% ===================== for RNN ===========================================
%         load('RNN/RNN_midLayer_conds_9.mat');
%         mask_data = reshape(conds_9(subs,:,:), size(conds_9,2), size(conds_9,3));
% ===================== for RNN =========================================== 
    for rois = 1
        mask_data           = zscore(mask_data);
        mask_data_all(subs,:,:) = mask_data;
% % %             set(gcf,'unit','pixels','position',[500 50 200 200])
% % %             neur_cond_RDM = pdist((mask_data)','Euclidean')';
% % %             imagesc(squareform((neur_cond_RDM)));
% % %             colormap(colorlist(1:end,:))
% % % %             axis off
% % %             print(['NeeuralRDM_' num2str(subs)], '-dsvg','-painters');
        neur_cond_RDM       = pdist((mask_data)', 'mahal',covdiag((mask_data)'));
        neur_RDM((subs-1)*36+1:subs*36,1) = (neur_cond_RDM');
    end %concact the RDMs of each subject [subs*lenth(RDM),1]
end


%% data diagnosis ---------------------------------------------------------
% % % % % % VarNames = {'x','y','subj'};
% % % % % % tbl = array2table([beha_d_all(:,190),neur_RDM,subj_c],'VariableNames',VarNames);
% % % % % % y=neur_RDM; x=beha_d_all(:,190);
% % % % % % 
% % % % % % figure; set(gcf,'unit','pixels','position',[500 50 350 1200]); hold on
% % % % % % g = unique(subj_c);
% % % % % % for i = 1:numel(g)
% % % % % %     idx = subj_c==g(i);
% % % % % %     scatter(x(idx), y(idx), 30, 'filled'); % scatter by each group
% % % % % %     % fit linear regression
% % % % % %     if sum(idx)>1
% % % % % %         b = polyfit(x(idx), y(idx), 1);
% % % % % %         xx = linspace(min(x(idx)), max(x(idx)),100)';
% % % % % %         yy = polyval(b, xx);
% % % % % %         plot(xx, yy, 'LineWidth',1.5);
% % % % % %     end
% % % % % % end
% % % % % % hold off
% % % % % % xlabel('X'); ylabel('Y'); legend(num2str(g),'location','best');
% % % % % % title('scatter for subs and fitting');
% % % % % % % 2. fit LME
% % % % % % lmeREML = fitlme(tbl, 'y ~ x + (1|subj)', 'FitMethod','REML');
% % % % % % % 3. residual diagnosis
% % % % % % fittedVals = fitted(lmeREML); % size N×1
% % % % % % residVals = residuals(lmeREML); % size N×1
% % % % % % 
% % % % % % figure;
% % % % % % scatter(fittedVals, residVals, 12, 'filled');
% % % % % % refline(0,0);
% % % % % % xlabel('Fitted value');
% % % % % % ylabel('residual');
% % % % % % title('fit value vs residual for unequal variance');
% % % % % % 
% % % % % % % 3.2 QQ plot
% % % % % % figure;       
% % % % % % qqplot(residVals);
% % % % % % title('QQ plot for normal');
% % % % % % % 3. model comparation
% % % % % % lme1 = fitlme(tbl, 'y ~ 1+(1|subj)', 'FitMethod','ML','CovariancePattern','Isotropic');
% % % % % % lme2 = fitlme(tbl, 'y ~ 1+x + (1+x|subj)', 'FitMethod','ML','CovariancePattern','Isotropic');
% % % % % % cmp = compare(lme1, lme2);
% % % % % % disp(cmp(:,{'Model','AIC','BIC','pValue'}))
% data diagnosis ---------------------------------------------------------

%% run the LME
% for subs=1:29
parfor iters = 1:length(iters_all)
       beha_d = (beha_d_all(:,iters));
       neur_cond_RDM = (neur_RDM(:,1));
        % testing the behavioral strategy model:
        % *hierarchical model
        % *win stay lose swi: theta=0, scale=1/4; (iters_all(26,:)) ???
        % *delayed swi (N-ers based swi):theta=90, scale=2; (iters_all(18091,:)) ???
        % *pcpt based swi: theta=0, scale=0; (iters_all(1,:)) ???
       
        % run linear mix model
        VarNames = {'beha_d','neur_cond_RDM','subj_c'};
        tbl = array2table([beha_d,neur_cond_RDM,subj_c],'VariableNames',VarNames);
        formula = 'neur_cond_RDM ~ 1 + beha_d + (1 + beha_d | subj_c)';
        % fixed intercept for fixed effects, random intercept and slope for rand
        % vars
        fit_results = fitlme(tbl,formula);%******
        RMSE{iters,1} = sqrt(mean(residuals(fit_results).^2));

        LogLikelihood{iters,1} = -fit_results.ModelCriterion.LogLikelihood;
        beta{iters,1}          = fit_results.Coefficients.Estimate(2);
        pvalue{iters,1}        = fit_results.Coefficients.pValue(2);
        AIC{iters,1}           = fit_results.ModelCriterion.AIC;
        BIC{iters,1}           = fit_results.ModelCriterion.BIC;
end
% end

LogLikelihood =cell2mat(LogLikelihood)
beta          =cell2mat(beta); 
pvalue        =cell2mat(pvalue)
AIC           =cell2mat(AIC);
BIC           =cell2mat(BIC);
RMSE          =cell2mat(RMSE);

%% plot the dynamic of theta & scale from neural data
%find min LogLikelihood
% for subs =1:29
    addpath(genpath('/opt/matlab-special-heatmap'));
    set(0, 'DefaultAxesFontName', 'Helvetica'); % Fonttype for axis
    set(0, 'DefaultTextFontName', 'Helvetica'); % Fonttype for text
    set(gcf,'unit','pixels','position',[500 50 900 200])
    
colorlist = flipud(slanCM(21));
title_rdm = {'Orthogonal','Errors-based','Difficulty-based'}
for idx = 1:3
    subplot(1,4,idx)
    set(gca,'Fontsize',12)
       if idx==3
            d = pdist([1 2 3 1 2 3 1 2 3]', 'euclidean');
       else
            [X, d]= SolveCosFunction_all(iters_all(idx,:));
       end
        imagesc(squareform(d)); axis off;
        colormap(colorlist(1:end,:))
%         colorbar
        hold on
        title(title_rdm{idx},'Fontsize',12,'Fontweight','normal');
end
hold off
% print(['Geometry_fMRI_Conf_dACCmvpa_modelCom1'], '-dsvg','-painters');

% MDS----------------------------------------------------------------------
% % % clear Y
% % % mask_data_all_b= zscore(mask_data_all,0,3);
% % % data_RDM = reshape(mean(mask_data_all_b,1),size(mask_data_all_b,3),1*size(mask_data_all_b,2));
% % % 
% % %     rng(4) %rng(4)
% % % DisMatrix = squareform(pdist(data_RDM, 'mahal',covdiag(data_RDM)));
% % % [Y, stress,dispa] = mdscale(DisMatrix,3,'criterion','metricsstress','start','random');
% % % 
% % % subplot(2,3,4)
% % % for p_i = 1:size(Y,1)
% % % s=scatter(Y(p_i,1),Y(p_i,2),'filled');
% % % s.SizeData = 300; s.MarkerFaceColor=colorlist(p_i*15,:); hold on
% % % text(Y(p_i,1),Y(p_i,2),sprintf('%.0f %.0f',round(p_i)))
% % % end
% % % % axis equal
% % % plot(Y([3 2 1],1),Y([3 2 1],2),'-','LineWidth',3','Color',colorlist(60,:)); 
% % % plot(Y([4 6 5],1),Y([4 6 5],2),'-','LineWidth',3','Color',colorlist(110,:)); 
% % % plot(Y([7 8 9],1),Y([7 8 9],2),'-','LineWidth',3','Color',colorlist(180,:));
% % % 
% % % plot(Y([2 6 8],1),Y([2 6 8],2),'--','LineWidth',3','Color',colorlist(10,:));
% % % grid off; xlabel('MDS1'); ylabel('MDS2');
% % % hold off
% MDS----------------------------------------------------------------------

% % % subplot(2,3,5)
% % %     hbar = bar(AIC(1:3)); 
% % %     hbar.FaceColor= colorlist([170],:);
% % %     hbar.EdgeColor= colorlist([170],:);
% % %     hbar.FaceAlpha= 1;
% % %     ylabel('AIC'); xlabel('Models'); xticklabels({'Orth','Ners','Difs'})
% % %     title('Model Compare'); 
subplot(1,4,4)
    set(gca,'Fontsize',12)
    hbar = bar(AIC(1:3)-AIC(4)); 
    hbar.FaceColor= colorlist([190],:);
    hbar.EdgeColor= colorlist([190],:);
    hbar.FaceAlpha= 1;
    ylabel('Delta AIC'); xlabel('Models'); xticklabels({'Orth','Err','Diff'});
%     xtickangle(30)
    title('Model Compare to(84\circ, 0.9)','Fontsize',12,'Fontweight','normal');
    box off
print(['Geometry_fMRI_Conf_dACCmvpa_modelCom2'], '-dsvg','-painters');






% subplot(2,3,5)

%% MDS
clear Y
%%%%%% mds1
for subs=1:29
    data_RDM = reshape(mask_data_all(subs,:,:),size(mask_data_all,2),1*size(mask_data_all,3));
%     [coeff,Ysub, ~,~,~] = pca(data_RDM); %pca(data_RDM,'NumComponents',3)
%     data_RDM(:,1:3)=-(data_RDM(:,1:3));
    DisMatrix1(subs,:,:) = squareform(pdist(data_RDM', 'mahal',covdiag(data_RDM')));
end
DisMatrix       = reshape(mean(DisMatrix1,1),size(DisMatrix1,2),1*size(DisMatrix1,3));%([3:6 9 13 16:17 19:21 24 26:27 29],:,:)
rng(3) %rng(4)
[Y, stress,dispa] = mdscale(DisMatrix,3,'criterion','sstress',...
    'start','random');


%%%%%% mds2
% % % mask_data_all_b= zscore(mask_data_all,0,3);
% % % data_RDM = reshape(mean(mask_data_all,1),size(mask_data_all,2),1*size(mask_data_all,3));
% % % 
% % % DisMatrix = squareform(pdist(data_RDM', 'mahal',covdiag(data_RDM')));
% % % % DisMatrix = squareform((pdist(data_RDM, 'Spearman'))); %Spearman
% % % rng(4) %rng(4)
% % % [Y, stress,dispa] = mdscale(DisMatrix,3,'criterion','sammon',...
% % %     'start','random');

%%%%%%%%%%% for PCA
% % %     [coeff,Y, ~,~,~] = pca(data_RDM'); %pca(data_RDM,'NumComponents',3)

%%%%%%for orther toolbox
% % % RDMs.RDM = DisMatrix;
% % % RDMs.name = [];
% % % RDMs.color = [0 0 1];
% % % RDMs.name = 'rsa';
% % %   userOptions.analysisName = 'rsa';
% % %   userOptions.rootPath = '/mnt/HR_project_SZU/preprocess_by_spm12/tRSA_confidence_9cond_fan_searchlight_10mm_p';
% % % localOptions.dotColours = [0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;0 0 1;0 0 1]
% % % 
% % % MDSConditions(RDMs,userOptions,localOptions);

%%%%%%for setting start points
% % %     x1 = linspace(-1, 1, ceil(sqrt(9)));
% % %     y1 = linspace(-1, 1, ceil(sqrt(9)));
% % %     [x2,y2] = meshgrid(x1, y1);
% % %     custom_start = [x2(:), y2(:)];
% % %     custom_start = custom_start(1:9,:);
% % % % theta = linspace(0, 2*pi, 9)'  
% % % % custom_start = [cos(theta), sin(theta)];
% % % [Y] = mdscale(DisMatrix,2,'start',custom_start);


% % % % % % for best-------------------------------------------------------
DisMatrix = squareform(SolveCosFunction([98*pi/180, 1])); 
% % % RNN[98*pi/180, 1]; dACC[84*pi/180, 0.9] 
[Y] = mdscale(DisMatrix,3);


% Y=zscore(Y);
figure; set(gcf,'unit','pixels','position',[500 50 1200 200])
% 3-D
colorlist_ners = flip(slanCM(2));%ner: 2 dif :77
colorlist_difs = flip(slanCM(77));%ner: 2 dif :77
colorlist_difs=[colorlist_difs(20,:); colorlist_difs(50,:); colorlist_difs(90,:);...
    colorlist_difs(20,:); colorlist_difs(50,:); colorlist_difs(90,:);...
    colorlist_difs(20,:); colorlist_difs(50,:); colorlist_difs(90,:);];

subplot(1,4,1)
    set(gca,'Fontsize',12)
% axis equal
points_lab = {'E0D1','E0D2','E0D3','E1D1','E1D2','E1D3','E2D1','E2D2','E2D3'}
for p_i = 1:size(Y,1)
s=scatter(Y(p_i,1),Y(p_i,2),'filled');
s.SizeData = 100; s.MarkerFaceColor=colorlist_difs(p_i,:)
hold on
        if p_i==1 | p_i==4 | p_i==7
            text_dis1= -0.4; text_dis2 = 0.4;
        else
            text_dis1= -0.4; text_dis2 = -0.4;
        end
    text(Y(p_i,1)+text_dis1,Y(p_i,2)+text_dis2, points_lab{p_i},'Fontsize',8)
end
grid off

% % % fill3([Y(1,1),Y(2,1),Y(3,1)],...
% % % [Y(1,2),Y(2,2),Y(3,2)],...
% % % [Y(1,3),Y(2,3),Y(3,3)],colorlist_ners(20,:),'EdgeColor','none')
% % % fill3([Y(4,1),Y(5,1),Y(6,1)],...
% % % [Y(4,2),Y(5,2),Y(6,2)],...
% % % [Y(4,3),Y(5,3),Y(6,3)],colorlist_ners(50,:),'EdgeColor','none')
% % % fill3([Y(7,1),Y(8,1),Y(9,1)],...
% % % [Y(7,2),Y(8,2),Y(9,2)],...
% % % [Y(7,3),Y(8,3),Y(9,3)],colorlist_ners(90,:),'EdgeColor','none')

plot(Y([3 2 1],1),Y([3 2 1],2),'-','LineWidth',2','Color',colorlist_ners(20,:)); 
plot(Y([6 5 4],1),Y([6 5 4],2),'-','LineWidth',2','Color',colorlist_ners(50,:)); 
plot(Y([9 8 7],1),Y([9 8 7],2),'-','LineWidth',2','Color',colorlist_ners(90,:));

% plot(Y([2 5 8],1),Y([2 5 8],2),'--','LineWidth',2','Color',colorlist_difs(2,:));
% plot(Y([3 6 9],1),Y([3 6 9],2),'--','LineWidth',2','Color',colorlist_difs(3,:));
% plot(Y([1 4 7],1),Y([1 4 7],2),'--','LineWidth',2','Color',colorlist_difs(1,:));
grid off

xlabel('MDS dimension1'); ylabel('MDS dimension2');
xlim([-2 2]);ylim([-2 2])
view(2)
hold off


% 2-D
subplot(1,4,2)
for p_i = 1:size(Y,1)
s=scatter(Y(p_i,1),Y(p_i,2),'filled');
s.SizeData = 100; s.MarkerFaceColor=colorlist_difs(p_i,:)
hold on
text(Y(p_i,1),Y(p_i,2),sprintf('%.0f %.0f',round(p_i)))
end
% axis equal
plot(Y([3 2 1],1),Y([3 2 1],2),'-','LineWidth',2','Color',colorlist_ners(20,:)); 
plot(Y([6 5 4],1),Y([6 5 4],2),'-','LineWidth',2','Color',colorlist_ners(50,:)); 
plot(Y([9 8 7],1),Y([9 8 7],2),'-','LineWidth',2','Color',colorlist_ners(90,:));

plot(Y([2 5 8],1),Y([2 5 8],2),'--','LineWidth',2','Color',colorlist_difs(2,:));
grid off
% legend('Hard','Mid','Easy','0Er','1E','2E','Diff','Ners')
xlabel('MDS1'); ylabel('MDS2');
hold off

% 2-D
subplot(1,4,3)
for p_i = 1:size(Y,1)
s=scatter(Y(p_i,1),Y(p_i,3),'filled');
s.SizeData = 100; s.MarkerFaceColor=colorlist_difs(p_i,:)
hold on
text(Y(p_i,1),Y(p_i,3),sprintf('%.0f %.0f',round(p_i)))
end
% axis equal
plot(Y([3 2 1],1),Y([3 2 1],3),'-','LineWidth',2','Color',colorlist_ners(20,:)); 
plot(Y([6 5 4],1),Y([6 5 4],3),'-','LineWidth',2','Color',colorlist_ners(50,:)); 
plot(Y([9 8 7],1),Y([9 8 7],3),'-','LineWidth',2','Color',colorlist_ners(90,:));

plot(Y([2 5 8],1),Y([2 5 8],3),'--','LineWidth',2','Color',colorlist_difs(2,:));
grid on
xlabel('MDS1'); ylabel('MDS3');
hold off

% 2-D
subplot(1,4,4)
for p_i = 1:size(Y,1)
s=scatter(Y(p_i,2),Y(p_i,3),'filled');
s.SizeData = 100; s.MarkerFaceColor=colorlist_difs(p_i,:)
hold on
text(Y(p_i,2),Y(p_i,3),sprintf('%.0f %.0f',round(p_i)))
end

plot(Y([3 2 1],2),Y([3 2 1],3),'-','LineWidth',2','Color',colorlist_ners(20,:)); 
plot(Y([6 5 4],2),Y([6 5 4],3),'-','LineWidth',2','Color',colorlist_ners(50,:)); 
plot(Y([9 8 7],2),Y([9 8 7],3),'-','LineWidth',2','Color',colorlist_ners(90,:));

plot(Y([2 5 8],2),Y([2 5 8],3),'--','LineWidth',2','Color',colorlist_difs(2,:));
grid on
xlabel('MDS2'); ylabel('MDS3');
hold off

print(['GLM2_Dev_Ners_Geome_MDS_fMRI_dACC_mvpa-bestRNN'], '-dsvg','-painters');

%% for cal mahalonobis
 function sigma=covdiag(x)

    [t,n]=size(x);
    meanx=mean(x);
    x=x-meanx(ones(t,1),:);

    sample=(1/t).*(x'*x);

    prior=diag(diag(sample));

    d=1/n*norm(sample-prior,'fro')^2;
    y=x.^2;
    r2=1/n/t^2*sum(sum(y'*y))-1/n/t*sum(sum(sample.^2));

    shrinkage=max(0,min(1,r2/d));
    sigma=shrinkage*prior+(1-shrinkage)*sample;
 end
 
 function rdms = create_rdms_from_data(data)
    % data: 9×45×10  (条件×特征×被试)
    [n_conditions, n_features, n_subjects] = size(data);
    rdms_cell = cell(n_subjects, 1);
    
    for sub = 1:n_subjects
        % 提取当前被试的数据: 9条件 × 45特征
        subject_data = squeeze(data(:, :, sub)); % 9×45
        
        % 计算RDM矩阵 (9×9)
        rdm_matrix = pdist(subject_data, 'correlation'); % 或 'euclidean'
        rdm_matrix = squareform(rdm_matrix);
        
        % 创建RDM对象
        rdm_obj = rsa.rdm.RDMs();
        rdm_obj.RDM = rdm_matrix;
        rdm_obj.name = sprintf('Subject_%02d', sub);
        rdm_obj.color = [0, 0, 0];
        
        rdms_cell{sub} = rdm_obj;
    end
    
    % 合并所有被试的RDMs
    rdms = rsa.rdm.concat(rdms_cell);
end