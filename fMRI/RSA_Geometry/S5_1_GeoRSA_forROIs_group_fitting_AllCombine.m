%% setup
clear; clc
restoredefaultpath
addpath('/opt/spm12')
%% set dir
matdir=dir('/mnt/HR_project_SZU/preprocess_by_spm12/tRSA_confidence_9cond_fan_searchlight_10mm_p/data'); %read the list of subjects
subdir={matdir([matdir.isdir]).name};
subdir=subdir(~ismember(subdir,{'.','..'}));

% ROI_name = arrayfun(@(i) sprintf('Schaefer%03d',i),1:400,'UniformOutput',false);
ROI_list = dir('mask/m01_mu_p_FWEc227_mask.nii'); %mu_p_FWEc227_mask mask_mvpa_ner_dif
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

for theta = 0:1*pi/180:pi %[-pi, -90*pi/180, -45*pi/180, 0, 45*pi/180 ,90*pi/180, pi] % 0:pi/180:pi
for scale = 0:0.1:2 % scale = 0:0.01:2
        iters                = iters+1;
        iters_all(iters,1:2) = [theta scale];
        beha_d_single = SolveCosFunction([iters_all(iters,1),iters_all(iters,2)]);
        beha_d_all(:,iters)        = repmat((beha_d_single)',29,1);
        subj_c = (repelem(1:29,36)'); % 30 is the length(subs), 36 is length(RDM)
end
end

%% extract the neural data and construct the neural RDM
for subs = 1:length(subdir)
    subs
% ===================== for fMRI ==========================================
%         load([save_dir filesep subdir{subs} '_dACC_conf.mat']);
% ===================== for RNN ===========================================
        load('RNN/RNN_midLayer_conds_9.mat');
        mask_data = (reshape(conds_9(subs,:,:), size(conds_9,2), size(conds_9,3)))';%***
% ===================== for RNN ===========================================        
    for rois = 1
        mask_data = zscore(mask_data);
%         [coeff, scores, latent, ~, explained]=pca(mask_data);
%         mask_data = scores(:,1:9); 
%         neur_cond_RDM = pdist((mask_data)','correlation')';
        neur_cond_RDM = pdist((mask_data'), 'mahal',covdiag((mask_data')));
        imagesc(squareform((neur_cond_RDM)));
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
       iters
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

LogLikelihood =cell2mat(LogLikelihood); 
beta          =cell2mat(beta); 
pvalue        =cell2mat(pvalue);
AIC           =cell2mat(AIC);
BIC           =cell2mat(BIC);
RMSE          =cell2mat(RMSE);

%% plot the dynamic of theta & scale from neural data
%find min LogLikelihood
% for subs =1:29

min_RDMpredict_log{1} = find( min((LogLikelihood4(:,1)))==(LogLikelihood(:,1)));
min_RDMpredict_p(1)   = (min(pvalue(min_RDMpredict_log{1},1)));
min((LogLikelihood(:,1)))

min_RDMpredict_p(min_RDMpredict_p==0) = nan;
min_RDMpredict_p_fdr = pval_adjust(min_RDMpredict_p, 'bonferroni')

if min_RDMpredict_p_fdr<0.05
angle_rois(1) = mean(iters_all([min_RDMpredict_log{1}],1))
scale_rois(1) = mean(iters_all([min_RDMpredict_log{1}],2))
else
    fprintf(":no significant results")
end


%%
    addpath(genpath('/opt/matlab-special-heatmap'));
    color_list = [colormap(flipud(slanCM(112))); colormap(flipud(slanCM(97)))];%97
    set(gcf,'unit','pixels','position',[500 50 400 350])
    set(0, 'DefaultAxesFontName', 'Helvetica'); % Fonttype for axis
    set(0, 'DefaultTextFontName', 'Helvetica'); % Fonttype for text
% subplot(1,1,1)
colorlist = flipud(slanCM(21));
for scale = mean(scale_rois)
    for theta = mean(angle_rois)
        [X, d]= SolveCosFunction_all([theta, scale]);
        imagesc(squareform(d));
        set(gca,'Fontsize',12)
        colormap(colorlist(1:end,:))
        colorbar
        yticks([1:1:9]); yticklabels({'E0D1','E0D2','E0D3','E1D1','E1D2','E1D3','E2D1','E2D2','E2D3'});
        xticks([1:1:9]); xticklabels({'E0D1','E0D2','E0D3','E1D1','E1D2','E1D3','E2D1','E2D2','E2D3'});
        xtickangle(45)
        hold on
        title(['scale=' num2str(round(scale,2)) ';' 'Angle=' num2str(round(theta.*180/pi,2))],...
        'Fontsize',14,'Fontweight','normal');
    end
end
hold off

print(['Geometry_ROIs_RNN_Conf_Mah.svg'], '-dsvg','-painters');

%%
z=[];
iters=0;
for theta = 1:size(0:1*pi/180:pi,2) %x axis
for scale = 1:size(0:0.1:2,2) %y axis
    iters=iters+1;
    z(theta,scale)=LogLikelihood(iters);
end
end
z = imresize(z, 5, 'bilinear');

figure
set(gcf,'unit','pixels','position',[500 50 500 350])
set(0, 'DefaultAxesFontName', 'Helvetica'); % Fonttype for axis
set(0, 'DefaultTextFontName', 'Helvetica'); % Fonttype for text
    imagesc(z)
    set(gca,'Fontsize',12)
    hold on
    xticks(1:10*5:21*5); xticklabels({'0','1.57','\pi'}); xlabel('Theta');
    yticks(1:90*5:181*5); yticklabels({'0','1','2'}); ylabel('Scale')
    set(gca,'LineWidth',1.5)
    colormap(flip(slanCM(9)))
    c=colorbar; caxis([min(z,[],'all') max(z,[],'all')])
    c.Ticks = [min(z,[],'all') max(z,[],'all')]; c.TickLabels = {'Min','Max'};
    title('Fitting results (-log-likelihood)','Fontsize',14,'Fontweight','normal')
    
    annotation('arrow',[.92 .92], [.9 .17],'LineWidth',1.5) %arrows
    text(1.13, .5,'Better','LineWidth',1.5, 'Units','normalized',...
        'HorizontalAlignment','left','VerticalAlignment','top','rotation',90) %arrows 
    
hold off
print(['Geometry_parameters_RNN.svg'], '-dsvg','-painters');


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