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

params = {[90*pi/180, 1], [90*pi/180, 0], [90*pi/180, 1e+2],[43*pi/180, .7]};
for idx = 1:length(params) %[-pi, -90*pi/180, -45*pi/180, 0, 45*pi/180 ,90*pi/180, pi] % 0:pi/180:pi
        iters                = iters+1;
        iters_all(iters,1:2) = params{idx};
        if idx==3
            beha_d_single = pdist([1 2 3 1 2 3 1 2 3]', 'euclidean');
        else
            beha_d_single = SolveCosFunction([iters_all(iters,1),iters_all(iters,2)]);
        end
        imagesc(squareform(beha_d_single));
        beha_d_all(:,iters)        = repmat((beha_d_single)',1,1);
        subj_c = (repelem(1:1,36)'); % 30 is the length(subs), 36 is length(RDM)
end



%% extract the neural data and construct the neural RDM
for subs = 1:length(subdir)
    subs
% ===================== for fMRI ==========================================
        load([save_dir filesep subdir{subs} '_dACC_conf.mat']);
% ===================== for RNN ===========================================
%         load('RNN/RNN_midLayer_conds_9.mat');
%         mask_data = reshape(conds_9(subs,:,:), size(conds_9,3), size(conds_9,2));
% ===================== for RNN =========================================== 
    for rois = 1
        mask_data=zscore(mask_data);
%         neur_cond_RDM = pdist((mask_data)','correlation')';
        neur_cond_RDM = pdist((mask_data)', 'mahal',covdiag((mask_data)'));
        imagesc(squareform((neur_cond_RDM)));
        neur_RDM{subs}(:,1) = (neur_cond_RDM');
    end %concact the RDMs of each subject [subs*lenth(RDM),1]
end

%% run the LME
for subs=1:29
    subs
for iters = 1:length(iters_all)
       beha_d = (beha_d_all(:,iters));
       neur_cond_RDM = (neur_RDM{subs});
        % testing the behavioral strategy model:
        % *hierarchical model
        % *win stay lose swi: theta=0, scale=1/4; (iters_all(26,:)) ???
        % *delayed swi (N-ers based swi):theta=90, scale=2; (iters_all(18091,:)) ???
        % *pcpt based swi: theta=0, scale=0; (iters_all(1,:)) ???
       
        % run linear mix model
        VarNames = {'beha_d','neur_cond_RDM','subj_c'};
        tbl = array2table([beha_d,neur_cond_RDM,subj_c],'VariableNames',VarNames);
        formula = 'neur_cond_RDM ~ 1 + beha_d';
        % fixed intercept for fixed effects, random intercept and slope for rand
        % vars
        fit_results = fitlme(tbl,formula);%******
        RMSE{iters,subs} = sqrt(mean(residuals(fit_results).^2));

        LogLikelihood{iters,subs} = -fit_results.ModelCriterion.LogLikelihood;
        beta{iters,subs}          = fit_results.Coefficients.Estimate(2);
        pvalue{iters,subs}        = fit_results.Coefficients.pValue(2);
        AIC{iters,subs}           = fit_results.ModelCriterion.AIC;
        BIC{iters,subs}           = fit_results.ModelCriterion.BIC;
end
end

LogLikelihood =cell2mat(LogLikelihood)
beta          =cell2mat(beta); 
pvalue        =cell2mat(pvalue)
AIC           =cell2mat(AIC);
BIC           =cell2mat(BIC);
RMSE          =cell2mat(RMSE);


%% plot the theta & scale from neural data
%find min LogLikelihood
% for subs =1:29
    addpath(genpath('/opt/matlab-special-heatmap'));
    set(gcf,'unit','pixels','position',[500 50 850 550])
colorlist = flipud(slanCM(21));
title_rdm = {'Orth','Ners','Difs'}
for idx = 1:3
    subplot(2,3,idx)
        if idx==3
            d = pdist([1 2 3 1 2 3 1 2 3]', 'euclidean');
        else
            [X, d]= SolveCosFunction_all(iters_all(idx,:));
        end
        imagesc(squareform(d));
        colormap(colorlist(1:end,:))
%         colorbar
        hold on
        title(title_rdm{idx});
end

bar1 = subplot(2,3,4)
    hbar = bar(mean(beta(1:3,:),2)); 
    hbar.FaceColor= colorlist([10],:);
    hbar.EdgeColor= colorlist([10],:);
    hbar.FaceAlpha= 1;
    hold on
    scatter((1:29)*0.01+0.9,beta(1,:))
    scatter((1:29)*0.01+0.9*2,beta(2,:))
    scatter((1:29)*0.01+0.9*3,beta(3,:))
    ylabel('Beta'); xlabel('Models'); xticklabels({'Orth','Ners','Difs'})
    title('Neural representational strength','Fontsize',8); 
        set(bar1, 'Position', [0.15 0.1 0.25 0.35]); %[left bottom width height]
    hold off conf_fitmod.RMSE; %

bar2=subplot(2,3,5)
    hbar = bar(mean(beta(1:3,:),2)-mean(beta(4,:),2)); 
    hbar.FaceColor= colorlist([10],:);
    hbar.EdgeColor= colorlist([10],:);
    hbar.FaceAlpha= 1;
    hold on
    scatter((1:29)*0.01+0.9,beta(1,:)-beta(4,:))
    scatter((1:29)*0.01+0.9*2,beta(2,:)-beta(4,:))
    scatter((1:29)*0.01+0.9*3,beta(3,:)-beta(4,:))
    ylabel('Delta Beta'); xlabel('Models'); xticklabels({'Orth','Ners','Difs'})
    title('Neural representational strength','Fontsize',8);
    set(bar2, 'Position', [0.65 0.1 0.25 0.35]);
    
hold off
print(['Geometry_fMRI_Conf_dACCmvpa_modelCom_induvidual.tif'], '-dtiff', '-r300','-opengl');



%% correlation
% prepare the accuracy of subjects and plot the correlation between CosScale
% and accuracy
     beha_data=readtable('/mnt/HR_project_SZU/preprocess_by_spm12/data_29sub_raw.xlsx');
     sub_index=table2array(unique(beha_data(:,1)));
     load('/mnt/HR_project_SZU/preprocess_by_spm12/model_res_switch_fmri_m_4_29sub.mat');
    for subs=1:length(sub_index)
        indics=[];
            beha_data1        = beha_data(beha_data.idIdx == sub_index(subs) & beha_data.expIdx==1,:);%only inference task
            beha_data1        = rmmissing(beha_data1);%%******
            infer_acc(subs,1) = sum(beha_data1.feedback==1)/length(beha_data1.feedback);
            infer_err_swi(subs,1) = sum(beha_data1.ruleSwitch(1:end-1).* (beha_data1.feedback(2:end)))/sum(beha_data1.ruleSwitch(1:end-1));
            infer_err_swi_n(subs,1) = sum(beha_data1.ruleSwitch==1);
            infer_rule(subs,1)= sum(beha_data1.rulelist == beha_data1.ruleResp)/length(beha_data1.rulelist);  
            infer_pcpt(subs,1)= sum(beha_data1.answer_angle == beha_data1.pcptResp)/length(beha_data1.pcptResp);
            
            %confidence related
            table_d = table;
            CBMdata = Input(cell2mat({Input.ID}) == sub_index(subs));
            table_d.tdev_n = CBMdata.tDev(CBMdata.StateExp==1);
            table_d.tdev_n(table_d.tdev_n==0.1 | table_d.tdev_n==0.9) = 3;
            table_d.tdev_n(table_d.tdev_n==0.3 | table_d.tdev_n==0.7) = 2;
            table_d.tdev_n(table_d.tdev_n==0.5) = 1;
            
            table_d.nback_n = CBMdata.Nback_act(CBMdata.StateExp==1);
            table_d.nback_n(table_d.nback_n>=2) = 2;
            
            table_d.conf_n = CBMdata.pr_of_switch(CBMdata.StateExp==1);
            conf_fitmod = fitlm(table_d,'conf_n ~ tdev_n + nback_n + tdev_n * nback_n');
            conf_beta(subs,1) =  conf_fitmod.RMSE % ;conf_fitmod.Rsquared.Ordinary ;% % the fitting extent of error and difficulty on P(swi)
    end
% for theta value
% x_cor = (abs(beta(1,:))-abs(beta(4,:)))';
x_cor = (abs(beta(1,:))-abs(beta(4,:)))';
[rvalues(1), pvalues(1)] = corr(x_cor, infer_acc,'Type','Pearson')
% [rvalues(2), pvalues(2)] = corr(infer_err_swi, y_cor,'Type','Pearson')

close all
set(gcf,'unit','pixels','position',[500 50 550 450])
    set(0, 'DefaultAxesFontName', 'Helvetica'); % Fonttype for axis
    set(0, 'DefaultTextFontName', 'Helvetica'); % Fonttype for text

% subplot(1,2,1)    
    mdl = fitlm(x_cor, infer_acc);
    x_fit = linspace(min(x_cor), max(x_cor), 200)'; % for independent variable, x
    [y_fit, y_ci] = predict(mdl, x_fit);
    hold on;
    fill([x_fit; flipud(x_fit)], [y_ci(:,1); flipud(y_ci(:,2))], ...
    colorlist(170,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);

    scatter(x_cor, infer_acc, 60, colorlist(170,:), 'filled');
    plot(x_fit, y_fit, 'b-', 'LineWidth', 4, 'Color',colorlist(170,:));

    legend({'95% CI','Data', ['r: ', num2str(round(rvalues(1),3)),...
    ' p: ', num2str(round(pvalues(1),3))]},...
    'Location','SouthEast','FontName','Helvetica','Fontsize',12);
%     legend boxoff
    xticks(round([min(x_cor),min(x_cor)/2,0, max(x_cor)/2, max(x_cor)],2))
    ylim([min(infer_acc)-0.05 max(infer_acc)+0.02])
    yticks(round([min(infer_acc)-0.05:0.06: max(infer_acc)],2))
    
    set(gca,'TickDir','out','LineWidth',1.2,'FontName','Helvetica');
    xlabel('Beta: Orth - Group best (43, 0.7)','Fontsize',12); 
    ylabel('Inference task accuracy','Fontsize',12);
    ax=gca; ax.TickDir ='out'; box off; ax.FontSize=12
    title([])
    
% subplot(1,2,2)    
%     mdl = fitlm(infer_err_swi,y_cor);
%     x_fit = linspace(min(infer_err_swi), max(infer_err_swi), 200)';
%     [y_fit, y_ci] = predict(mdl, x_fit);
%     hold on;
%     fill([x_fit; flipud(x_fit)], [y_ci(:,1); flipud(y_ci(:,2))], ...
%     colorlist(170,:), 'EdgeColor', 'none', 'FaceAlpha', 0.15);
% 
%     scatter(infer_err_swi, y_cor, 40, colorlist(170,:), 'filled');
%     plot(x_fit, y_fit, 'b-', 'LineWidth', 3,'Color',colorlist(170,:));
% 
%     legend({'95% CI','Data',['r: ', num2str(round(rvalues(2),3)),...
%     '  p: ', num2str(round(pvalues(2),3))]},'Location','SouthEast','Fontsize',12);
%     yticks(round([min(y_cor),min(y_cor)/2,0, max(y_cor)/2, max(y_cor)],2))
%     xlim([min(infer_err_swi)-0.01 max(infer_err_swi)+0.1])
%     xticks(round([min(infer_err_swi):0.1: max(infer_err_swi) max(infer_err_swi)+0.1],2))
%     
%     set(gca,'TickDir','out','LineWidth',1.2,'FontName','Helvetica');
%     ylabel('beta: Orth - Group best (43, 0.7)','Fontsize',12); 
%     xlabel('Inference task accuracy','Fontsize',12);
%     ylabel('beta: Orth - Group best (43, 0.7)','Fontsize',12); 
%     xlabel('Rule switch accuracy','Fontsize',12);
%     ax=gca; ax.TickDir ='out'; box off; ax.FontSize=12
%     title([])

print(['Geometry_fMRI_Conf_dACCall_corr1'], '-dsvg','-painters');


% ================== Subject-level meadian effect =========================
% X = orthogonality angle
% M = confidence quality
% Y = behavioral accuracy
% Standard mediation:
%   path a : M ~ X
%   path c : Y ~ X
%   path b,c': Y ~ X + M
% indirect effect = a*b

% significance by bootstrap CI
    Xz = zscore(x_cor);      Mz = zscore(conf_beta);      Yz = zscore(infer_acc);       
    nSub = length(Xz);
% 1) PATH a : M ~ X
    Xa = [ones(nSub,1), Xz];
    [b_a,~,~,~,stats_a] = regress(Mz, Xa);
    a  = b_a(2);
    R2_a = stats_a(1);
    p_a  = round(stats_a(3), 3);
% 2) TOTAL EFFECT c : Y ~ X
    Xc = [ones(nSub,1), Xz];
    [b_c,~,~,~,stats_c] = regress(Yz, Xc);
    c_total = b_c(2);
    R2_c = stats_c(1);
    p_c  = round(stats_c(3), 3);
% 3) PATH b and DIRECT EFFECT c' : Y ~ X + M
    Xb = [ones(nSub,1), Xz, Mz];
    [b_b,~,~,~,stats_b] = regress(Yz, Xb);
    c_prime = b_b(2);   % direct effect of X controlling M
    b_path  = b_b(3);   % effect of M controlling X
    R2_b = stats_b(1);
    F_b  = stats_b(2);
    p_model_b = round(stats_b(3),3 );
% ---- get p-values for each coefficient in Y ~ X + M ----
    tbl_b = array2table([Xz, Mz, Yz], 'VariableNames', {'X','M','Y'});
    mdl_b = fitlm(tbl_b, 'Y ~ X + M');

    coefTable = mdl_b.Coefficients;
    p_c_prime = round(coefTable.pValue(strcmp(coefTable.Properties.RowNames, 'X')), 3);
    p_b_path  = round(coefTable.pValue(strcmp(coefTable.Properties.RowNames, 'M')), 3);
% 4) INDIRECT EFFECT = a*b
    ab = a * b_path;
% 5) BOOTSTRAP INDIRECT EFFECT
    nBoot = 10000;
    boot_ab = nan(nBoot,1);
    rng(123); % reproducible

for ib = 1:nBoot
    idx = randsample(nSub, nSub, true);
    X_boot = Xz(idx);
    M_boot = Mz(idx);
    Y_boot = Yz(idx);
    % a path
    Xa_boot = [ones(nSub,1), X_boot];
    ba_boot = regress(M_boot, Xa_boot);
    a_boot = ba_boot(2);
    % b path
    Xb_boot = [ones(nSub,1), X_boot, M_boot];
    bb_boot = regress(Y_boot, Xb_boot);
    b_boot = bb_boot(3);

    boot_ab(ib) = a_boot * b_boot;
end

    CI_ab = prctile(boot_ab, [2.5 97.5]);
    p_ab_twosided = 2 * min(mean(boot_ab >= 0), mean(boot_ab <= 0));

% 6) PROPORTION MEDIATED
% be cautious if c_total is very small or sign changes
    prop_mediated = ab / c_total;

% ==========================
% 9) PRINT RESULTS
% ==========================
fprintf('\n===== SUBJECT-LEVEL MEDIATION RESULTS =====\n');

fprintf('\nPath a: M ~ X\n');
fprintf('a = %.4f, p = %.6f, R2 = %.4f\n', a, p_a, R2_a);

fprintf('\nTotal effect c: Y ~ X\n');
fprintf('c = %.4f, p = %.6f, R2 = %.4f\n', c_total, p_c, R2_c);

fprintf('\nPath b and direct effect c'': Y ~ X + M\n');
fprintf('b  = %.4f, p = %.6f\n', b_path, p_b_path);
fprintf('c'' = %.4f, p = %.6f\n', c_prime, p_c_prime);
fprintf('Model R2 = %.4f, model p = %.6f\n', R2_b, p_model_b);

fprintf('\nIndirect effect:\n');
fprintf('a*b = %.4f\n', ab);
fprintf('Bootstrap 95%% CI = [%.4f, %.4f]\n', CI_ab(1), CI_ab(2));
fprintf('Bootstrap p (approx, two-sided) = %.6f\n', p_ab_twosided);

fprintf('\nProportion mediated (ab/c) = %.4f\n', prop_mediated);

if CI_ab(1) > 0 || CI_ab(2) < 0
    fprintf('Conclusion: significant mediation (bootstrap CI does not include 0).\n');
else
    fprintf('Conclusion: indirect effect not significant (bootstrap CI includes 0).\n');
end


%% 10) FIGURE ==============================================================
fig = figure('Color','w','Position',[80 100 1500 400]);
    set(0, 'DefaultAxesFontName', 'Helvetica'); % Fonttype for axis
    set(0, 'DefaultTextFontName', 'Helvetica'); % Fonttype for text
% ---------- style ----------
dotC   = colorlist(170,:); %[0.20 0.20 0.20];
lineC  = colorlist(170,:); lineC2 = [0.10 0.10 0.10];
shadeC = colorlist(170,:); %[0.85 0.85 0.85];
nsC    = colorlist(170,:); nsC2 =[0.50 0.50 0.50];

ms   = 60;
lw   = 1.6;
fsA  = 11;
fsT  = 11;
fsTk = 9;

% ---------- helper ----------
plot_reg_panel = @(ax,x,y,xlab,ylab,ttl) local_regpanel(ax,x,y,xlab,ylab,ttl,...
    dotC,lineC,shadeC,ms,lw,fsA,fsT,fsTk);

% ==========================
% Panel a: X -> M
% ==========================
ax1 = subplot(1,4,1);
plot_reg_panel(ax1, Xz, Mz, ...
    'Beta: Orth-Group best (z)', ...
    'Switch evidence RMSE (z)', ...
    sprintf('Path a: \\beta = %.3f, p = %.3g', a, p_a));

% ==========================
% Panel b: X -> Y
% ==========================
ax2 = subplot(1,4,2);
plot_reg_panel(ax2, Xz, Yz, ...
    'Beta: Orth-Group best (z)', ...
    'Inference task accuracy (z)', ...
    sprintf('Total effect c: \\beta = %.3f, p = %.3g', c_total, p_c));

% ==========================
% Panel c: partial effect of M on Y controlling X
% ==========================
mdl_M = fitlm(Xz, Mz);   M_res = mdl_M.Residuals.Raw;
mdl_Y = fitlm(Xz, Yz);   Y_res = mdl_Y.Residuals.Raw;

ax3 = subplot(1,4,3);
plot_reg_panel(ax3, M_res, Y_res, ...
    'Switch evidence RMSE residuals', ...
    'Accuracy residuals', ...
    sprintf('Path b: \\beta = %.3f, p = %.3g', b_path, p_b_path));

% ==========================
% Panel d: mediation diagram
% ==========================
ax4 = subplot(1,4,4);
axis(ax4,[0 1 0 1]); axis(ax4,'off'); hold(ax4,'on');

% node positions
boxX = [0.05 0.42 0.15 0.10]; %[x y w h]
boxM = [0.45 0.72 0.15 0.10];
boxY = [0.85 0.42 0.15 0.10];

% boxes
rectangle('Position',boxX,'Curvature',0.03,'FaceColor',[0.97 0.97 0.97],...
    'EdgeColor',[0.15 0.15 0.15],'LineWidth',2.0);
rectangle('Position',boxM,'Curvature',0.03,'FaceColor',[0.97 0.97 0.97],...
    'EdgeColor',[0.15 0.15 0.15],'LineWidth',2.0);
rectangle('Position',boxY,'Curvature',0.03,'FaceColor',[0.97 0.97 0.97],...
    'EdgeColor',[0.15 0.15 0.15],'LineWidth',2.0);

text(0.125,0.47,'X','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
text(0.125,0.39,'Neural orthogonality','HorizontalAlignment','center','FontSize',10);

text(0.525,0.77,'M','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
text(0.5,0.86,'Switch evidence RMSE','HorizontalAlignment','center','FontSize',10);

text(0.925,0.47,'Y','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
text(0.9,0.39,'Inference task accuracy','HorizontalAlignment','center','FontSize',10);

text(0.5,1.05,'Summary of mediation model','HorizontalAlignment','center','FontSize',12,'FontWeight','bold');

% significance colors
ca  = ternary(p_a < 0.05, lineC2, nsC2);
cb  = ternary(p_b_path < 0.05, lineC2, nsC2);
cab = ternary((CI_ab(1) > 0 || CI_ab(2) < 0), lineC2, nsC2);

% arrows
annotation(fig,'arrow',[0.80 0.845],[0.53 0.70],'Color',ca,'LineWidth',1.4,'HeadLength',7,'HeadWidth',7);
annotation(fig,'arrow',[0.885 0.930],[0.70 0.53],'Color',cb,'LineWidth',1.4,'HeadLength',7,'HeadWidth',7);
annotation(fig,'arrow',[0.81 0.92],[0.47 0.47],'Color',[0.2 0.2 0.2],'LineWidth',1.4,'HeadLength',7,'HeadWidth',7);

% path text
text(0.32,0.67,sprintf('a = %.3f\np = %.3g',a,p_a),'Rotation',45,...
    'HorizontalAlignment','center','FontSize',10,'Color',ca);

text(0.73,0.67,sprintf('b = %.3f\np = %.3g',b_path,p_b_path),'Rotation',-45,...
    'HorizontalAlignment','center','FontSize',10,'Color',cb);

text(0.50,0.29,sprintf('c = %.3f, p = %.3g\nc'' = %.3f, P = %.3g', ...
    c_total,p_c,c_prime,p_c_prime),...
    'HorizontalAlignment','center','FontSize',10,'Color',[0.15 0.15 0.15]);

text(0.50,0.10,sprintf('Indirect effect: a×b = %.3f\n95%% CI [%.3f, %.3f]', ...
    ab,CI_ab(1),CI_ab(2)),...
    'HorizontalAlignment','center','FontSize',10,'Color',cab);

% title('Mediation model','FontSize',fsT,'FontWeight','normal');

% ---------- panel letters ----------
annotation(fig,'textbox',[0.03 0.90 0.03 0.05],'String','A','EdgeColor','none','FontSize',12,'FontWeight','bold');
annotation(fig,'textbox',[0.27 0.90 0.03 0.05],'String','B','EdgeColor','none','FontSize',12,'FontWeight','bold');
annotation(fig,'textbox',[0.51 0.90 0.03 0.05],'String','C','EdgeColor','none','FontSize',12,'FontWeight','bold');
annotation(fig,'textbox',[0.75 0.90 0.03 0.05],'String','D','EdgeColor','none','FontSize',12,'FontWeight','bold');

% ---------- top summary ----------
% annotation(fig,'textbox',[0.24 0.92 0.52 0.05],...
%     'String',sprintf('Indirect effect a×b = %.3f, 95%% CI [%.3f, %.3f]',ab,CI_ab(1),CI_ab(2)),...
%     'EdgeColor','none','HorizontalAlignment','center','FontSize',10,'Color',[0.1 0.1 0.1]);
set(fig,'InvertHardcopy','off');
ax1.Position = [0.04 0.15 0.22 0.70]; %[left bottom width height];
ax2.Position = [0.28 0.15 0.22 0.70]; %[left bottom width height];
ax3.Position = [0.52 0.15 0.22 0.70]; %[left bottom width height];
ax4.Position = [0.76 0.15 0.20 0.70]; %[left bottom width height];

print(['Geometry_confi_performance'], '-dsvg','-painters');

print(['Geometry_confi_performance.png'], '-dpng','-r300');

%% ---------- local functions ----------
function local_regpanel(ax,x,y,xlab,ylab,ttl,dotC,lineC,shadeC,ms,lw,fsA,fsT,fsTk)
axes(ax); hold on;
scatter(x,y,ms,'MarkerFaceColor',dotC,'MarkerEdgeColor','none','MarkerFaceAlpha',0.7);

mdl = fitlm(x,y);
xg = linspace(min(x),max(x),200)';
[yg,yci] = predict(mdl,xg);

fill([xg; flipud(xg)], [yci(:,1); flipud(yci(:,2))], ...
    shadeC, 'EdgeColor','none', 'FaceAlpha',0.15);
plot(xg,yg,'-','Color',lineC,'LineWidth',lw);

xlabel(xlab,'FontSize',fsA);
ylabel(ylab,'FontSize',fsA);
title(ttl,'FontSize',fsT,'FontWeight','normal');

set(ax,'Box','off','TickDir','out','LineWidth',1.0,'FontSize',fsTk,'Layer','top');
axis square;
end

function out = ternary(cond,a,b)
if cond
    out = a;
else
    out = b;
end
end

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