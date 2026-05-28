%% Combined RNN behavior prediction + training geometry
% Keep SolveCosFunction.m in the same folder / MATLAB path.
% Minimal integration of your original S1 + S5 logic; MDS/PCA not included.
clear; clc; close all
restoredefaultpath

%% paths ------------------------------------------------------------------
data_dir      = '/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/DeepRL_RNN/data';
modelres_path = '/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/model_res_switch_fmri_m_4_29sub.mat';
func_dir      = '/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/MRI_behavioural-29subs/Hiearchica_Reasoning_model/depend/Functions';

addpath(genpath(func_dir));
addpath(pwd);   % for SolveCosFunction.m

final_mat = fullfile(data_dir,'update_data_switch_weight_diff0.mat');
epoch_files = dir(fullfile(data_dir,'update_data_switch_weight_diff_epoch*.mat'));
epoch_id = arrayfun(@(x) sscanf(x.name,'update_data_switch_weight_diff_epoch%d.mat'), epoch_files);
[epoch_id, ord] = sort(epoch_id); epoch_files = epoch_files(ord);

set(0,'DefaultAxesFontName','Helvetica');
set(0,'DefaultTextFontName','Helvetica');
c = {[0.4660 0.6740 0.1880], [0.2784 0.4471 0.451], [0.7412 0.3569 0.0235]};
id = 1:29;

%% behavior data from your original model result --------------------------
load(modelres_path,'res_o')
for iT = 1:3
    res_o_switch_3 = [res_o{iT}.res_switch(:,3), ...
                      res_o{iT}.res_switch(:,2)+res_o{iT}.res_switch(:,4), ...
                      res_o{iT}.res_switch(:,1)+res_o{iT}.res_switch(:,5)];
    res_o_count_3  = [res_o{iT}.res_count(:,3), ...
                      res_o{iT}.res_count(:,2)+res_o{iT}.res_count(:,4), ...
                      res_o{iT}.res_count(:,1)+res_o{iT}.res_count(:,5)];
    Pr_Sw_o_3{iT} = res_o_switch_3 ./ res_o_count_3;
end
for ki = 1:3
    Pr_Sw_o_3{ki} = Pr_Sw_o_3{ki}(1:29,:);
end
Pr_Sw_o_3_new{1} = Pr_Sw_o_3{3};
Pr_Sw_o_3_new{2} = Pr_Sw_o_3{1};
Pr_Sw_o_3_new{3} = Pr_Sw_o_3{2};
Pr_Sw_o_3 = Pr_Sw_o_3_new;

%% choose 3 RNN prediction files for subplot 1-3 --------------------------
% Each subplot uses the same plotting style as your original S1 figure.
if numel(epoch_files) >= 2
    mid_i = max(1, round(numel(epoch_files)/2));
    plot_files = {fullfile(epoch_files(1).folder,epoch_files(1).name), ...
                  fullfile(epoch_files(mid_i).folder,epoch_files(mid_i).name), ...
                  fullfile(epoch_files(end).folder,epoch_files(end).name)};
    plot_titles = {sprintf('Epoch %d',epoch_id(1)), ...
                   sprintf('Epoch %d',epoch_id(mid_i)), ...
                   sprintf('Epoch %d',epoch_id(end))};
else
    plot_files = {final_mat, final_mat, final_mat};
    plot_titles = {'RNN prediction','RNN prediction','RNN prediction'};
end

% R² is saved by Python into each epoch*.mat file; MATLAB only reads it here.
plot_r2 = nan(1,3);
for p = 1:3
    plot_r2(p) = read_r2_from_mat(plot_files{p});
end

%% geometry model grid: keep your original S5 construction ----------------
% NOTE: this intentionally follows your original assignment style:
% beha_d_all(:,iters) = repmat((beha_d_single)',29,1);
iters = 0;
for theta = 0:1*pi/180:pi
for scale = 0:0.1:2
    iters = iters + 1;
    iters_all(iters,1:2) = [theta scale];
    beha_d_single = SolveCosFunction([iters_all(iters,1),iters_all(iters,2)]);
    beha_d_all(:,iters) = repmat((beha_d_single)',29,1);
    subj_c = (repelem(1:29,36)');
end
end

% fixed orthogonal model: angle = 90 deg, scale = 1
fixed_it = find((iters_all(:,1) == 90*pi/180) & (iters_all(:,2) == 1));
if isempty(fixed_it)
    error('Cannot find fixed model: angle = 90 deg, scale = 1');
end
fixed_it = fixed_it(1);

%% fit angle/scale across training epochs --------------------------------
angle_epoch = nan(numel(epoch_files),1);
scale_epoch = nan(numel(epoch_files),1);
beta_fixed_epoch = nan(numel(epoch_files),1);
p_fixed_epoch = nan(numel(epoch_files),1);
train_loss_epoch = nan(numel(epoch_files),1);
val_loss_epoch = nan(numel(epoch_files),1);
test_loss_epoch = nan(numel(epoch_files),1);

% hidden-state perturbation before readout
r2_geo_rm_epoch = nan(numel(epoch_files),1);
test_loss_geo_rm_epoch = nan(numel(epoch_files),1);
delta_loss_geo_rm_epoch = nan(numel(epoch_files),1);

% covariance stretching / anisotropy diagnostics
sigma_mean_epoch    = nan(numel(epoch_files),1);
anisotropy_epoch    = nan(numel(epoch_files),1);
ellipse_epoch       = nan(numel(epoch_files),1);
sigma_mean_sub_all  = nan(numel(epoch_files),29);
lambda1_sub_all     = nan(numel(epoch_files),29);
lambda2_sub_all     = nan(numel(epoch_files),29);
anisotropy_sub_all  = nan(numel(epoch_files),29);
ellipse_sub_all     = nan(numel(epoch_files),29);

for e = 1:numel(epoch_files)
    fprintf('Fitting geometry for epoch %d (%d/%d)\n', epoch_id(e), e, numel(epoch_files));
    epoch_matfile = fullfile(epoch_files(e).folder,epoch_files(e).name);
    conds_9 = collect_conds9_from_RNNmat(epoch_matfile);
    train_loss_epoch(e) = read_scalar_from_mat(epoch_matfile,'train_loss_epoch');
    val_loss_epoch(e)   = read_scalar_from_mat(epoch_matfile,'val_loss_epoch');
    test_loss_epoch(e)  = read_scalar_from_mat(epoch_matfile,'test_loss_epoch');
    r2_geo_rm_epoch(e) = read_scalar_from_mat(epoch_matfile,'r2_geo_rm_epoch');
    test_loss_geo_rm_epoch(e) = read_scalar_from_mat(epoch_matfile,'test_loss_geo_rm_epoch');
    delta_loss_geo_rm_epoch(e) = read_scalar_from_mat(epoch_matfile,'delta_loss_geo_rm_epoch');

    sigma_mean_sub  = nan(29,1);
    lambda1_sub     = nan(29,1);
    lambda2_sub     = nan(29,1);
    anisotropy_sub  = nan(29,1);
    ellipse_sub     = nan(29,1);

    for subs = 1:29
        % exactly the RNN branch in your original S5 script
        mask_data = (reshape(conds_9(subs,:,:), size(conds_9,2), size(conds_9,3)))';
        mask_data = zscore(mask_data);

        Sigma_sub = covdiag((mask_data'));
        Sigma_sub = (Sigma_sub + Sigma_sub') ./ 2;

        eigvals = sort(eig(Sigma_sub), 'descend');
        eigvals(eigvals < 0) = 0;

        sigma_mean_sub(subs) = mean(diag(Sigma_sub), 'omitnan');
        lambda1_sub(subs)    = eigvals(1);
        lambda2_sub(subs)    = eigvals(2);
        anisotropy_sub(subs) = eigvals(1) ./ max(eigvals(2), eps);
        ellipse_sub(subs)    = (eigvals(1) - eigvals(2)) ./ max(eigvals(1) + eigvals(2), eps);

        neur_cond_RDM = pdist((mask_data'), 'mahal', Sigma_sub);
        neur_RDM(subs,:) = (neur_cond_RDM');
    end

    sigma_mean_epoch(e)       = mean(sigma_mean_sub, 'omitnan');
    anisotropy_epoch(e)       = mean(anisotropy_sub, 'omitnan');
    ellipse_epoch(e)          = mean(ellipse_sub, 'omitnan');
    sigma_mean_sub_all(e,:)   = sigma_mean_sub;
    lambda1_sub_all(e,:)      = lambda1_sub;
    lambda2_sub_all(e,:)      = lambda2_sub;
    anisotropy_sub_all(e,:)   = anisotropy_sub;
    ellipse_sub_all(e,:)      = ellipse_sub;

    boot_n = 1;
    for b = 1:boot_n
        bootSubIdx = randsample(29, 29, true);
        neur_RDM_b = neur_RDM(bootSubIdx,:);
        neur_RDM_b = reshape(neur_RDM_b.', [], 1);
        if b == boot_n
            neur_RDM_b = reshape(neur_RDM.', [], 1);
        end
        neur_cond_RDM = neur_RDM_b;

        LogLikelihood = nan(size(iters_all,1),1);
        pvalue = nan(length(iters_all),1);

        parfor it = 1:size(iters_all,1)
            beha_d = beha_d_all(:,it);
            VarNames = {'beha_d','neur_cond_RDM','subj_c'};
            tbl = array2table([beha_d,neur_cond_RDM,subj_c],'VariableNames',VarNames);
            formula = 'neur_cond_RDM ~ 1 + beha_d + (1 + beha_d | subj_c)';
            fit_results = fitlme(tbl,formula);
            LogLikelihood(it) = -fit_results.ModelCriterion.LogLikelihood;
            pvalue(it)        = fit_results.Coefficients.pValue(2);
        end

        min_RDMpredict_log{1} = find(min((LogLikelihood))==(LogLikelihood));
        min_RDMpredict_p(1)   = min(pvalue(min_RDMpredict_log{1}));
        angle_rois(1) = mean(iters_all([min_RDMpredict_log{1}],1));
        scale_rois(1) = mean(iters_all([min_RDMpredict_log{1}],2));

        % beta for the fixed orthogonal model: angle = 90 deg, scale = 1
        beha_d_fixed = beha_d_all(:,fixed_it);
        VarNames = {'beha_d','neur_cond_RDM','subj_c'};
        tbl_fixed = array2table([beha_d_fixed,neur_cond_RDM,subj_c],'VariableNames',VarNames);
        formula = 'neur_cond_RDM ~ 1 + beha_d + (1 + beha_d | subj_c)';
        fit_fixed = fitlme(tbl_fixed,formula);
        beta_fixed_epoch(e) = fit_fixed.Coefficients.Estimate(2);
        p_fixed_epoch(e)    = fit_fixed.Coefficients.pValue(2);
    end

    angle_epoch(e) = angle_rois(1);
    scale_epoch(e) = scale_rois(1);
end




%% plot -------------------------------------------------------------------

% unified soft palette, matched to the reference figure
colBlue   = [147 203 237] ./ 255;   % #93CBED
colPink   = [247 193 193] ./ 255;   % #F7C1C1
colOrange = [249 199 128] ./ 255;   % #F9C780
colTeal   = [174 219 210] ./ 255;   % #AEDBD2
colPurple = [188 156 200] ./ 255;   % #BC9CC8
colGreen  = [213 233 198] ./ 255;   % #D5E9C6
colPale   = [228 240 247] ./ 255;   % #E4F0F7
colAxis   = [0.22 0.22 0.22];
colZero   = [0.70 0.70 0.70];

% aligned manual layout: A-C / D-F / G-I / J-K
pos = [0.070 0.800 0.250 0.140;
       0.385 0.800 0.250 0.140;
       0.700 0.800 0.250 0.140;
       0.070 0.565 0.230 0.155;
       0.395 0.565 0.250 0.155;
       0.700 0.565 0.250 0.155;
       0.070 0.330 0.250 0.145;
       0.385 0.330 0.250 0.145;
       0.700 0.330 0.250 0.145;
       0.070 0.050 0.565 0.195;
       0.700 0.050 0.250 0.195];

text_t = {'A','B','C','D','E'};
h_main = figure('Color','w','unit','pixels','position',[100 50 900 1050]);

% A-C: keep original plotting logic; only harmonize axes/border style
for p = 1:3
    subplot(4,3,p,'Position',pos(p,:)); hold on
    [~, Pr_Sw_inp_3_rnn] = collect_behavior_from_RNNmat(plot_files{p});

    for ki = 1:3
        ste = nanstd(Pr_Sw_o_3{ki})./sqrt((length(id)-sum(isnan(Pr_Sw_o_3{ki}))));
        be(ki)=errorbar(1:3, nanmean(Pr_Sw_o_3{ki}), ste, 'ko', ...
            'LineWidth',2,'MarkerFaceColor',c{ki},'MarkerSize',10);
    end

    rn(1)=fillsteplottusc(Pr_Sw_inp_3_rnn{2},3); hold on
    rn(2)=fillsteplottust(Pr_Sw_inp_3_rnn{3},3); hold on
    rn(3)=fillsteplotgreen(Pr_Sw_inp_3_rnn{1},3); hold on

    set(gca,'XTiCk',1:1:3);
    set(gca,'XTickLabel',{'Hard (D1)' 'Middle (D2)' 'Easy (D3)'});
    axis([0.5 3.5 0 1]); yticks([0:.25:1])
    ylabel('Pr(Switch)','fontsize',12); xlabel('Difficulties','fontsize',12)
    if isnan(plot_r2(p))
        title(plot_titles{p},'fontsize',13,'FontWeight','Normal')
    else
        title(sprintf('%s, R^2 = %.3f', plot_titles{p}, plot_r2(p)), ...
            'fontsize',13,'FontWeight','Normal')
    end
    if p == 3
        lgd = legend([be(3),be(2),be(1),rn(2),rn(1),rn(3)], ...
            {'E2','E1','E0',...
            'E2 RNN','E1 RNN','E0 RNN'}, ...
            'Orientation','horizontal'); legend boxoff;
        lgd.Units = 'normalized';
        lgd.Position = [pos(2,1)+0.08, pos(3,2)+0.18, 0.10, 0.01];
    end
    text(-0.12, 1.18, text_t{p},'Units','normalized', ...
        'FontSize',18, 'FontWeight','bold', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top');
    format_axis(gca,colAxis);
end

% D: angle and scale
subplot(4,3,4,'Position',pos(4,:)); hold on
if ~isempty(epoch_files)
    yyaxis left
    plot(epoch_id, angle_epoch.*180/pi, 'o-', 'LineWidth',2.2, ...
        'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colBlue);
    ylabel('Angle (degree)'); ylim([85 105])

    yyaxis right
    plot(epoch_id, scale_epoch, 's--', 'LineWidth',2.2, ...
        'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colPurple);
    ylabel('Scale');

    xlabel('Epoch')
    title('Geometry across training','fontsize',14,'FontWeight','Normal')
    legend({'Angle','Scale'},'Location','best','Box','off','FontSize',11)
    ax = gca;
    ax.YAxis(1).Color = colBlue;
    ax.YAxis(2).Color = colPurple;
    ax.XColor = colAxis;
    ax.LineWidth = 1.5;
    ax.TickDir = 'out';
    ax.Box = 'off';
else
    text(.5,.5,'No epoch*.mat files found','HorizontalAlignment','center'); axis off
end
text(-0.12, 1.18, 'D','Units','normalized', ...
    'FontSize',18, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');

% E: fixed orthogonal model under three distance metrics
subplot(4,3,5,'Position',pos(5,:)); hold on
hE(1)=plot(epoch_id, beta_fixed_epoch, 'o-', 'LineWidth',2.2, ...
    'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colBlue);
yline(0,'--','LineWidth',1.1,'Color',colZero);
xlabel('Epoch');
ylabel('\beta for model: [90^\circ, 1]');
legend(hE, {'Mahalanobis'}, ...
    'Location','East','Orientation','vertical'); legend boxoff;
title('Different distance metrics','fontsize',14,'FontWeight','Normal');
text(-0.12, 1.18, 'E','Units','normalized', ...
    'FontSize',18, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');
format_axis(gca,colAxis);

% F: MSE loss across training
subplot(4,3,6,'Position',pos(6,:)); hold on
hLoss = gobjects(0); lossNames = {};
% if any(~isnan(train_loss_epoch))
%     hLoss(end+1) = plot(epoch_id, train_loss_epoch, 'o-', 'LineWidth',2.2, ...
%         'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colBlue);
%     lossNames{end+1} = 'Train';
% end
% if any(~isnan(val_loss_epoch))
%     hLoss(end+1) = plot(epoch_id, val_loss_epoch, 's-', 'LineWidth',2.2, ...
%         'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colTeal);
%     lossNames{end+1} = 'Val';
% end
if any(~isnan(test_loss_epoch))
    hLoss(end+1) = plot(epoch_id, test_loss_epoch, '^-', 'LineWidth',2.2, ...
        'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colPurple);
    lossNames{end+1} = 'Test';
end
ylabel('MSE loss'); xlabel('Epoch');
title('Prediction loss','fontsize',14,'FontWeight','Normal');
if ~isempty(hLoss), legend(hLoss, lossNames, 'Location','best','Box','off'); end
text(-0.12, 1.18, 'F','Units','normalized', ...
    'FontSize',18, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');
format_axis(gca,colAxis);

% covariance diagnostics
anis_mean = mean(anisotropy_sub_all, 2, 'omitnan');
anis_sem  = std(anisotropy_sub_all, 0, 2, 'omitnan') ./ sqrt(size(anisotropy_sub_all,2));
elli_mean = mean(ellipse_sub_all, 2, 'omitnan');
elli_sem  = std(ellipse_sub_all, 0, 2, 'omitnan') ./ sqrt(size(ellipse_sub_all,2));
sigma_mean = mean(sigma_mean_sub_all, 2, 'omitnan');
sigma_sem  = std(sigma_mean_sub_all, 0, 2, 'omitnan') ./ sqrt(size(sigma_mean_sub_all,2));

subplot(4,3,7,'Position',pos(7,:)); hold on
fill([epoch_id(:); flipud(epoch_id(:))], ...
     [sigma_mean - sigma_sem; flipud(sigma_mean + sigma_sem)], ...
     colBlue, 'EdgeColor','none', 'FaceAlpha',0.22);
plot(epoch_id, sigma_mean, 'o-', 'LineWidth',2.2, ...
    'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colBlue);
xlabel('Epoch'); ylabel('Diagonal covariance');
title('Covariance magnitude','fontsize',14,'FontWeight','Normal');
text(-0.12, 1.28, 'G','Units','normalized', ...
    'FontSize',18, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');
format_axis(gca,colAxis);

subplot(4,3,8,'Position',pos(8,:)); hold on
fill([epoch_id(:); flipud(epoch_id(:))], ...
     [anis_mean - anis_sem; flipud(anis_mean + anis_sem)], ...
     colTeal, 'EdgeColor','none', 'FaceAlpha',0.25);
plot(epoch_id, anis_mean, 'o-', 'LineWidth',2.2, ...
    'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colTeal);
xlabel('Epoch'); ylabel('\lambda_1 / \lambda_2');
title('Covariance anisotropy','fontsize',14,'FontWeight','Normal');
text(-0.12, 1.28, 'H','Units','normalized', ...
    'FontSize',18, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');
format_axis(gca,colAxis);

subplot(4,3,9,'Position',pos(9,:)); hold on
fill([epoch_id(:); flipud(epoch_id(:))], ...
     [elli_mean - elli_sem; flipud(elli_mean + elli_sem)], ...
     colPurple, 'EdgeColor','none', 'FaceAlpha',0.24);
plot(epoch_id, elli_mean, 'o-', 'LineWidth',2.2, ...
    'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colPurple);
xlabel('Epoch'); ylabel('(\lambda_1 - \lambda_2) / (\lambda_1 + \lambda_2)');
title('Ellipse index','fontsize',14,'FontWeight','Normal');
text(-0.12, 1.28, 'I','Units','normalized', ...
    'FontSize',18, 'FontWeight','bold', ...
    'HorizontalAlignment','left', 'VerticalAlignment','top');
format_axis(gca,colAxis);

% hidden-state geometry perturbation before readout
if any(~isnan(test_loss_geo_rm_epoch)) || any(~isnan(delta_loss_geo_rm_epoch))
    subplot(4,3,10:11,'Position',pos(10,:)); hold on
    hP = gobjects(0); legP = {};
    if any(~isnan(test_loss_epoch))
        hP(end+1) = plot(epoch_id, test_loss_epoch, 'o-', 'LineWidth',2.2, ...
            'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colBlue);
        legP{end+1} = 'Original hidden';
    end
    if any(~isnan(test_loss_geo_rm_epoch))
        hP(end+1) = plot(epoch_id, test_loss_geo_rm_epoch, 's-', 'LineWidth',2.2, ...
            'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colPink);
        legP{end+1} = 'Geometry-removed hidden';
    end
    xlabel('Epoch'); ylabel('MSE loss');
    title('Prediction after hidden-state perturbation','FontSize',14,'FontWeight','Normal');
    if ~isempty(hP), legend(hP, legP, 'Location','best','Box','off'); end
    text(-0.07, 1.22, 'J','Units','normalized', ...
        'FontSize',18, 'FontWeight','bold', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top');
    format_axis(gca,colAxis);

    subplot(4,3,12,'Position',pos(11,:)); hold on
    plot(epoch_id, delta_loss_geo_rm_epoch, 'o-', 'LineWidth',2.2, ...
        'MarkerFaceColor','w', 'MarkerSize',5.5, 'Color',colOrange);
    yline(0,'--','LineWidth',1.1,'Color',colZero);
    xlabel('Epoch'); ylabel('\Delta MSE loss');
    title('Geometry contribution','FontSize',14,'FontWeight','Normal');
    text(-0.12, 1.22, 'K','Units','normalized', ...
        'FontSize',18, 'FontWeight','bold', ...
        'HorizontalAlignment','left', 'VerticalAlignment','top');
    format_axis(gca,colAxis);
end

print('RNN_prediction_training_geometry_combined_style.svg','-dsvg','-painters');

save('RNN_training_angle_scale.mat', ...
     'epoch_id','angle_epoch','scale_epoch', ...
     'beta_fixed_epoch','p_fixed_epoch', ...
     'train_loss_epoch','val_loss_epoch','test_loss_epoch', ...
     'r2_geo_rm_epoch','test_loss_geo_rm_epoch','delta_loss_geo_rm_epoch', ...
     'sigma_mean_epoch','anisotropy_epoch','ellipse_epoch', ...
     'sigma_mean_sub_all','lambda1_sub_all','lambda2_sub_all', ...
     'anisotropy_sub_all','ellipse_sub_all', ...
     'plot_r2')







%% local functions --------------------------------------------------------

function format_axis(ax, colAxis)
set(ax, 'FontSize',12, ...
        'LineWidth',1.5, ...
        'TickDir','out', ...
        'Box','off', ...
        'XColor',colAxis, ...
        'YColor',colAxis);
end

function r2_epoch = read_r2_from_mat(matfile)
r2_epoch = nan;
info = whos('-file', matfile);
names = {info.name};
if any(strcmp(names,'r2_epoch'))
    S = load(matfile,'r2_epoch');
    r2_epoch = double(S.r2_epoch);
    r2_epoch = r2_epoch(1);
end
end

function val = read_scalar_from_mat(matfile,varname)
val = nan;
info = whos('-file', matfile);
names = {info.name};
if any(strcmp(names,varname))
    S = load(matfile,varname);
    val = double(S.(varname));
    val = val(1);
end
end

function [Pr_Sw_o_3_baye, Pr_Sw_inp_3_rnn] = collect_behavior_from_RNNmat(matfile)
load(matfile,'infer_data')
for subs = 1:length(infer_data)
    row = getrow(infer_data,subs);
    tDev = double(row.tDev(:))';
    Nback = double(row.Nback(:))';
    pr_of_switch = double(row.pr_of_switch(:))';
    rnn_pre = double(row.rnn_pre(:))';
    for con = 1:3
        Pr_Sw_o_3_baye{con}(subs,1) = mean(pr_of_switch(tDev==1 & Nback==con-1));
        Pr_Sw_o_3_baye{con}(subs,2) = mean(pr_of_switch(tDev==2 & Nback==con-1));
        Pr_Sw_o_3_baye{con}(subs,3) = mean(pr_of_switch(tDev==3 & Nback==con-1));

        Pr_Sw_inp_3_rnn{con}(subs,1) = mean(rnn_pre(tDev==1 & Nback==con-1));
        Pr_Sw_inp_3_rnn{con}(subs,2) = mean(rnn_pre(tDev==2 & Nback==con-1));
        Pr_Sw_inp_3_rnn{con}(subs,3) = mean(rnn_pre(tDev==3 & Nback==con-1));
    end
end
for ki = 1:3
    Pr_Sw_o_3_baye{ki}  = Pr_Sw_o_3_baye{ki}(1:29,:);
    Pr_Sw_inp_3_rnn{ki} = Pr_Sw_inp_3_rnn{ki}(1:29,:);
end
end

function conds_9 = collect_conds9_from_RNNmat(matfile)
load(matfile,'infer_data')
for subs = 1:length(infer_data)
    row = getrow(infer_data,subs);
    tDev = double(row.tDev(:))';
    Nback = double(row.Nback(:))';
    hidden_state = double(row.hidden_state);
    if size(hidden_state,1) ~= length(tDev)
        hidden_state = hidden_state';
    end

    % exactly the condition order from your original S1 script
    conds_9(subs,1,:) = mean(hidden_state(tDev==1 & Nback==0,:));
    conds_9(subs,2,:) = mean(hidden_state(tDev==2 & Nback==0,:));
    conds_9(subs,3,:) = mean(hidden_state(tDev==3 & Nback==0,:));
    conds_9(subs,4,:) = mean(hidden_state(tDev==1 & Nback==1,:));
    conds_9(subs,5,:) = mean(hidden_state(tDev==2 & Nback==1,:));
    conds_9(subs,6,:) = mean(hidden_state(tDev==3 & Nback==1,:));
    conds_9(subs,7,:) = mean(hidden_state(tDev==1 & Nback>=2,:));
    conds_9(subs,8,:) = mean(hidden_state(tDev==2 & Nback>=2,:));
    conds_9(subs,9,:) = mean(hidden_state(tDev==3 & Nback>=2,:));
end
end

function row = getrow(infer_data,subs)
if iscell(infer_data)
    row = infer_data{subs};
else
    row = infer_data(subs);
end
end

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
