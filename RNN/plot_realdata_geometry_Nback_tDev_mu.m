function plot_realdata_geometry_Nback_tDev_mu()
% Real-data 2D geometry plot:
%   x = Nback (error count)
%   y = tDev  (difficulty)
%   color = mu_switch_estimated (confidence)
%
% Input: a .mat file that contains variable "infer_data" (1 x nSub struct array)
% Output: figures saved to ./FIG_REAL_GEOM/

clc; close all;

%% ===================== USER PATHS =====================
% >>> Set this to your uploaded file path (local or HPC path).
% If you are running inside the same environment as your previous code, use:
infer_mat_path = '/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/DeepRL_RNN/data/update_data_switch_weight_diff0.mat';


out_dir = fullfile(pwd, 'FIG_REAL_GEOM');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% ===================== LOAD =====================
S = load(infer_mat_path);
if ~isfield(S, 'infer_data')
    error('Cannot find variable "infer_data" in %s', infer_mat_path);
end
infer_data = S.infer_data;
nSub = numel(infer_data);

%% ===================== EXTRACT TRIAL-LEVEL DATA =====================
subj = [];
Nback = [];
tDev  = [];
mu    = [];

for s = 1:nSub
    % Each infer_data(s) is a struct with fields; values are typically 1xT.
    nb = double(infer_data{s}.Nback(:));
    nb(nb>=2)=2;
    td = double(infer_data{s}.tDev(:));
    m  = double(infer_data{s}.mu_switch_estimated(:));

    % basic sanity
    T = numel(nb);
    if numel(td) ~= T || numel(m) ~= T
        error('Subject %d: length mismatch among Nback/tDev/mu.', s);
    end

    % remove NaNs if any
    ok = ~(isnan(nb) | isnan(td) | isnan(m));
    nb = nb(ok); td = td(ok); m = m(ok);

    subj = [subj; s * ones(numel(nb),1)];
    Nback = [Nback; nb];
    tDev  = [tDev;  td];
    mu    = [mu;    m];
end

% Ranges (discrete levels)
uN = unique(Nback(:))';
uD = unique(tDev(:))';

fprintf('Loaded %d subjects, pooled trials = %d\n', nSub, numel(mu));
fprintf('Nback levels: %s\n', mat2str(uN));
fprintf('tDev levels : %s\n', mat2str(uD));

%% ===================== CONDITION GRID: mean mu per (Nback,tDev) =====================
% grid means and counts
meanMu = nan(numel(uD), numel(uN));   % rows: tDev, cols: Nback
nCell  = zeros(numel(uD), numel(uN));

for iD = 1:numel(uD)
    for iN = 1:numel(uN)
        idx = (tDev==uD(iD)) & (Nback==uN(iN));
        nCell(iD,iN) = sum(idx);
        if nCell(iD,iN) > 0
            meanMu(iD,iN) = mean(mu(idx));
        end
    end
end

%% ===================== FIT TWO MODELS (REAL DATA) =====================
% Additive: mu ~ 1 + Nback + tDev
% Multiplicative/interaction: mu ~ 1 + Nback + tDev + Nback*tDev
tbl = table(mu, Nback, tDev);

mdl_add  = fitlm(tbl, 'mu ~ 1 + Nback + tDev');
mdl_mult = fitlm(tbl, 'mu ~ 1 + Nback + tDev + Nback:tDev');

fprintf('\n=== Fit summary (pooled trials) ===\n');
fprintf('Additive:  R2=%.4f, adjR2=%.4f\n', mdl_add.Rsquared.Ordinary, mdl_add.Rsquared.Adjusted);
fprintf('Mult(Int): R2=%.4f, adjR2=%.4f\n', mdl_mult.Rsquared.Ordinary, mdl_mult.Rsquared.Adjusted);

%% ===================== MAKE FIGURE =====================
fig = figure('Color','w','Position',[100 100 1350 420]);

% ---------- Panel A: trial-level scatter (jittered tDev) ----------
subplot(1,3,1); hold on;

% jitter tDev to avoid overplot (since tDev is discrete)
rng(1);
jit = (rand(size(tDev))-0.5)*0.18;
yJ = tDev + jit;

sc = scatter(Nback, yJ, 14, mu, 'filled');
set(gca,'YDir','normal');
xlabel('Nback (error count)');
ylabel('tDev (difficulty)');
title('Real data (trial-level)');
cb = colorbar; ylabel(cb,'\mu_{switch\_estimated} (confidence)');
grid on; box on;
set(gca,'FontSize',11);

% show discrete lines for tDev levels
for iD = 1:numel(uD)
    yline(uD(iD), ':', 'LineWidth', 1);
end
ylim([min(uD)-0.6, max(uD)+0.6]);

% ---------- Panel B: condition grid heatmap (mean mu) ----------
subplot(1,3,2); hold on;
imagesc(uN, uD, meanMu);
set(gca,'YDir','normal');
xlabel('Nback (error count)');
ylabel('tDev (difficulty)');
title('Condition means: E[N\mu | Nback, tDev]');
cb = colorbar; ylabel(cb,'Mean confidence');
axis tight; box on;
set(gca,'FontSize',11);

% annotate each cell with mean and n
for iD = 1:numel(uD)
    for iN = 1:numel(uN)
        if ~isnan(meanMu(iD,iN))
            text(uN(iN), uD(iD), sprintf('%.2f\nn=%d', meanMu(iD,iN), nCell(iD,iN)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',9, 'Color','k');
        else
            text(uN(iN), uD(iD), sprintf('n=%d', nCell(iD,iN)), ...
                'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
                'FontSize',9, 'Color',[0.3 0.3 0.3]);
        end
    end
end

% ---------- Panel C: model contours (additive vs multiplicative) ----------
subplot(1,3,3); hold on;

% create continuous grid for smooth contours
xg = linspace(min(uN), max(uN), 120);
yg = linspace(min(uD), max(uD), 120);
[XG,YG] = meshgrid(xg,yg);

% predictions
mu_add_pred  = predict(mdl_add,  table(XG(:), YG(:), 'VariableNames',{'Nback','tDev'}));
mu_mult_pred = predict(mdl_mult, table(XG(:), YG(:), 'VariableNames',{'Nback','tDev'}));
MU_ADD  = reshape(mu_add_pred,  size(XG));
MU_MULT = reshape(mu_mult_pred, size(XG));

% background: observed condition means as points
% (plot means at discrete grid points)
for iD = 1:numel(uD)
    for iN = 1:numel(uN)
        if ~isnan(meanMu(iD,iN))
            plot(uN(iN), uD(iD), 'ko', 'MarkerSize', 5, 'MarkerFaceColor','w');
        end
    end
end

% contour lines
nLevels = 8;
contour(XG, YG, MU_ADD,  nLevels, 'LineWidth', 1.4);               % additive
contour(XG, YG, MU_MULT, nLevels, '--', 'LineWidth', 1.4);          % multiplicative

xlabel('Nback (error count)');
ylabel('tDev (difficulty)');
title(sprintf('Model contours (pooled)\nAdd: adjR^2=%.3f | Mult: adjR^2=%.3f', ...
    mdl_add.Rsquared.Adjusted, mdl_mult.Rsquared.Adjusted));

legend({'Observed condition means','Additive fit (solid)','Multiplicative fit (dashed)'}, ...
    'Location','best', 'FontSize',9);

grid on; box on;
set(gca,'FontSize',11);
axis([min(uN) max(uN) min(uD) max(uD)]);

sgtitle('Real-data 2D geometry: Nback × tDev → confidence (\mu_{switch\_estimated})', 'FontSize', 14);

%% ===================== SAVE =====================
out_png = fullfile(out_dir, 'REAL_geom_Nback_tDev_mu.png');
out_pdf = fullfile(out_dir, 'REAL_geom_Nback_tDev_mu.pdf');

exportgraphics(fig, out_png, 'Resolution', 300);
exportgraphics(fig, out_pdf, 'ContentType','vector');

fprintf('\nDONE. Figures saved to:\n  %s\n  %s\n', out_png, out_pdf);

end