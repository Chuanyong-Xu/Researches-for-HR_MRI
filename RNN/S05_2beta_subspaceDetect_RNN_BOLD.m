%% ===================== RNN TDR (ROBUST) + Permutation + Paper Figures (A–D) + Angles =====================
% One-file pipeline:
%   Load -> TDR ridge -> orthogonality perm test -> readout alignment perm test
%   -> SW error patterns -> save results -> paper figure panels A–D -> angles figure
%
% Data: *.mat
% infer_data: 1x29 struct array
% Required fields per subject:
%   BOLD [nTrial x nUnit]
%   Nback [1 x nTrial]   (error proxy; >=2 capped to 2)
%   tDev  [1 x nTrial]   (difficulty)
%   Cue [1 x nTrial]
%   RuleChoice [1 x nTrial]
%   SW [nTrial x 1]      (subject switches on next trial; user-defined)
%   mu_switch_estimated [nTrial x 1] (RNN predicted switch probability)

clear; clc; close all;

%% -------------------- Load --------------------
%====================BOLD====================
dataPath = '/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/infer_data_all.mat';
S = load(dataPath);
assert(isfield(S,'infer_data_all'), 'cannot find infer_data');
infer_data = S.infer_data_all;
nSub = numel(infer_data);
%==================== RNN ====================
% % % dataPath = '/Users/vincentxu/Desktop/Model/IOM-CBM-5points-MRI/RNN-29subs/DeepRL_RNN/data/update_data_switch_weight_diff0.mat';
% % % S = load(dataPath);
% % % assert(isfield(S,'infer_data'), 'cannot find infer_data');
% % % infer_data = S.infer_data;
% % % nSub = numel(infer_data);
% % % for subi = 1:length(infer_data)
% % %     infer_data{subi}.BOLD = infer_data{subi}.hidden_state
% % %     infer_data{subi}.mu_switch_estimated = infer_data{subi}.mu_switch_estimated %r[SE]
% % %     infer_data{subi}.pr_of_switch = infer_data{subi}.rnn_pre %s[SE]
% % % end
%==================== RNN ====================

fprintf('Loaded infer_data: nSub=%d\n', nSub);

%% -------------------- Settings --------------------
cap_Nback_at2 = true;          % Nback>=2 -> 2
zscore_hidden = true;          % z-score each unit across trials
zscore_vars   = true;          % z-score predictors (excluding intercept)

% Robust regression for TDR (H ~ X)
tdr_lambda = 1;                % ridge to prevent singular
min_levels = 2;                % predictor must have >=2 unique values to be "usable"

% Readout regression for y ~ hidden
use_ridge_readout = true;
readout_lambda = 10;

% Permutation
nPerm = 10000;                 % increase to 5000-10000 for final
rng(123);                        % reproducible

% SW error threshold (for premature/delayed)
thr = 0.5;
use_prob_only = false;         % true -> only probabilistic metrics (bias,brier), no threshold split

% Diagnostics
print_diag = true;

% Figure export names
outResultsMat = 'RNN_TDR_results.mat';
outFigABCD_png = 'Fig_RNN_TDR_ABCD.png';
outFigAngles_png = 'Fig_Angles_error_diff_rnnpre.png';

%% -------------------- Pre-allocate --------------------
% Geometry / permutation
angle_ED        = nan(nSub,1);   % angle(betaE,betaD), only when both exist
p_ED            = nan(nSub,1);   % perm p for orthogonality (deviation from 90)
valid_ED        = false(nSub,1); % true if both predictors usable

angle_w_plane   = nan(nSub,1);   % angle(w, outcome subspace)
proj_ratio      = nan(nSub,1);   % projection energy ratio of w onto outcome subspace
p_align         = nan(nSub,1);   % perm p for alignment (proj larger than null)

% Proxy performance
acc_rule        = nan(nSub,1);   % mean(RuleChoice==Cue)

% Error patterns
prem_rate       = nan(nSub,1);
delay_rate      = nan(nSub,1);
err_rate        = nan(nSub,1);
bias_prob       = nan(nSub,1);
brier           = nan(nSub,1);

% Save beta and w for angle figures
% nUnit determined from first subject
nUnit = size(infer_data{1}.BOLD, 2);
betaE_all = nan(nSub, nUnit);
betaD_all = nan(nSub, nUnit);
w_all     = nan(nSub, nUnit);
useE      = false(nSub,1);
useD      = false(nSub,1);

bad_constant = 0;

%% ===================== Main loop =====================
for s = 1:nSub
    D = infer_data{s};
    assert(isfield(D,'BOLD') && isfield(D,'Nback') && isfield(D,'tDev') && isfield(D,'mu_switch_estimated') ...
        && isfield(D,'Cue') && isfield(D,'RuleChoice') && isfield(D,'SW'), ...
        'infer_data{%d} 缺少必要字段', s);

    H = double(D.BOLD);          % [trial x unit]
    y = double(D.mu_switch_estimated(:));            % [trial x 1]
    N = double(D.Nback(:));              % [trial x 1]
    T = double(D.tDev(:));               % [trial x 1]
    Cue = double(D.Cue(:));
    RC  = double(D.RuleChoice(:));
    SW  = double(D.SW(:));               % [trial x 1]

    % Select trials

    idx = true(size(y));

    H = H(idx,:); y = y(idx); N = N(idx); T = T(idx);
    Cue = Cue(idx); RC = RC(idx); SW = SW(idx);

    % Clean non-finite
    good = all(isfinite(H),2) & isfinite(y) & isfinite(N) & isfinite(T) & isfinite(Cue) & isfinite(RC) & isfinite(SW);
    H = H(good,:); y = y(good); N = N(good); T = T(good);
    Cue = Cue(good); RC = RC(good); SW = SW(good);

    nTrial = size(H,1);
    if nTrial < 10
        warning('Sub %02d: too few trials (%d). Skip.', s, nTrial);
        continue;
    end

    % Nback cap
    if cap_Nback_at2
        N(N>=2) = 2;
    end

    % Proxy performance
    % acc_rule(s) = mean(RC == Cue); % rule selection accuracy
    acc_rule(s) = sum(double(D.TF(:))==1)/length(double(D.TF(:))); % performance(rule & perception)

    % Z-score hidden
    if zscore_hidden
        H = zscore(H,0,1);
    end

    % Predictor levels
    nLevelN = numel(unique(N));
    nLevelT = numel(unique(T));
    useN = (nLevelN >= min_levels);
    useT = (nLevelT >= min_levels);

    % If both constant -> cannot define outcome subspace
    if ~useN && ~useT
        bad_constant = bad_constant + 1;
        warning('Sub %02d: Nback & tDev both constant (nTrial=%d). Skip sub.', s, nTrial);
        continue;
    end

    % Build design matrix X: intercept + available predictors
    Xraw = [ones(nTrial,1), N, T];
    cols = 1;
    if useN, cols(end+1) = 2; end
    if useT, cols(end+1) = 3; end
    X = Xraw(:, cols); %select the cols for latter use

    if zscore_vars
        for c = 2:size(X,2)
            X(:,c) = zscore(X(:,c));
        end
    end

    if print_diag
        rX = rank(X); %if rX == cols, then the design maxtix is full rank / independent
        cX = cond(X'*X); % if cX < 10, the X-matrix is stable and singular / low collinear
        eigvals = eig(X' * X);
        fprintf('Sub %02d: nTrial=%d | levels N=%d T=%d | useN=%d useT=%d | Xcols=%d rank=%d cond=%.2e\n eig=%d ',...
            s, nTrial, nLevelN, nLevelT, useN, useT, size(X,2), rX, cX, min(eigvals));
    end

    %% ---------- TDR ridge: H ~ X ----------
    pX = size(X,2);
    B = (X' * X + tdr_lambda * eye(pX)) \ (X' * H);   % [pX x nUnit]

    % % % % % % % ---------- TDR OLS: H ~ X OLS via least-squares solver----------
    % % % % % % pX = size(X,2);
    % % % % % % % 
    % % % % % % B = X \ H;                 % [pX x nUnit]

    % Extract beta vectors in unit-space
    betaE = []; betaD = [];
    row = 2;
    if useN
        betaE = B(row,:)'; row = row + 1;
        betaE_all(s,:) = betaE(:)'; useE(s) = true;
    end
    if useT
        betaD = B(row,:)';
        betaD_all(s,:) = betaD(:)'; useD(s) = true;
    end

    % Outcome basis for subspace
    if useN && useT
        Bbasis = [betaE, betaD];
        valid_ED(s) = true;
        angle_ED(s) = angle_vec_vec(betaE, betaD);
    elseif useN
        Bbasis = betaE;
        valid_ED(s) = false;
    else
        Bbasis = betaD;
        valid_ED(s) = false;
    end

    %% ---------- Readout: mu_switch_estimated ~ hidden ----------
    % % % % % % yclip = min(max(y,0),1); % restrict the range of values to 0-1
    % % % % % % y_z = zscore(yclip);
    yclip = y; y_z = zscore(yclip);

    Hdm = [ones(nTrial,1), H];
    pH = size(Hdm,2);
    if use_ridge_readout
        w_full = (Hdm' * Hdm + readout_lambda * eye(pH)) \ (Hdm' * y_z);
    else
        w_full = (Hdm' * Hdm) \ (Hdm' * y_z);
    end
    w = w_full(2:end);
    w_all(s,:) = w(:)';

    % Alignment: w vs outcome subspace
    geom = angle_vec_to_subspace_robust(w, Bbasis);
    angle_w_plane(s) = geom.angle_deg;
    proj_ratio(s)    = geom.proj_energy_ratio;

    %% ---------- Permutation test 1: Orthogonality (only if 2D basis exists) ----------
    if valid_ED(s) && isfinite(angle_ED(s))
        dev_obs = abs(angle_ED(s) - 90);
        null_dev = nan(nPerm,1);

        for ppp = 1:nPerm
            perm = randperm(nTrial);
            Np = N(perm);
            Tp = T(perm);

            Xp = [ones(nTrial,1), Np, Tp];
            if zscore_vars
                Xp(:,2) = zscore(Xp(:,2));
                Xp(:,3) = zscore(Xp(:,3));
            end

            Bp = (Xp' * Xp + tdr_lambda * eye(3)) \ (Xp' * H);
            bE = Bp(2,:)';
            bD = Bp(3,:)';
            null_dev(ppp) = abs(angle_vec_vec(bE,bD) - 90);
        end

        p_ED(s) = (sum(null_dev >= dev_obs) + 1) / (nPerm + 1);
    end

    %% ---------- Permutation test 2: Readout alignment (projection) ----------
    null_proj = nan(nPerm,1);
    for ppp = 1:nPerm
        yperm = y_z(randperm(nTrial));
        if use_ridge_readout
            w_full_p = (Hdm' * Hdm + readout_lambda * eye(pH)) \ (Hdm' * yperm);
        else
            w_full_p = (Hdm' * Hdm) \ (Hdm' * yperm);
        end
        w_p = w_full_p(2:end);
        geom_p = angle_vec_to_subspace_robust(w_p, Bbasis);
        null_proj(ppp) = geom_p.proj_energy_ratio;
    end
    p_align(s) = (sum(null_proj >= proj_ratio(s)) + 1) / (nPerm + 1);  % one-sided: observed proj larger -> not in null

    %% ---------- SW error patterns ----------
    SWb = double(SW(:) > 0);
    brier(s) = mean((yclip - SWb).^2);
    bias_prob(s) = mean(yclip) - mean(SWb);

    if ~use_prob_only
        predSwitch = double(yclip >= thr);
        prem = (predSwitch==1 & SWb==0);
        dely = (predSwitch==0 & SWb==1);

        prem_rate(s)  = mean(prem);
        delay_rate(s) = mean(dely);
        err_rate(s)   = mean(predSwitch ~= SWb);
    end

    fprintf('  -> Sub%02d: angleED=%5.1f pED=%6.3g | proj=%.3f pAlign=%6.3g | prem=%.3f del=%.3f bias=%.3f brier=%.3f | acc=%.3f\n',...
        s, angle_ED(s), p_ED(s), proj_ratio(s), p_align(s), prem_rate(s), delay_rate(s), bias_prob(s), brier(s), acc_rule(s));
end

fprintf('\nSkipped subs (both predictors constant): %d\n', bad_constant);

%% ===================== Save results =====================
results = struct();
results.settings = struct( ...
    'cap_Nback_at2',cap_Nback_at2,'zscore_hidden',zscore_hidden,'zscore_vars',zscore_vars, ...
    'tdr_lambda',tdr_lambda,'readout_lambda',readout_lambda,'nPerm',nPerm,'thr',thr,'use_prob_only',use_prob_only);

results.angle_ED      = angle_ED;
results.p_ED          = p_ED;
results.valid_ED      = valid_ED;

results.angle_w_plane = angle_w_plane;
results.proj_ratio    = proj_ratio;
results.p_align       = p_align;

results.acc_rule      = acc_rule;

results.prem_rate     = prem_rate;
results.delay_rate    = delay_rate;
results.err_rate      = err_rate;
results.bias_prob     = bias_prob;
results.brier         = brier;

results.betaE_all     = betaE_all;
results.betaD_all     = betaD_all;
results.w_all         = w_all;
results.useE          = useE;
results.useD          = useD;

save(outResultsMat, 'results');
fprintf('Saved results to %s\n', outResultsMat);

%% ===================== Paper figures (Panel A–D) + Angles =====================
make_paper_figures_from_results(results, outFigABCD_png, outFigAngles_png);

%% ===================== Done =====================
fprintf('\nDone. Exported figures:\n  %s\n  %s\n  %s\n  %s\n', outFigABCD_png, outFigAngles_png);





%% ===================== Local functions =====================
function make_paper_figures_from_results(results, outABCD_png, outAngles_png)
angle_ED      = results.angle_ED(:);
p_ED          = results.p_ED(:);
valid_ED      = results.valid_ED(:);

angle_w_plane = results.angle_w_plane(:);
proj_ratio    = results.proj_ratio(:);
p_align       = results.p_align(:);

acc_rule      = results.acc_rule(:);
prem_rate     = results.prem_rate(:);
delay_rate    = results.delay_rate(:);

betaE_all = results.betaE_all;
betaD_all = results.betaD_all;
w_all     = results.w_all;
useE      = results.useE(:);
useD      = results.useD(:);

nSub = numel(acc_rule);

% Additional angles: w vs E/D
angle_w_E = nan(nSub,1);
angle_w_D = nan(nSub,1);

for s = 1:nSub
    w = w_all(s,:)';
    if any(~isfinite(w)) || norm(w)<1e-12, continue; end
    if useE(s)
        bE = betaE_all(s,:)';
        if all(isfinite(bE)) && norm(bE)>1e-12
            angle_w_E(s) = angle_vec_vec(w, bE);
        end
    end
    if useD(s)
        bD = betaD_all(s,:)';
        if all(isfinite(bD)) && norm(bD)>1e-12
            angle_w_D(s) = angle_vec_vec(w, bD);
        end
    end
end

% Indices
idxED    = valid_ED & isfinite(angle_ED) & isfinite(p_ED);
idxAlign = isfinite(proj_ratio) & isfinite(p_align) & isfinite(angle_w_plane);
idxGeom  = idxED & isfinite(prem_rate) & isfinite(delay_rate) & isfinite(acc_rule);

% FDR
sig_ED_fdr    = false(nSub,1);
sig_align_fdr = false(nSub,1);
if any(idxED)
    sig_ED_fdr(idxED) = fdr_bh(p_ED(idxED), 0.05);
end
if any(idxAlign)
    sig_align_fdr(idxAlign) = fdr_bh(p_align(idxAlign), 0.05);
end

% -------- print individual-level & group-level stats --------
fprintf('\n=== Individual-level statistics ===\n');
for s = 1:nSub
    if idxED(s)
        fprintf('Sub%02d | angle(E,D)=%.2f deg | p_ED=%.4g | FDR_ED=%d\n', ...
            s, angle_ED(s), p_ED(s), sig_ED_fdr(s));
    end
    if idxAlign(s)
        fprintf('Sub%02d | proj_ratio=%.3f | angle(w,span)=%.2f deg | p_align=%.4g | FDR_align=%d\n', ...
            s, proj_ratio(s), angle_w_plane(s), p_align(s), sig_align_fdr(s));
    end
end

[p_group_ED, chi2_ED]       = fisher_combine_p(p_ED(idxED));
[p_group_align, chi2_align] = fisher_combine_p(p_align(idxAlign));

fprintf('\n=== Group-level statistics ===\n');
fprintf('[Orthogonality] valid N=%d | FDR sig=%d | angle(E,D)=%.2f ± %.2f deg | Fisher chi2=%.2f, p=%.4g\n', ...
    sum(idxED), sum(sig_ED_fdr), mean(angle_ED(idxED),'omitnan'), sem(angle_ED(idxED)), chi2_ED, p_group_ED);

fprintf('[Alignment]     valid N=%d | FDR sig=%d | proj_ratio=%.3f ± %.3f | angle(w,span)=%.2f ± %.2f deg | Fisher chi2=%.2f, p=%.4g\n', ...
    sum(idxAlign), sum(sig_align_fdr), mean(proj_ratio(idxAlign),'omitnan'), sem(proj_ratio(idxAlign)), ...
    mean(angle_w_plane(idxAlign),'omitnan'), sem(angle_w_plane(idxAlign)), chi2_align, p_group_align);

% Style
FS = 11; LW = 1.2;
set(groot, 'defaultAxesFontName','Arial');
set(groot, 'defaultTextFontName','Arial');

%% ---------------- Figure: overlap----------------
fig1 = figure('Color','w','Position',[100 100 1100 750]);
nexttile; hold on

xB = proj_ratio(idxAlign);  yB = angle_w_plane(idxAlign);

yyaxis left
histogram(xB,'BinWidth',0.05,'FaceAlpha',0.45,'FaceColor',[0.2784    0.4471    0.4510])
xlim([0 1.35]); xticks(0:0.2:1); ylabel('Count')
xlabel('Readout alignment to span(error,difficulty)')
title('B  Readout is not in null space')

yyaxis right
scatter(1.12+0.02*randn(size(yB)), yB, 100, 'filled', 'MarkerFaceAlpha',0.45,'MarkerFaceColor',[0.2784    0.4471    0.4510])
ylabel('Angle(w, span) [deg]')
xline(1.05,':')   % 分隔线

text(0.02,0.98,sprintf('n=%d (FDR sig=%d)\nproj=%.2f\\pm%.2f\nangle=%.1f\\pm%.1f', ...
    numel(xB), sum(sig_align_fdr), mean(xB,'omitnan'), sem(xB), mean(yB,'omitnan'), sem(yB)), ...
    'Units','normalized','VerticalAlignment','top','FontSize',FS)

set(gca,'FontSize',FS,'LineWidth',LW); box on


exportgraphics(fig1, outABCD_png, 'Resolution', 300);

%% ---------------- Figure: Angles among error, difficulty, mu_switch_estimated (readout) ----------------
fig2 = figure('Color','w','Position',[150 150 1100 320]);
t2 = tiledlayout(1,4,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
histogram(angle_ED(idxED), 'BinWidth', 3, 'FaceAlpha',0.45,'FaceColor',[0.2784    0.4471    0.4510]);
xline(90,'--','LineWidth',LW);
xlabel('Angle(E,D) [deg]'); ylabel('Count'); title('Angle(E,D)');
set(gca,'FontSize',FS,'LineWidth',LW); box on;

nexttile; hold on;
idxWE = isfinite(angle_w_E);
histogram(angle_w_E(idxWE), 'BinWidth', 3, 'FaceAlpha',0.45,'FaceColor',[0.2784    0.4471    0.4510]);
xlabel('Angle(w,E) [deg]'); ylabel('Count'); title('Angle(readout, error)');
set(gca,'FontSize',FS,'LineWidth',LW); box on;

nexttile; hold on;
idxWD = isfinite(angle_w_D);
histogram(angle_w_D(idxWD), 'BinWidth', 3, 'FaceAlpha',0.45,'FaceColor',[0.2784    0.4471    0.4510]);
xlabel('Angle(w,D) [deg]'); ylabel('Count'); title('Angle(readout, difficulty)');
set(gca,'FontSize',FS,'LineWidth',LW); box on;

nexttile; hold on;
histogram(angle_w_plane(idxAlign), 'BinWidth', 3, 'FaceAlpha',0.45,'FaceColor',[0.2784    0.4471    0.4510]);
xlabel('Angle(w, span) [deg]'); ylabel('Count'); title('Angle(readout, span(E,D))');
set(gca,'FontSize',FS,'LineWidth',LW); box on;

exportgraphics(fig2, outAngles_png, 'Resolution', 300);
end




%% ===================== q + r : FIXED LAYOUT (NO extra axes) =====================
% REQUIRE: infer_data, betaE_all, betaD_all, useE, useD, readout_lambda
% This version DOES NOT create extra axes, so tiledlayout won't break.




plot_RS_geometry_realdata(infer_data, betaE_all, betaD_all, useE, useD, readout_lambda, results);




%% ===== 2D Geometric Schematic for Paper =====
% Schematic of orthogonal representational geometry
% and readout alignment in the Nback–Difficulty subspace

% ---- Define example vectors (near orthogonal) ----
theta_deg = mean(results.angle_ED ); % your empirical mean angle
theta = deg2rad(theta_deg);

vN = [1; 0];                          % Nback axis
vD = [cos(theta); sin(theta)];        % Difficulty axis

% normalize
vN = vN / norm(vN);
vD = vD / norm(vD);

% Example readout (linear combination of both)
w = 0.7*vN + 0.5*vD;
w = w / norm(w);

% ---- Plot ----
figure('Color','w','Position',[300 300 600 600]); hold on;
axis equal
axis([0 1.4 -0.02 1.1])
box off
axis off

% Subspace shading
patch([0 vN(1) vD(1)], ...
      [0 vN(2) vD(2)], ...
      [0.85 0.9 1], ...
      'FaceAlpha',0.2,'EdgeColor','none');

% Plot vectors
quiver(0,0,vN(1),vN(2),0,'LineWidth',2,'Color',[0 0.4 0.9]);
quiver(0,0,vD(1),vD(2),0,'LineWidth',2,'Color',[0.9 0 0]);
quiver(0,0,w(1),w(2),0,'LineWidth',4,'Color',[0 0 0]);

% Labels
text(vN(1)*1.05,vN(2)*1.05,'Number of Errors',...
    'FontSize',12,'Color',[0 0.4 0.9]);

text(vD(1)*1.05,vD(2)*1.05,'Difficulty',...
    'FontSize',12,'Color',[0.9 0 0]);

text(w(1)*1.05,w(2)*1.05,'Readout: Confidence',...
    'FontSize',12,'Color',[0 0 0]);

% Angle arc
ang = linspace(0,theta,100);
r = 0.4;
plot(r*cos(ang), r*sin(ang),'k','LineWidth',1.5);
text(r*cos(theta/2)*0.9, r*sin(theta/2)*1.3,...
    sprintf('\\theta = %.1f^\\circ',theta_deg),...
    'FontSize',12);

xlabel('Neural Dimension 1')
ylabel('Neural Dimension 2')
title('Orthogonal Representational Geometry Supporting Inference')
set(gca,'FontSize',12)












%% ---------------------------------- functions --------------------------
function theta_deg = angle_vec_vec(v1, v2)
v1 = v1(:); v2 = v2(:);
den = norm(v1)*norm(v2);
if den < eps, theta_deg = NaN; return; end
c = dot(v1,v2)/den;
c = max(min(c,1),-1);
theta_deg = acosd(c);
end

function out = angle_vec_to_subspace_robust(v, B)
% Robust: removes NaN/Inf columns and handles rank deficiency via QR.
out.angle_deg = NaN;
out.proj_energy_ratio = NaN;

v = v(:);
if any(~isfinite(v)) || norm(v) < 1e-12
    return;
end
if isempty(B)
    return;
end
B = B(:,:);
goodCol = all(isfinite(B),1) & (vecnorm(B,2,1) > 1e-12);
B = B(:,goodCol);
if isempty(B)
    return;
end

[Q,R] = qr(B,0);
d = abs(diag(R));
if isempty(d)
    return;
end
r = sum(d > 1e-10 * max(d));
Q = Q(:,1:r);
if isempty(Q)
    return;
end

vproj = Q*(Q'*v);
c = norm(vproj)/(norm(v)+eps);
c = max(min(c,1),-1);
out.angle_deg = acosd(c);
out.proj_energy_ratio = (norm(vproj)^2)/(norm(v)^2 + eps);
end

function s = sem(x)
x = x(isfinite(x));
s = std(x)/sqrt(max(1,numel(x)));
end

function sig = fdr_bh(pvals, q)
p = pvals(:);
[ps,ord] = sort(p);
m = numel(p);
thr = (1:m)'/m*q;
k = find(ps<=thr, 1, 'last');
sig = false(size(p));
if ~isempty(k)
    sig(ord(1:k)) = true;
end
sig = reshape(sig, size(pvals));
end

function [p_group, chi2stat] = fisher_combine_p(pvals)
pvals = pvals(isfinite(pvals) & pvals > 0 & pvals <= 1);
if isempty(pvals)
    p_group = NaN;
    chi2stat = NaN;
    return;
end
chi2stat = -2 * sum(log(pvals));
df = 2 * numel(pvals);
p_group = 1 - chi2cdf(chi2stat, df);
end