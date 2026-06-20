function results = run_Nback_tDev_ortho_perm(matPath, nPerm)
% RUN_NBACK_TDEV_ORTHO_PERM
% 目标：
%   两个维度：
%     1) 错误次数：Nback（>=2 归为 2）
%     2) 难度：tDev
%   从 hidden_state (trial × units) 中用多元回归提取两条表征向量：
%     vN: Nback 维度向量
%     vD: tDev 维度向量
%   然后检验：
%     (A) vN 与 vD 是否正交（cos/angle）+ trial-level permutation
%     (B) readout 向量 w 是否落在 span(vN,vD) 子空间内（alignment）+ y-permutation
%   并做：
%     (C) group-level sign-flip
%     (D) performance × orthogonality（performance=accuracy，用 RuleChoice vs Cue 仅用于计算表现，不进入维度定义）
%     (E) group plots
%
% 用法：
%   results = run_Nback_tDev_ortho_perm('/mnt/data/update_data_switch_weight_diff0.mat', 5000);

if nargin < 2 || isempty(nPerm), nPerm = 5000; end

S = load(matPath);
infer_data = S.infer_data;
nSub = numel(infer_data);

subj(nSub) = struct();
fprintf('Loaded %d subjects from %s\n', nSub, matPath);

for s = 1:nSub
    D = infer_data{s};

    % --- trial variables ---
    Nback = double(D.Nback(:));        % T×1
    tDev  = double(D.tDev(:));         % T×1
    H     = double(D.hidden_state);    % T×units

    % --- Nback>=2 -> 2 ---
    NbackCat = Nback;
    NbackCat(NbackCat >= 2) = 2;       % {0,1,2}

    % --- standardize predictors ---
    xN = zscore(NbackCat);
    xD = zscore(tDev);

    % ============ 1) 提取 vN, vD（两维度回归）============
    % H = b0 + bN*xN + bD*xD + noise
    X = [ones(size(H,1),1), xN, xD];   % T×3
    B = X \ H;                         % 3×units

    vN = B(2,:)';   % units×1
    vD = B(3,:)';   % units×1

    % normalize to unit length
    vN = vN / norm(vN);
    vD = vD / norm(vD);

    cosND_obs = dot(vN, vD);
    angND_obs = acosd(max(-1,min(1,cosND_obs)));

    % ============ 2) readout 向量 w 与 span(vN,vD) 的对齐 ============
    % 默认用 pr_of_switch 当 readout，你也可改成 mu_switch_estimated
    y = double(D.pr_of_switch(:));
    y = zscore(y);

    % y = a0 + H*w + noise
    Z = [ones(size(H,1),1), zscore(H,0,1)]; % T×(1+units)
    w = (Z \ y);
    w = w(2:end);
    w = w / norm(w);

    Und = orth([vN vD]); % units×k (k<=2)
    proj_w = Und*(Und'*w);
    align_obs = (norm(proj_w)^2) / (norm(w)^2); % [0,1]

    % ============ 3) 置换检验 ============
    cosND_perm = zeros(nPerm,1);
    align_perm = zeros(nPerm,1);

    Tn = size(H,1);
    for p = 1:nPerm
        % ---- A: 打乱 (xN,xD) 与 trial 的对应（检验“几何关系是否只是偶然”）----
        idx = randperm(Tn);
        Xp = [ones(Tn,1), xN(idx), xD(idx)];
        Bp = Xp \ H;

        vNp = Bp(2,:)'; vDp = Bp(3,:)';
        vNp = vNp / norm(vNp);
        vDp = vDp / norm(vDp);

        cosND_perm(p) = dot(vNp, vDp);

        % ---- B: 打乱 readout y（检验 alignment 是否超出机会）----
        idy = randperm(Tn);
        yp = y(idy);

        wp = (Z \ yp);
        wp = wp(2:end);
        wp = wp / norm(wp);

        proj_wp = Und*(Und'*wp);
        align_perm(p) = (norm(proj_wp)^2) / (norm(wp)^2);
    end

    % p-values
    % 非正交性：|cos| 越大越“不正交”
    p_nonorth = (sum(abs(cosND_perm) >= abs(cosND_obs)) + 1) / (nPerm + 1);

    % 更正交：|cos| 越小越“更正交”
    p_orth = (sum(abs(cosND_perm) <= abs(cosND_obs)) + 1) / (nPerm + 1);

    % alignment：越大越“在子空间里”
    p_align = (sum(align_perm >= align_obs) + 1) / (nPerm + 1);

    % store
    subj(s).ID         = D.ID;
    subj(s).cosND_obs  = cosND_obs;
    subj(s).angND_obs  = angND_obs;
    subj(s).p_nonorth  = p_nonorth;
    subj(s).p_orth     = p_orth;

    subj(s).align_obs  = align_obs;
    subj(s).p_align    = p_align;

    subj(s).cosND_perm = cosND_perm;
    subj(s).align_perm = align_perm;

    fprintf('Sub %02d: cosND=%.3f (angle=%.1f°), p_nonorth=%.4f, p_orth=%.4f, align=%.3f, p_align=%.4f\n', ...
        s, cosND_obs, angND_obs, p_nonorth, p_orth, align_obs, p_align);
end

% ============ 4) 组水平 sign-flip ============
cos_vec   = [subj.cosND_obs]';
align_vec = [subj.align_obs]';

nPermG = 20000;

% (i) cosND 是否偏离 0（双侧）
mu_cos_obs = mean(cos_vec);
mu_cos_perm = zeros(nPermG,1);
for p = 1:nPermG
    signs = (rand(nSub,1) > 0.5)*2 - 1;
    mu_cos_perm(p) = mean(cos_vec .* signs);
end
p_group_cos_twoside = (sum(abs(mu_cos_perm) >= abs(mu_cos_obs)) + 1) / (nPermG + 1);

% (ii) align 是否高于各自置换基线（更稳健）
align_base = zeros(nSub,1);
for s = 1:nSub
    align_base(s) = mean(subj(s).align_perm);
end
align_delta = align_vec - align_base;

mu_align_obs = mean(align_delta);
mu_align_perm = zeros(nPermG,1);
for p = 1:nPermG
    signs = (rand(nSub,1) > 0.5)*2 - 1;
    mu_align_perm(p) = mean(align_delta .* signs);
end
p_group_align_oneside = (sum(mu_align_perm >= mu_align_obs) + 1) / (nPermG + 1);

fprintf('\n=== Group-level ===\n');
fprintf('Mean cosND = %.4f, p(two-sided signflip)=%.5f\n', mu_cos_obs, p_group_cos_twoside);
fprintf('Mean (align - permBaseline) = %.4f, p(one-sided signflip)=%.5f\n', mu_align_obs, p_group_align_oneside);

results.subj = subj;
results.group.mu_cos_obs = mu_cos_obs;
results.group.p_cos_twoside = p_group_cos_twoside;
results.group.mu_align_delta_obs = mu_align_obs;
results.group.p_align_oneside = p_group_align_oneside;

% ============ 5) performance × orthogonality ============
% 注意：performance 需要“正确/错误”的定义才能算，这里只用于表现，不进入 vN/vD
perf = zeros(nSub,1);
for s = 1:nSub
    D = infer_data{s};
    Cue        = double(D.Cue(:));
    RuleChoice = double(D.RuleChoice(:));
    errTrial = double(RuleChoice ~= Cue); % 仅用于 accuracy
    perf(s) = 1 - mean(errTrial);
end

absCos  = abs(cos_vec);
ang_vec = [subj.angND_obs]';
dev90   = abs(ang_vec - 90);

% Spearman
[r_absCos, p_absCos] = corr(perf, absCos, 'Type','Spearman', 'Rows','complete');
[r_dev90,  p_dev90 ] = corr(perf, dev90,  'Type','Spearman', 'Rows','complete');

% permutation for correlation
nPermCorr = 20000;
rperm_absCos = zeros(nPermCorr,1);
rperm_dev90  = zeros(nPermCorr,1);
for p = 1:nPermCorr
    idx = randperm(nSub);
    perfP = perf(idx);
    rperm_absCos(p) = corr(perfP, absCos, 'Type','Spearman', 'Rows','complete');
    rperm_dev90(p)  = corr(perfP, dev90,  'Type','Spearman', 'Rows','complete');
end
pperm_absCos = (sum(abs(rperm_absCos) >= abs(r_absCos)) + 1) / (nPermCorr + 1);
pperm_dev90  = (sum(abs(rperm_dev90)  >= abs(r_dev90 )) + 1) / (nPermCorr + 1);

fprintf('\n=== Perf × Orthogonality (Spearman + permutation) ===\n');
fprintf('Perf vs |cosND|:   rho = %.4f, p(param)=%.4g, p(perm)=%.5f\n', r_absCos, p_absCos, pperm_absCos);
fprintf('Perf vs |ang-90|:  rho = %.4f, p(param)=%.4g, p(perm)=%.5f\n', r_dev90,  p_dev90,  pperm_dev90);

results.behavior.perf_accuracy = perf;
results.behavior.corr.perf_vs_absCos = struct('rho', r_absCos, 'p_param', p_absCos, 'p_perm', pperm_absCos);
results.behavior.corr.perf_vs_dev90  = struct('rho', r_dev90,  'p_param', p_dev90,  'p_perm', pperm_dev90);

% ============ 6) group plots ============
alignDelta = align_delta;

figure('Name','Group: cosND distribution','Color','w');
histogram(cos_vec, 12);
xlabel('cos(v_{Nback}, v_{tDev})'); ylabel('Count');
title(sprintf('Group cosND (mean=%.3f)', mean(cos_vec)));

figure('Name','Group: |cosND| distribution','Color','w');
histogram(absCos, 12);
xlabel('|cos(v_{Nback}, v_{tDev})| (deviation from orthogonality)'); ylabel('Count');
title(sprintf('Group |cosND| (mean=%.3f)', mean(absCos)));

figure('Name','Group: angle(Nback,tDev) distribution','Color','w');
histogram(ang_vec, 12);
xlabel('angle(v_{Nback}, v_{tDev}) (deg)'); ylabel('Count');
title(sprintf('Group angle (mean=%.1f°)', mean(ang_vec)));

figure('Name','Group: alignDelta distribution','Color','w');
histogram(alignDelta, 12);
xlabel('align(obs) - mean(align(perm))'); ylabel('Count');
title(sprintf('Group alignDelta (mean=%.4f)', mean(alignDelta)));

figure('Name','Perf vs |cosND|','Color','w');
scatter(absCos, perf, 40, 'filled'); hold on;
lsline;
xlabel('|cos(v_{Nback}, v_{tDev})| (smaller = more orthogonal)');
ylabel('Performance (accuracy)');
title(sprintf('Spearman rho=%.3f, p(perm)=%.5f', r_absCos, pperm_absCos));
grid on;

figure('Name','Perf vs |angle-90|','Color','w');
scatter(dev90, perf, 40, 'filled'); hold on;
lsline;
xlabel('|angle(v_{Nback}, v_{tDev}) - 90°| (smaller = more orthogonal)');
ylabel('Performance (accuracy)');
title(sprintf('Spearman rho=%.3f, p(perm)=%.5f', r_dev90, pperm_dev90));
grid on;




%% ============ (ADD-ON) Performance × Orthogonality + Group plots ============
% Paste this block at the end of your script (after group-level fprintf), before the final "end".
% Assumes variables exist in workspace:
%   infer_data, subj, nSub
%   (optional) results struct already exists; if not, code will still run and just not store.

% -------- 1) performance (accuracy) --------
perf = zeros(nSub,1);
for s = 1:nSub
    D = infer_data{s};
    Cue        = double(D.Cue(:));
    RuleChoice = double(D.RuleChoice(:));
    errTrial = double(RuleChoice ~= Cue);   % only for performance
    perf(s) = 1 - mean(errTrial);           % accuracy
end

% -------- 2) orthogonality metrics --------
cosND  = [subj.cosND_obs]';           % cos between vNback and vtDev
absCos = abs(cosND);                  % deviation from orthogonality (smaller=more orthogonal)
angND  = [subj.angND_obs]';           % degrees
dev90  = abs(angND - 90);             % smaller=more orthogonal

% -------- 3) alignment delta (obs - perm baseline) --------
alignObs  = [subj.align_obs]';
alignBase = zeros(nSub,1);
for s = 1:nSub
    alignBase(s) = mean(subj(s).align_perm);
end
alignDelta = alignObs - alignBase;

% -------- 4) performance × orthogonality correlations --------
[r_absCos, p_absCos] = corr(perf, absCos, 'Type','Spearman', 'Rows','complete');
[r_dev90,  p_dev90 ] = corr(perf, dev90,  'Type','Spearman', 'Rows','complete');

% permutation test for correlation (shuffle perf across subjects)
nPermCorr = 20000;
rperm_absCos = zeros(nPermCorr,1);
rperm_dev90  = zeros(nPermCorr,1);
for p = 1:nPermCorr
    idx = randperm(nSub);
    perfP = perf(idx);
    rperm_absCos(p) = corr(perfP, absCos, 'Type','Spearman', 'Rows','complete');
    rperm_dev90(p)  = corr(perfP, dev90,  'Type','Spearman', 'Rows','complete');
end
pperm_absCos = (sum(abs(rperm_absCos) >= abs(r_absCos)) + 1) / (nPermCorr + 1);
pperm_dev90  = (sum(abs(rperm_dev90)  >= abs(r_dev90 )) + 1) / (nPermCorr + 1);

fprintf('\n=== Perf × Orthogonality (Spearman + permutation) ===\n');
fprintf('Perf vs |cosND|:   rho = %.4f, p(param)=%.4g, p(perm)=%.5f\n', r_absCos, p_absCos, pperm_absCos);
fprintf('Perf vs |ang-90|:  rho = %.4f, p(param)=%.4g, p(perm)=%.5f\n', r_dev90,  p_dev90,  pperm_dev90);

% optionally store into results if it exists
if exist('results','var')
    results.behavior.perf_accuracy = perf;
    results.behavior.alignDelta = alignDelta;
    results.behavior.corr.perf_vs_absCos = struct('rho', r_absCos, 'p_param', p_absCos, 'p_perm', pperm_absCos);
    results.behavior.corr.perf_vs_dev90  = struct('rho', r_dev90,  'p_param', p_dev90,  'p_perm', pperm_dev90);
end

%% -------- 5) group plots --------
% (a) cosND distribution
figure('Name','Group: cosND distribution','Color','w');
histogram(cosND, 12);
xlabel('cos(v_{Nback}, v_{tDev})'); ylabel('Count');
title(sprintf('Group cosND (mean=%.3f)', mean(cosND)));

% (b) |cosND| distribution
figure('Name','Group: |cosND| distribution','Color','w');
histogram(absCos, 12);
xlabel('|cos(v_{Nback}, v_{tDev})| (deviation from orthogonality)'); ylabel('Count');
title(sprintf('Group |cosND| (mean=%.3f)', mean(absCos)));

% (c) angle distribution
figure('Name','Group: angle(Nback,tDev) distribution','Color','w');
histogram(angND, 12);
xlabel('angle(v_{Nback}, v_{tDev}) (deg)'); ylabel('Count');
title(sprintf('Group angle (mean=%.1f°)', mean(angND)));

% (d) alignDelta distribution
figure('Name','Group: alignDelta distribution','Color','w');
histogram(alignDelta, 12);
xlabel('align(obs) - mean(align(perm))'); ylabel('Count');
title(sprintf('Group alignDelta (mean=%.4f)', mean(alignDelta)));

% (e) perf vs |cosND|
figure('Name','Perf vs |cosND|','Color','w');
scatter(absCos, perf, 40, 'filled'); hold on;
lsline;
xlabel('|cos(v_{Nback}, v_{tDev})| (smaller = more orthogonal)');
ylabel('Performance (accuracy)');
title(sprintf('Spearman rho=%.3f, p(perm)=%.5f', r_absCos, pperm_absCos));
grid on;

% (f) perf vs |angle-90|
figure('Name','Perf vs |angle-90|','Color','w');
scatter(dev90, perf, 40, 'filled'); hold on;
lsline;
xlabel('|angle(v_{Nback}, v_{tDev}) - 90°| (smaller = more orthogonal)');
ylabel('Performance (accuracy)');
title(sprintf('Spearman rho=%.3f, p(perm)=%.5f', r_dev90, pperm_dev90));
grid on;

end

