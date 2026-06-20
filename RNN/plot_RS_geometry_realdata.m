function plot_RS_geometry_realdata(infer_data, betaE_all, betaD_all, useE, useD, readout_lambda, results)
% r/s style semi-circle plots (REAL DATA)
% r: SE = readout for rnn_pre
% s: SE = readout for SW
% Curves: NULL distributions from trial-level permutation
% Dots: observed subject angles
% Output: Fig_RS_geometry_REALDATA.svg

thetaMax = 90; nbins = 45;
nPerm_null = 10000; nPerm_stat = 10000;
pmin = 1e-4; rng(1);

% ---- compute r & s ----
[obsE_r, obsD_r, nullE_r, nullD_r] = collect_angles_trialperm( ...
    infer_data, betaE_all, betaD_all, useE, useD, readout_lambda, nPerm_null, 'mu_switch_estimated');
[obsE_s, obsD_s, nullE_s, nullD_s] = collect_angles_trialperm( ...
    infer_data, betaE_all, betaD_all, useE, useD, readout_lambda, nPerm_null, 'P(switch)');

[thetaC, pE_r] = prob_curve(nullE_r, thetaMax, nbins); [~, pD_r] = prob_curve(nullD_r, thetaMax, nbins);
[~,      pE_s] = prob_curve(nullE_s, thetaMax, nbins); [~, pD_s] = prob_curve(nullD_s, thetaMax, nbins);

rE_r = log_radius(pE_r, pmin); rD_r = log_radius(pD_r, pmin);
rE_s = log_radius(pE_s, pmin); rD_s = log_radius(pD_s, pmin);

rObsE_r = 0.90 + 0.10*rand(size(obsE_r)); rObsD_r = 0.90 + 0.10*rand(size(obsD_r));
rObsE_s = 0.90 + 0.10*rand(size(obsE_s)); rObsD_s = 0.90 + 0.10*rand(size(obsD_s));

muE_r = mean(obsE_r,'omitnan'); muD_r = mean(obsD_r,'omitnan');
muE_s = mean(obsE_s,'omitnan'); muD_s = mean(obsD_s,'omitnan');

p_pair_r = paired_perm_p(obsE_r - obsD_r, nPerm_stat);
p_pair_s = paired_perm_p(obsE_s - obsD_s, nPerm_stat);

cE = [250/255 216/255 36/255]; cD = [167/255 207/255 140/255];

% ---- fprintf notes (替代原左下角 annotation panel) ----
fprintf('\n=== RS geometry (trial-permutation NULL) ===\n');
fprintf('Curves: NULL P(theta) from trial permutation (nPerm_null=%d)\n', nPerm_null);
fprintf('Dots: observed subject theta(O,SE). Radial ticks: relative prob (normalized), pmin=%g\n', pmin);
fprintf('[r] SE=confidence | N=%d | mean theta_E=%.1f, theta_D=%.1f | paired p=%.4g\n', numel(obsE_r), muE_r, muD_r, p_pair_r);
fprintf('[s] SE=P(switch)      | N=%d | mean theta_E=%.1f, theta_D=%.1f | paired p=%.4g\n', numel(obsE_s), muE_s, muD_s, p_pair_s);
% is the muD_r larger than permu distri? is the muE_r larger than permu distri?
p_diff_vs_null = (sum(nullD_r >= muD_r) + 1) / (numel(nullD_r) + 1);
p_err_vs_null = (sum(nullE_r <= muE_r) + 1) / (numel(nullE_r) + 1);
fprintf('[r] SE=confidence | permP(nullD_r >= muD_r)=%.4f, permP(nullE_r <= muE_r)=%.4g\n', p_diff_vs_null, p_err_vs_null);


% ---- 1×3 layout: q | r | s ----
fig = figure('Color','w','Position',[80 80 1000 320]);
% tl = tiledlayout(fig,1,5,'Padding','compact','TileSpacing','compact');
tl = tiledlayout(fig,6,5,'Padding','compact','TileSpacing','compact'); % (rows, columns)
set(0, 'DefaultAxesFontName','Helvetica');
% set(gca,'FontSize',12)


% q (示意图：用 r 的均值角度定义 error/difficulty 轴，避免你之前的“不一致”)
axQ = nexttile(tl,6,[5 1]); cla(axQ); hold(axQ,'on'); axis(axQ,'equal'); axis(axQ,'off');
draw_q_cartesian_on(axQ, muE_r, muD_r);
title('A', 'FontSize',16,'Position',[0, 1.42, 0])

% r
axR = nexttile(tl,7,[5 2]); cla(axR); hold(axR,'on'); axis(axR,'equal'); axis(axR,'off');
draw_semicircle_axes_on(axR, thetaMax, pmin);
% plot_polarcurve_on(axR, thetaC, rE_r, cE, 3);  plot_polarcurve_on(axR, thetaC, rD_r, cD, 3);
plot_polarpoints_on(axR, obsE_r, rObsE_r, cE, 28); plot_polarpoints_on(axR, obsD_r, rObsD_r, cD, 28);
plot_ray_on(axR, muE_r, cE, '--', 3); plot_ray_on(axR, muD_r, cD, '--', 3);
annotate_panel(axR, '', 'SE = readout for confidence', numel(obsE_r), muE_r, muD_r, p_pair_r, cE, cD);
title('B', 'FontSize',16,'Position',[-1, 1.35, 0])

% s
% % % % % % axS = nexttile(tl,4,[1 2]); cla(axS); hold(axS,'on'); axis(axS,'equal'); axis(axS,'off');
% % % % % % draw_semicircle_axes_on(axS, thetaMax, pmin);
% % % % % % plot_polarcurve_on(axS, thetaC, rE_s, cE, 3);  plot_polarcurve_on(axS, thetaC, rD_s, cD, 3);
% % % % % % plot_polarpoints_on(axS, obsE_s, rObsE_s, cE, 28); plot_polarpoints_on(axS, obsD_s, rObsD_s, cD, 28);
% % % % % % plot_ray_on(axS, muE_s, cE, '--', 3); plot_ray_on(axS, muD_s, cD, '--', 3);
% % % % % % annotate_panel(axS, 's', 'SE = readout for P(switch)', numel(obsE_s), muE_s, muD_s, p_pair_s, cE, cD);

% % % % % % print(fig,'Fig_RS_geometry_BOLD.svg','-dsvg','-painters');
% % % % % % fprintf('Saved: Fig_RS_geometry_REALDATA.svg\n');
hold on


%% figure for overlap rate
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

% Style
FS = 11; LW = 1.2;
% ---------------- Figure: overlap----------------
axS = nexttile(tl,9,[5 2]); hold on

xB = proj_ratio(idxAlign);  yB = angle_w_plane(idxAlign);

histogram(xB,'BinWidth',0.05,'FaceAlpha',0.55,'FaceColor',[0.0671 0.1420 0.4189],...
    'EdgeColor',[1 1 1], 'EdgeAlpha', 0.55);
% yticks(0:2:8);% xlim([0 0.7]); xticks(0:0.2:0.4); 
yticks(0:4:16);

ylabel('Number of subjects')
xlabel('Readout alignment to span(error,difficulty)')

text(0.02,0.98,sprintf('n=%d (FDR sig=%d)\nproj=%.2f\\pm%.2f\nangle=%.1f\\pm%.1f', ...
    numel(xB), sum(sig_align_fdr), mean(xB,'omitnan'), sem(xB), mean(yB,'omitnan'), sem(yB)), ...
    'Units','normalized','VerticalAlignment','top','FontSize',FS)
set(gca,'FontSize',FS,'LineWidth',LW); box off

title('C','FontSize',16,'Position',[-0.05, 8.5, 0]) %==========BOLD=====
% title('C','FontSize',16,'Position',[-0.07, 6.5, 0])   %==========RNN=====

print(fig,'Fig_RS_geometry_BOLD.svg','-dsvg','-painters'); %==========BOLD=====
print(fig,'Fig_RS_geometry_BOLD.png','-dpng','-r300'); %==========BOLD=====

% print(fig,'Fig_RS_geometry_RNN.svg','-dsvg','-painters'); %==========RNN=====
% print(fig,'Fig_RS_geometry_RNN.png','-dpng','-r300'); %==========RNN=====
end







%% ============================ subfunctions ============================
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

function s = sem(x)
x = x(isfinite(x));
s = std(x)/sqrt(max(1,numel(x)));
end

function [obsE, obsD, nullE, nullD] = collect_angles_trialperm(infer_data, betaE_all, betaD_all, useE, useD, readout_lambda, nPerm_null, mode)
% mode: 'rnn_pre' or 'SW'

nSub = numel(infer_data);
obsE=[]; obsD=[]; nullE=[]; nullD=[];

for s = 1:nSub
    if ~(useE(s) && useD(s)), continue; end
    bE = betaE_all(s,:)'; bD = betaD_all(s,:)';
    if any(~isfinite(bE)) || any(~isfinite(bD)) || norm(bE)<1e-10 || norm(bD)<1e-10, continue; end

    D = infer_data{s};
    H = double(D.BOLD);

    if strcmpi(mode,'mu_switch_estimated')
        y = double(D.mu_switch_estimated(:));
        y = min(max(y,0),1);
    else
        % y = double(D.pr_of_switch(:) > 0); %******
        y = double(D.pr_of_switch(:));
    end


    idx = true(size(y));


    H = H(idx,:); y = y(idx);
    good = all(isfinite(H),2) & isfinite(y);
    H = H(good,:); y = y(good);

    if size(H,1) < 20, continue; end

    H = zscore(H,0,1);
    y = zscore(double(y));

    X = [ones(size(H,1),1), H];
    p = size(X,2);

    w_full = (X' * X + readout_lambda * eye(p)) \ (X' * y);
    SE = w_full(2:end);
    if any(~isfinite(SE)) || norm(SE)<1e-10, continue; end

    obsE(end+1,1) = acute_angle(bE, SE);
    obsD(end+1,1) = acute_angle(bD, SE);

    for k = 1:nPerm_null
        yperm = y(randperm(numel(y)));
        w_full_p = (X' * X + readout_lambda * eye(p)) \ (X' * yperm);
        SEp = w_full_p(2:end);
        nullE(end+1,1) = acute_angle(bE, SEp);
        nullD(end+1,1) = acute_angle(bD, SEp);
    end
end
end

function th = acute_angle(a,b)
a=a(:); b=b(:);
c = dot(a,b)/(norm(a)*norm(b)+eps);
c = max(min(c,1),-1);
th0 = acosd(c);
th = min(th0,180-th0);
end

function [thetaC, p] = prob_curve(angles, thetaMax, nbins)
angles = angles(isfinite(angles));
edges = linspace(0, thetaMax, nbins+1);
p = histcounts(angles, edges, 'Normalization','probability');
thetaC = (edges(1:end-1)+edges(2:end))/2;
if max(p)>0, p = p/max(p); end
p = max(p, 1e-4);
end

function r = log_radius(p, pmin)
p = max(min(p,1), pmin);
r = (log10(p) - log10(pmin)) / (0 - log10(pmin));
r = max(min(r,1),0);
end

function p = paired_perm_p(delta, nPerm)
delta = delta(isfinite(delta));
n = numel(delta);
if n<2, p=NaN; return; end
obs = mean(delta);
permDist = zeros(nPerm,1);
for k = 1:nPerm
    signFlip = (rand(n,1)>0.5)*2 - 1;
    permDist(k) = mean(delta .* signFlip);
end
p = mean(abs(permDist) >= abs(obs));
end

function annotate_panel(ax, letter, subtitle, N, muE, muD, p_pair, cE, cD)
text(ax, -1.1, 1.08, sprintf('%s   P(\\theta(O,SE))   (N=%d)', letter, N), 'FontSize', 13, 'FontWeight','bold');
text(ax, -1.1, 0.98, subtitle, 'FontSize', 11, 'FontWeight','bold');
% text(ax, -0.98, 0.86, 'Curves: NULL (trial permutation)', 'FontSize', 10);
text(ax, -1.1, 0.78, 'Dots: observed subject angles', 'FontSize', 10);
text(ax, -1.1, 0.66, sprintf('Mean \\theta_{errors}=%.1f°', muE), 'FontSize', 10, 'Color', cE);
text(ax, -1.1, 0.58, sprintf('Mean \\theta_{difficulty}=%.1f°', muD), 'FontSize', 10, 'Color', cD);
text(ax, -1.1, 0.42, sprintf('Paired perm p(\\theta_E-\\theta_D)=%.4f', p_pair), 'FontSize', 10, 'FontWeight','bold');

% small legend
plot(ax, [-1.1 -0.80], [0.34 0.34], '-', 'Color', cE, 'LineWidth', 4); text(ax,-0.78,0.34,'errors','FontSize',10);
plot(ax, [-1.1 -0.80], [0.26 0.26], '-', 'Color', cD, 'LineWidth', 4); text(ax,-0.78,0.26,'difficulty','FontSize',10);
end

function draw_q_cartesian_on(ax, muE, muD)
SE = [1 0];
E  = [cosd(muE) sind(muE)];
D  = [cosd(muD) sind(muD)];

plot(ax, [0 1.25],[0 0],'k-','LineWidth',1.5);
plot(ax, [0 0],[0 1.00],'k-','LineWidth',1.2);

quiver(ax,0,0,SE(1),SE(2),0,'LineWidth',2.8,'MaxHeadSize',0.25,'Color',[255/255 51/255 51/255]);
quiver(ax,0,0,E(1),E(2),0,'LineWidth',2.8,'MaxHeadSize',0.25,'Color',[250/255 216/255 36/255]);
quiver(ax,0,0,D(1),D(2),0,'LineWidth',2.8,'MaxHeadSize',0.25,'Color',[167/255 207/255 140/255]);

text(ax,1.02,0.03,'SE','Color',[255/255 51/255 51/255],'FontSize',12,'FontWeight','bold');
text(ax,E(1)+0.03,E(2)+0.03,'errors','Color',[250/255 216/255 36/255],'FontSize',12,'FontWeight','bold');
text(ax,D(1)+0.03,D(2)-0.06,'difficulty','Color',[167/255 207/255 140/255],'FontSize',12,'FontWeight','bold');

draw_arc_on(ax, 0.38, 0, muE, [250/255 216/255 36/255]);
draw_arc_on(ax, 0.52, 0, muD, [167/255 207/255 140/255]);

text(ax,0.54,0.20,sprintf('\\theta(errors,SE)=%.1f^\\circ',muE),'Color',[250/255 216/255 36/255],'FontSize',11);
text(ax,0.30,0.5,sprintf('\\theta(difficulty,SE)=%.1f^\\circ',muD),'Color',[167/255 207/255 140/255],'FontSize',11);

text(ax,0.00,1.3,'Angle definition','FontSize',13,'FontWeight','bold');
text(ax,0.00,1.15,'SE = readout for confidence','FontSize',12,'FontWeight','bold');

xlim(ax, [-0.001 1.03]); ylim(ax, [-0.001 1.10]);
end

function draw_arc_on(ax, r, ang1, ang2, col)
tt = linspace(deg2rad(ang1), deg2rad(ang2), 120);
plot(ax, r*cos(tt), r*sin(tt), '-', 'Color', col, 'LineWidth', 2.8);
end

function draw_semicircle_axes_on(ax, thetaMax, pmin)
th = linspace(0, deg2rad(thetaMax), 260);
plot(ax, cos(th), sin(th), 'Color', [0 0 0 0.25], 'LineWidth', 1.4);
plot(ax, [0 1], [0 0], 'Color', [0 0 0 0.25], 'LineWidth', 1.4);

angTicks = [0 30 60 90];
for a = angTicks
    plot(ax, [0 cosd(a)], [0 sind(a)], 'Color', [0 0 0 0.12], 'LineWidth', 1);
    text(ax, 1.06*cosd(a), 1.06*sind(a), sprintf('%d°',a), 'FontSize', 11, 'HorizontalAlignment','center');
end

tickP = [0.01 0.1 1];
tickP = max(tickP, pmin);
tickR = (log10(tickP)-log10(pmin)) / (0-log10(pmin));

for i = 1:numel(tickR)
    plot(ax, tickR(i)*cos(th), tickR(i)*sin(th), ':', 'Color', [0 0 0 0.18], 'LineWidth', 1.2);
    text(ax, -0.06, tickR(i), sprintf('%.2g', tickP(i)), 'FontSize', 10, 'HorizontalAlignment','right');
end

text(ax, -0.06, 0, 'P(\theta) (relative, normalized)', 'FontSize', 11, 'HorizontalAlignment','right', 'FontWeight','bold');
xlim(ax, [-1.05 1.12]); ylim(ax, [-0.10 1.12]);
end

function plot_polarcurve_on(ax, thetaDeg, r, col, lw)
x = r .* cosd(thetaDeg);
y = r .* sind(thetaDeg);
plot(ax, x, y, '-', 'Color', col, 'LineWidth', lw);
end

function plot_polarpoints_on(ax, thetaDeg, r, col, ms)
x = r .* cosd(thetaDeg);
y = r .* sind(thetaDeg);
scatter(ax, x, y, ms, 'MarkerFaceColor', col, 'MarkerEdgeColor','none', 'MarkerFaceAlpha',0.9);
end

function plot_ray_on(ax, thetaDeg, col, ls, lw)
plot(ax, [0 cosd(thetaDeg)], [0 sind(thetaDeg)], ls, 'Color', col, 'LineWidth', lw);
end