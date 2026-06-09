%% FIG_STUART_LANDAU.M  -  Stuart-Landau TRIM
%  Output: ../Plots/fig_SL.pdf/png

clc; clear; close all;

results_path = fullfile(fileparts(mfilename('fullpath')), 'results_sl.mat');
if ~exist(results_path,'file')
    error('results_sl.mat not found.');
end
load(results_path);

%% ---- Colour palette (Wong 2011, CBF-safe) --------------------------------
C.blue   = [0.000 0.447 0.698];
C.orange = [0.902 0.624 0.000];
C.red    = [0.835 0.369 0.000];
C.green  = [0.000 0.620 0.451];
C.grey   = [0.600 0.600 0.600];

set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',9);

%% Figure
fig = figure('Units','centimeters','Position',[2 2 17.4 11],'Color','w');
tl  = tiledlayout(2,2,'TileSpacing','compact','Padding','tight');

% ---- (a) Limit-cycle oscillation
ax1 = nexttile(tl);
plot(t_s, ts(1,idx_s), 'Color',C.blue,   'LineWidth',0.9,'DisplayName','$\mathrm{Re}(A)$'); hold on;
plot(t_s, ts(2,idx_s), 'Color',C.red, 'LineWidth',0.9,'DisplayName','$\mathrm{Im}(A)$');
hold off;
xlabel('$t$','Interpreter','latex'); ylabel('Amplitude','Interpreter','none');
title('(a) Limit-cycle oscillation $A(t)$','Interpreter','latex');
leg1 = legend('Box','off','FontSize',8,'Location','northeast'); leg1.Interpreter='latex';
box off; grid on;
ylim([-0.6 0.6]);

% ---- (b) Three candidate regulator signals
ax2 = nexttile(tl);
plot(t_s, Z_wn(idx_s)/sc,  '-', 'Color',C.grey,   'LineWidth',0.6,'DisplayName','White noise'); hold on;
plot(t_s, Z_U0(idx_s)/sc,  '-', 'Color',C.green,   'LineWidth',1.2,'DisplayName','$U_0$');
plot(t_s, Z_ou(idx_s)/sc,  '-', 'Color',C.orange,  'LineWidth',0.9,'DisplayName','Indep.\ OU');
hold off;
xlabel('$t$','Interpreter','latex'); ylabel('$Z/\sigma_{U_0}$','Interpreter','latex');
title('(b) Candidate regulator signals (normalised)','Interpreter','Latex');
leg2 = legend('Box','off','FontSize',8,'Location','Best'); leg2.Interpreter='latex';
ylim([-4 4]); box off; grid on;

% ---- (c) Null distribution: U0 (passes) vs OU (fails)
ax3 = nexttile(tl);
all_sig = [Sig_wn; SigRan_raw(:); Sig_ou; Sig_raw];
edges_c = linspace(min(all_sig)*0.95, max(all_sig)*1.05, 22);
[hc, he] = histcounts(SigRan_raw, edges_c, 'Normalization','probability');
bar(0.5*(he(1:end-1)+he(2:end)), hc, 1, ...
    'FaceColor',[0.82 0.82 0.82],'EdgeColor','none','FaceAlpha',0.85); hold on;
xline(Sig_wn,  '-', 'Color',C.grey,   'LineWidth',2.5, ...
    'Label','White noise', ...
    'Interpreter','none','FontSize',8,'LabelVerticalAlignment','bottom');
xline(Sig_ou,  '-', 'Color',C.orange, 'LineWidth',2.5, ...
    'Label','Indep.\ OU', ...
    'Interpreter','latex','FontSize',8,'LabelVerticalAlignment','bottom');
xline(Sig_raw, '-', 'Color',C.green,  'LineWidth',2.5, ...
    'Label','$U_0$', ...
    'Interpreter','latex','FontSize',8,'LabelVerticalAlignment','bottom');
hold off;
xlabel('$\Sigma$ (observed)','Interpreter','latex');
ylabel('Probability','Interpreter','none');
title('(c) Null test: $U_0$ rejects $H_0$; OU and noise do not','Interpreter','latex');
xlim([edges_c(1) edges_c(end)]);
box off; grid on;

% ---- (d) Theta summary bar: U0 > Indep. OU > White noise
ax4 = nexttile(tl);
hold on;
bar(1, Theta_scan(1), 0.6, 'FaceColor',C.green,  'EdgeColor','none','FaceAlpha',0.9);
bar(2, Theta_scan(3), 0.6, 'FaceColor',C.orange, 'EdgeColor','none','FaceAlpha',0.9);
bar(3, Theta_scan(2), 0.6, 'FaceColor',C.grey,   'EdgeColor','none','FaceAlpha',0.9);
yline(3,'--k','LineWidth',1.0,'Label','$\Theta=3$','Interpreter','latex','FontSize',8);
hold off;
xticks(1:3);
xticklabels({'$U_0$','Indep.\ OU','White noise'});
ax4.XAxis.TickLabelInterpreter = 'latex';
ylabel('$\Theta_\Sigma$ (shuffle null)','Interpreter','latex');
title('(d) TRIM: genuine $\gg$ correlated $\gg$ uncorrelated','Interpreter','latex');
ylim([0, max(Theta_scan)*1.45]);
box off; grid on;

%% ---- Save ----------------------------------------------------------------
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Plots');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
print(fig, fullfile(out_dir,'fig_SL.pdf'), '-dpdf', '-vector', '-r300');
print(fig, fullfile(out_dir,'fig_SL.png'), '-dpng', '-r300');
fprintf('Saved: Plots/fig_SL\n');
