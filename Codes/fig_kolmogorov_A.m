%% FIG_KOLMOGOROV_A.M  -  Kolmogorov: detection and regulatory identification
%  Requires: results_kolmogorov.mat
%  Output:   ../Plots/fig_kolmogorov_A.pdf/png

clc; clear; close all;

results_path = fullfile(fileparts(mfilename('fullpath')), 'results_kolmogorov.mat');
if ~exist(results_path,'file')
    error('results_kolmogorov.mat not found.');
end
load(results_path);
mode_tex      = {'$a_{0,k_f}$','$a_{1,0}$','$a_{1,k_f}$','$a_{2,0}$', ...
                 '$a_{0,2k_f}$','$a_{1,2k_f}$','$a_{2,k_f}$'};
mode_tex_bare = {'a_{0,k_f}','a_{1,0}','a_{1,k_f}','a_{2,0}', ...
                 'a_{0,2k_f}','a_{1,2k_f}','a_{2,k_f}'};

fprintf('Loaded: %s\n', results_path);
fprintf('Best-pair (%s,%s)|Z=%s  tau*=%.1f  Theta=%.0f  p=%.4f\n', ...
    mode_tex_bare{xi_51},mode_tex_bare{yi_51},mode_tex_bare{i_large}, ...
    tau_star_bc,Theta_bc,pval_bc);

%% ---- Colour palette (Wong 2011, CBF-safe) --------------------------------
C.blue   = [0.000 0.447 0.698];
C.orange = [0.902 0.624 0.000];
C.red    = [0.835 0.369 0.000];
C.green  = [0.000 0.620 0.451];
C.grey   = [0.600 0.600 0.600];

div_cm = interp1(linspace(0,1,5), ...
    [C.blue;[0.3,0.6,0.9];[0.97,0.97,0.97];[0.99,0.68,0.25];C.red], linspace(0,1,256));

set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',9);

%% Figure A
fig = figure('Units','centimeters','Position',[2 2 17.4 12],'Color','w');
tl  = tiledlayout(2,3,'TileSpacing','compact','Padding','tight');

%% ---- (a) Energy dissipation with extreme events -------------------------
ax_a = nexttile(tl,1);
plot(t_vec, D_vec, '-', 'Color',C.grey, 'LineWidth',0.4); hold on;
scatter(t_vec(is_extreme), D_vec(is_extreme), 6, C.red, 'filled','MarkerFaceAlpha',0.7);
yline(D_extreme,'--','Color',C.red,'LineWidth',1.0, ...
    'Label','$D_e$','Interpreter','latex','FontSize',8,'LabelVerticalAlignment','bottom');
hold off;
xlabel('$t$ (t.u.)','Interpreter','latex');
ylabel('$D(t)$','Interpreter','latex');
title('(a) Energy dissipation','Interpreter','latex','FontSize',9.5);
xlim([t_vec(1) t_vec(end)]); box on; grid on;

%% ---- (b) Theta regulatory gain spectrum (bars + 90% CI) ----------------
ax_b = nexttile(tl,2);
bar_clr = repmat(C.grey,n_modes,1); bar_clr(i_large,:) = C.green;
hold on;
for m=1:n_modes
    bar(m,Theta_reg(m),0.68,'FaceColor',bar_clr(m,:),'EdgeColor','none','FaceAlpha',0.9);
end
errorbar(1:n_modes, Theta_reg, max(0,Theta_reg-Theta_CI_lo), max(0,Theta_CI_hi-Theta_reg), ...
    'k.','LineWidth',0.9,'CapSize',4);
yline(3,'--k','LineWidth',1,'Label','$\Theta\!=\!3$','Interpreter','latex','FontSize',8);
hold off;
xticks(1:n_modes); xticklabels(mode_tex);
ax_b.XAxis.TickLabelInterpreter='latex'; xtickangle(30);
ylabel('$\max\,\Theta_\Sigma$  (90\% CI)','Interpreter','latex');
title('(b) Regulatory gain spectrum','Interpreter','latex','FontSize',9.5);
ylim([0, max(Theta_reg)*1.35]); box on; grid on;

%% ---- (c) Bootstrap rank stability ---------------------------------------
ax_c = nexttile(tl,3);
rank1_per_mode = zeros(1,n_modes);
for m=1:n_modes
    rank1_per_mode(m) = mean(boot_ranks(:,1)==m)*100;
end
bar_clr2 = repmat(C.grey,n_modes,1); bar_clr2(i_large,:) = C.green;
hold on;
for m=1:n_modes
    bar(m,rank1_per_mode(m),0.68,'FaceColor',bar_clr2(m,:),'EdgeColor','none','FaceAlpha',0.9);
end
text(i_large, rank1_per_mode(i_large)+5.5, sprintf('%.0f%%',rank1_per_mode(i_large)), ...
    'HorizontalAlignment','center','FontSize',8.5,'Color',C.green,'FontWeight','bold', ...
    'Interpreter','none');
hold off;
xticks(1:n_modes); xticklabels(mode_tex);
ax_c.XAxis.TickLabelInterpreter='latex'; xtickangle(30);
ylabel('Ranked \#1 in bootstrap (\%)','Interpreter','latex');
title('(c) Bootstrap rank stability','Interpreter','latex','FontSize',9.5);
ylim([0 115]); box on; grid on;

%% ---- (d) MIz at tau* - best pair, nodes_spod ---------------------------
ax_d = nexttile(tl,4);
P_v    = length(MIz_bc);
m_pct  = (1:P_v)'/P_v*100;
MI_dev = MIz_bc - mean(MIz_bc);
scatter(m_pct, MIz_bc, 22, MI_dev,'filled','MarkerFaceAlpha',0.9); hold on;
yline(mean(MIz_bc),'--','Color',C.grey,'LineWidth',0.9);
plot(m_pct, movmean(MIz_bc,max(3,floor(P_v/5))),'-k','LineWidth',1.4);
colormap(ax_d,div_cm); clim([-max(abs(MI_dev))-eps, max(abs(MI_dev))+eps]);
cb = colorbar('Location','west');
pos = cb.Position;
pos(4) = 0.5*pos(4);   % reduce height to 50%
pos(2) = pos(2) + 0.14; % recenter vertically
cb.Position = pos;
hold off;
xlabel('Percentile of $|a_{1,0}|$','Interpreter','latex');
ylabel('$\mathrm{MI}(X,Y\,|\,z_m)$','Interpreter','latex');
title(sprintf('(d) TRIM: $(%s,\\,%s)\\,|\\,Z\\!=\\!a_{1,0}$', ...
    mode_tex_bare{xi_51}, mode_tex_bare{yi_51}), ...
    'Interpreter','latex','FontSize',9.5);
text(0.24,0.91,sprintf('$\\Theta_\\Sigma=%.0f$,  $p=%.3f$',Theta_bc,pval_bc), ...
    'Units','norm','Interpreter','latex','FontSize',9,'Color',C.green, ...
    'BackgroundColor','w','Margin',1.5);
box on; grid on;

%% ---- (e) TRIM MI profile for Z = a_{1,kf} ---------------------------------
ax_e = nexttile(tl,5);
xi_e  = dir_triple_idx(2,1);
yi_e  = dir_triple_idx(2,2);
zi_e  = dir_triple_idx(2,3);
MIz_e    = MIz_dir_cell{2};
P_e      = length(MIz_e);
m_pct_e  = (1:P_e)'/P_e*100;
MI_dev_e = MIz_e - mean(MIz_e);
% Th_e, pv_e pre-computed by precompute_fig_data.m
scatter(m_pct_e, MIz_e, 22, MI_dev_e,'filled','MarkerFaceAlpha',0.9); hold on;
yline(mean(MIz_e),'--','Color',C.grey,'LineWidth',0.9);
plot(m_pct_e, movmean(MIz_e,max(3,floor(P_e/5))),'-k','LineWidth',1.4);
colormap(ax_e, div_cm);
clim(ax_e,[-max(abs(MI_dev_e))-eps, max(abs(MI_dev_e))+eps]);
cb_e = colorbar(ax_e,'Location','east');
pos_e = cb_e.Position; pos_e(4)=0.5*pos_e(4); pos_e(2)=pos_e(2)+0.12;
cb_e.Position=pos_e;
hold off;
xlabel(sprintf('Percentile of $|%s|$',mode_tex_bare{zi_e}),'Interpreter','latex');
ylabel('$\mathrm{MI}(X,Y\,|\,z_m)$','Interpreter','latex');
title(sprintf('(e) TRIM: $(%s,\\,%s)\\,|\\,Z\\!=\\!%s$', ...
    mode_tex_bare{xi_e},mode_tex_bare{yi_e},mode_tex_bare{zi_e}), ...
    'Interpreter','latex','FontSize',9.5);
text(0.24,0.91,sprintf('$\\Theta_\\Sigma=%.0f$,  $p=%.3f$',Th_e,pv_e), ...
    'Units','norm','Interpreter','latex','FontSize',9,'Color',C.orange, ...
    'BackgroundColor','w','Margin',1.5);
box on; grid on;

%% ---- (f) Sigma(tau) lag sweep - both triads ----------------------------
ax_f = nexttile(tl,6);
h_fs = plot(tau_vals, Sigma_sweep,'-o','Color',C.orange,'LineWidth',1.2, ...
    'MarkerSize',4,'MarkerFaceColor',C.orange, ...
    'DisplayName',sprintf('$(%s,\\,%s)$',mode_tex_bare{i_mean},mode_tex_bare{i_forced}));
hold on;
h_bp = plot(tau_vals, Sigma_sweep_bc,'-o','Color',C.blue,'LineWidth',1.2, ...
    'MarkerSize',4,'MarkerFaceColor',C.blue, ...
    'DisplayName',sprintf('$(%s,\\,%s)$',mode_tex_bare{xi_51},mode_tex_bare{yi_51}));
hold off;
xlabel('$\tau$ (t.u.)','Interpreter','latex');
ylabel('$\Sigma(\tau)$','Interpreter','latex');
title('(f) Lag sweep  $[\,\cdot\,,\,\cdot\,|\,Z_\tau\!=\!a_{1,0}]$', ...
    'Interpreter','latex','FontSize',9.5);
legend([h_fs,h_bp],'Box','off','FontSize',7.5,'Location','northeast','Interpreter','latex');
box on; grid on;

%% ---- Save ---------------------------------------------------------------
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Plots');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
print(fig, fullfile(out_dir,'fig_kolmogorov_A.pdf'), '-dpdf', '-vector', '-r300');
print(fig, fullfile(out_dir,'fig_kolmogorov_A.png'), '-dpng', '-r300');
fprintf('Saved: Plots/fig_kolmogorov_A\n');
