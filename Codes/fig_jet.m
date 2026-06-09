%% FIG_JET_COMBINED.M - Turbulent jet: single combined analysis figure
%  Requires: results_jet.mat  (run_analysis_jet.m)
%  Output:   ../Plots/fig_jet_combined.pdf/png

clc; clear; close all;

results_path = fullfile(fileparts(mfilename('fullpath')), 'results_jet.mat');
if ~exist(results_path,'file')
    error('results_jet.mat not found. Run run_analysis_jet.m first.');
end
load(results_path);
fprintf('Loaded: %s\n', results_path);

%% ---- Colour palette (Wong 2011, CBF-safe) ---------------------------------
C.blue   = [0.000 0.447 0.698];
C.orange = [0.902 0.624 0.000];
C.red    = [0.835 0.369 0.000];
C.green  = [0.000 0.620 0.451];
C.grey   = [0.600 0.600 0.600];

div_cm = interp1(linspace(0,1,5), ...
    [C.blue;[0.3,0.6,0.9];[0.97,0.97,0.97];[0.99,0.68,0.25];C.red], linspace(0,1,256));
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',9);

%% ---- Triad 1: (f_X, f_Y) | Z = f_Z  [bxi, byi | bz] ----------------------
MIz_1   = res.MIz_cell{bl, bz};
P_1     = length(MIz_1);
m_pct_1 = (1:P_1)' / P_1 * 100;
Th_1    = res.Theta(bl, bz);
dir_1   = ternary(mean(MIz_1(1:5)) < mean(MIz_1(end-4:end)), 'INC', 'DEC');

X1 = res.nodes(:,bxi);  Y1 = res.nodes(:,byi);  Z1 = res.nodes(:,bz);
t_sp  = res.t_vec(:);
q25_1 = prctile(Z1,25);  q75_1 = prctile(Z1,75);
lo1 = Z1 <= q25_1;  hi1 = Z1 >= q75_1;
r_lo1 = corr(X1(lo1), Y1(lo1));  r_hi1 = corr(X1(hi1), Y1(hi1));
p_lo1 = polyfit(X1(lo1), Y1(lo1), 1);
p_hi1 = polyfit(X1(hi1), Y1(hi1), 1);
xy_max1 = max([X1; Y1]) * 1.05;
x_line1 = linspace(0, xy_max1, 80);
fprintf('Triad 1 (%s,%s|Z=%s):  Theta=%.1f  %s\n', ...
    f_maths{bxi}, f_maths{byi}, f_maths{bz}, Th_1, dir_1);
fprintf('  r(Q1)=%+.3f  r(Q4)=%+.3f\n', r_lo1, r_hi1);

%% ---- Triad 2: (f_X, f_Z) | Z = f_Y  [pA, pB | byi] -----------------------
pA = min(bxi,bz);  pB = max(bxi,bz);
edge_l2 = find(res.edges(:,1)==pA & res.edges(:,2)==pB, 1);
zi_2    = byi;
MIz_2   = res.MIz_cell{edge_l2, zi_2};
P_2     = length(MIz_2);
m_pct_2 = (1:P_2)' / P_2 * 100;
Th_2    = res.Theta(edge_l2, zi_2);
dir_2   = ternary(mean(MIz_2(1:5)) < mean(MIz_2(end-4:end)), 'INC', 'DEC');

X2 = res.nodes(:,pA);  Y2 = res.nodes(:,pB);  Z2 = res.nodes(:,zi_2);
q25_2 = prctile(Z2,25);  q75_2 = prctile(Z2,75);
lo2 = Z2 <= q25_2;  hi2 = Z2 >= q75_2;
r_lo2 = corr(X2(lo2), Y2(lo2));  r_hi2 = corr(X2(hi2), Y2(hi2));
p_lo2 = polyfit(X2(lo2), Y2(lo2), 1);
p_hi2 = polyfit(X2(hi2), Y2(hi2), 1);
xy_max2 = max([X2; Y2]) * 1.05;
x_line2 = linspace(0, xy_max2, 80);
fprintf('Triad 2 (%s,%s|Z=%s):  Theta=%.1f  %s\n', ...
    f_maths{pA}, f_maths{pB}, f_maths{zi_2}, Th_2, dir_2);
fprintf('  r(Q1)=%+.3f  r(Q4)=%+.3f\n', r_lo2, r_hi2);

%% Figure
fig = figure('Units','centimeters','Position',[2 2 17.4 13],'Color','w');
tl  = tiledlayout(3,8,'TileSpacing','compact','Padding','tight');

%% (a) SPOD eigenspectrum  [row 1, cols 1-4] --------------------------------
ax_a = nexttile(1,[1 2]);
semilogy(ax_a, res.f_spod, res.L(:,1)+1e-20, '-', 'Color',C.blue, 'LineWidth',1.2); hold on;
if size(res.L,2) > 1
    semilogy(ax_a, res.f_spod, res.L(:,2)+1e-20, '-', 'Color',C.grey, 'LineWidth',0.7);
end
ymax_L = max(res.L(:,1)) * 2;
for fi_ = 1
    xline(ax_a, res.f_targets(fi_), '--', 'Color',f_colors{fi_}, 'LineWidth',0.9);
    text(ax_a, res.f_targets(fi_)+0.008, 0.1, f_labels{fi_}, ...
        'Interpreter','latex','FontSize',7.5,'Color',f_colors{fi_},'Rotation',90);
end
for fi_ = 2:nF
    xline(ax_a, res.f_targets(fi_), '--', 'Color',f_colors{fi_}, 'LineWidth',0.9);
    text(ax_a, res.f_targets(fi_)+0.008, ymax_L*0.38^fi_, f_labels{fi_}, ...
        'Interpreter','latex','FontSize',7.5,'Color',f_colors{fi_},'Rotation',90);
end
hold off;
xlabel(ax_a,'Frequency (St)','Interpreter','none');
ylabel(ax_a,'SPOD eigenvalue $\lambda$','Interpreter','latex');
title(ax_a,'(a) SPOD eigenspectrum','Interpreter','Latex');
legend(ax_a,{'Mode 1','Mode 2'},'Box','off','FontSize',7.5,'Location','southwest');
xlim(ax_a,[0 1.1]);  box(ax_a,'on');  grid(ax_a,'on');

%% (b) BMD bispectrum  [row 1, cols 5-8] ------------------------------------
ax_b = nexttile(3, [1 3]);
if has_bmd_full
    L_full = abs(L_bmd_full);  L_full(isnan(L_bmd_full)) = 0;
    [f1g, f2g] = ndgrid(f_bmd_full/f0, f_bmd_full/f0);
    L_log  = log10(L_full + eps);
    L_log(isnan(L_bmd_full)) = NaN;
    mask_r1 = f1g >= 0 & f2g >= 0 & f1g >= f2g;
    L_disp  = L_log;  L_disp(~mask_r1 | L_full==0) = NaN;
    pcolor(ax_b, f1g, f2g, L_disp);  shading flat;  hold on;
    for fi_ = 1:nF
        fn = res.f_targets(fi_) / f0;
        plot(ax_b, [fn fn],[0 fn],'--','Color',f_colors{fi_},'LineWidth',0.6);
        plot(ax_b, [0 fn],[fn fn],'--','Color',f_colors{fi_},'LineWidth',0.6);
    end
    hold off;
    colormap(ax_b, parula);
    cb_b = colorbar(ax_b,'Location','west');
    ws = warning('off','MATLAB:callback:error');
    cb_b.Label.Interpreter = 'latex';
    cb_b.Label.String = '$\log_{10}|\lambda_\mathrm{BMD}|$';
    cb_b.FontSize = 7;
    warning(ws);
    xlim(ax_b,[0 3]);  ylim(ax_b,[0 3]);
    set(ax_b,'XTick',[1 2 3],'XTickLabel',{'$f_0$','$2f_0$','$3f_0$'}, ...
             'YTick',[1 2 3],'YTickLabel',{'$f_0$','$2f_0$','$3f_0$'}, ...
             'TickLabelInterpreter','latex');
    mask_vis = mask_r1 & f1g <= 3 & f2g <= 3 & L_full > 0;
    clim(ax_b,[max(L_log(mask_vis))-3, max(L_log(mask_vis))]);
end
xlabel(ax_b,'$f_1$','Interpreter','latex');
ylabel(ax_b,'$f_2$','Interpreter','latex');
title(ax_b,'(b) BMD bispectrum','Interpreter','Latex');
box(ax_b,'on');  grid(ax_b,'on');

%% (c) TRIM heatmap  [row 1, cols 9-12] -------------------------------------
ax_dt = nexttile(6,[1,3]);
Th_max_d = max(Th_disp(~isnan(Th_disp)));
Th_show  = Th_disp;  Th_show(isnan(Th_show) | Th_show < 0) = 0;
imagesc(ax_dt, Th_show);
colormap(ax_dt, parula);
clim(ax_dt,[0, Th_max_d]);
hold(ax_dt,'on');
for l = 1:nEdge
    for zi = 1:nF
        if isnan(Th_disp(l,zi))
            fill(ax_dt, [zi-0.5,zi+0.5,zi+0.5,zi-0.5,zi-0.5], ...
                        [l-0.5,  l-0.5,  l+0.5,  l+0.5,  l-0.5], ...
                [0.82 0.82 0.82], 'EdgeColor','none');
            continue
        end
        tc = ternary(Th_disp(l,zi) > 0.6*Th_max_d, 'k', 'w');
        text(ax_dt, zi, l, sprintf('%.1f',Th_disp(l,zi)), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',9,'Color',tc,'FontWeight','bold','Interpreter','none');
    end
end
hold(ax_dt,'off');
set(ax_dt,'XTick',1:nF,'XTickLabel',f_labels,'TickLabelInterpreter','latex', ...
          'YTick',1:nEdge,'YTickLabel',pair_lbl,'TickLabelInterpreter','latex', ...
          'FontSize',9);
xlabel(ax_dt,'Regulator $f_Z$','Interpreter','latex');
title(ax_dt,'(c) TRIM $\Theta_\Sigma$ per triad','Interpreter','latex');
axis(ax_dt,'normal');  % release imagesc square-pixel constraint
%box(ax_dt,'on');

%% ---- (d) MI profile - Triad 1  [row2, cols1-3] ---------------------------
ax_d = nexttile(tl,9,[1,4]);
MI_dev1 = MIz_1 - mean(MIz_1);
clim_d  = max(abs(MI_dev1)) + eps;
scatter(ax_d, m_pct_1, MIz_1, 20, MI_dev1,'filled','MarkerFaceAlpha',0.85); hold on;
yline(mean(MIz_1),'--','Color',C.grey,'LineWidth',0.9,'HandleVisibility','off');
plot(ax_d, m_pct_1, movmean(MIz_1,max(3,floor(P_1/5))),'-k','LineWidth',1.6);
hold off;
colormap(ax_d, div_cm); clim(ax_d,[-clim_d,clim_d]);
xlabel(ax_d,sprintf('Percentile of $|a_{%s}|$',f_maths{bz}),'Interpreter','latex');
ylabel(ax_d,'$\mathrm{MI}(X,Y\,|\,z_m)$','Interpreter','latex');
title(ax_d,sprintf('(d) $(%s,\\,%s)\\,|\\,Z\\!=\\!%s$', ...
    f_maths{bxi},f_maths{byi},f_maths{bz}),'Interpreter','latex','FontSize',9.5);
text(ax_d,0.53,0.85,sprintf('$\\Theta_\\Sigma=%.0f$  (%s)',Th_1,dir_1), ...
    'Units','norm','Interpreter','latex','FontSize',9,'Color',C.green, ...
    'FontWeight','bold','BackgroundColor','w','Margin',1.5);
box on; grid on;

%% ---- (e) MI profile - Triad 2  [row2, cols4-6] ---------------------------
ax_e = nexttile(tl,13,[1,4]);
MI_dev2 = MIz_2 - mean(MIz_2);
clim_e  = max(abs(MI_dev2)) + eps;
scatter(ax_e, m_pct_2, MIz_2, 20, MI_dev2,'filled','MarkerFaceAlpha',0.85); hold on;
yline(mean(MIz_2),'--','Color',C.grey,'LineWidth',0.9,'HandleVisibility','off');
plot(ax_e, m_pct_2, movmean(MIz_2,max(3,floor(P_2/5))),'-k','LineWidth',1.6);
hold off;
colormap(ax_e, div_cm); clim(ax_e,[-clim_e,clim_e]);
xlabel(ax_e,sprintf('Percentile of $|a_{%s}|$',f_maths{zi_2}),'Interpreter','latex');
%set(ax_e,'YTickLabel',[]);
title(ax_e,sprintf('(e) $(%s,\\,%s)\\,|\\,Z\\!=\\!%s$', ...
    f_maths{pA},f_maths{pB},f_maths{zi_2}),'Interpreter','latex','FontSize',9.5);
text(ax_e,0.03,0.85,sprintf('$\\Theta_\\Sigma=%.0f$  (%s)',Th_2,dir_2), ...
    'Units','norm','Interpreter','latex','FontSize',9,'Color',C.orange, ...
    'FontWeight','bold','BackgroundColor','w','Margin',1.5);
box on; grid on;

%% ---- (f) |Z1(t)| time series with Q1/Q4 coloring  [row3, cols1-2] -------
ax_f = nexttile(tl,17,[1,2]);
plot(ax_f, t_sp, Z1,'-','Color',C.grey,'LineWidth',0.4); hold on;
plot(ax_f, t_sp(lo1), Z1(lo1),'.','MarkerSize',3,'Color',C.blue);
plot(ax_f, t_sp(hi1), Z1(hi1),'.','MarkerSize',3,'Color',C.red);
yline(q25_1,'--','Color',C.blue,'LineWidth',0.9,'HandleVisibility','off');
yline(q75_1,'--','Color',C.red, 'LineWidth',0.9,'HandleVisibility','off');
hold off;
xlabel(ax_f,'$t\,D/U_j$','Interpreter','latex');
ylabel(ax_f,sprintf('$|a_{%s}|$',f_maths{bz}),'Interpreter','latex');
title(ax_f,sprintf('(f) $Z\\!=\\!%s$ amplitude',f_maths{bz}),'Interpreter','latex','FontSize',9.5);
xlim(ax_f,[t_sp(1) t_sp(end)]); ylim(ax_f,[0 max(Z1)*1.05]); box on; grid on;
text(ax_f,0.80,0.15,'Q1','Units','norm','Color',C.blue,'FontSize',8,'FontWeight','bold','Interpreter','none');
text(ax_f,0.01,0.88,'Q4','Units','norm','Color',C.red, 'FontSize',8,'FontWeight','bold','Interpreter','none');

%% ---- (g) Slope lines only - Triad 1  [row3, col3] -----------------------
ax_g = nexttile(tl,19,[1 2]);
x_g  = linspace(0, xy_max1, 80);
plot(ax_g, x_g, max(0,polyval(p_lo1,x_g)),'-','Color',C.blue,'LineWidth',2.2); hold on;
plot(ax_g, x_g, max(0,polyval(p_hi1,x_g)),'-','Color',C.red, 'LineWidth',2.2);
hold off;
text(ax_g,0.05,0.88,sprintf('Q1: $r\\!=\\!%.2f$',r_lo1), ...
    'Units','norm','Interpreter','latex','FontSize',8,'Color',C.blue,'FontWeight','bold');
text(ax_g,0.05,0.68,sprintf('Q4: $r\\!=\\!%.2f$',r_hi1), ...
    'Units','norm','Interpreter','latex','FontSize',8,'Color',C.red,'FontWeight','bold');
xlabel(ax_g,sprintf('$|a_{%s}|$',f_maths{bxi}),'Interpreter','latex');
ylabel(ax_g,sprintf('$|a_{%s}|$',f_maths{byi}),'Interpreter','latex');
title(ax_g,'(g)','Interpreter','Latex','FontSize',9.5);
xlim(ax_g,[0 xy_max1]); ylim(ax_g,[0 xy_max1]);
grid on;

%% ---- (h) |Z2(t)| time series with Q1/Q4 coloring  [row3, cols4-5] -------
ax_h = nexttile(tl,21,[1,2]);
plot(ax_h, t_sp, Z2,'-','Color',C.grey,'LineWidth',0.4); hold on;
plot(ax_h, t_sp(lo2), Z2(lo2),'.','MarkerSize',3,'Color',C.blue);
plot(ax_h, t_sp(hi2), Z2(hi2),'.','MarkerSize',3,'Color',C.red);
yline(q25_2,'--','Color',C.blue,'LineWidth',0.9,'HandleVisibility','off');
yline(q75_2,'--','Color',C.red, 'LineWidth',0.9,'HandleVisibility','off');
hold off;
xlabel(ax_h,'$t\,D/U_j$','Interpreter','latex');
ylabel(ax_h,sprintf('$|a_{%s}|$',f_maths{zi_2}),'Interpreter','latex');
title(ax_h,sprintf('(h) $Z\\!=\\!%s$ amplitude',f_maths{zi_2}),'Interpreter','latex','FontSize',9.5);
xlim(ax_h,[t_sp(1) t_sp(end)]); ylim(ax_h,[0 max(Z2)*1.05]); box on; grid on;
text(ax_h,0.80,0.15,'Q1','Units','norm','Color',C.blue,'FontSize',8,'FontWeight','bold','Interpreter','none');
text(ax_h,0.01,0.88,'Q4','Units','norm','Color',C.red, 'FontSize',8,'FontWeight','bold','Interpreter','none');

%% ---- (i) Slope lines only - Triad 2  [row3, col6] -----------------------
ax_i = nexttile(tl,23,[1,2]);
x_i  = linspace(0, xy_max2, 80);
plot(ax_i, x_i, max(0,polyval(p_lo2,x_i)),'-','Color',C.blue,'LineWidth',2.2); hold on;
plot(ax_i, x_i, max(0,polyval(p_hi2,x_i)),'-','Color',C.red, 'LineWidth',2.2);
hold off;
text(ax_i,0.05,0.88,sprintf('Q1: $r\\!=\\!%.2f$',r_lo2), ...
    'Units','norm','Interpreter','latex','FontSize',8,'Color',C.blue,'FontWeight','bold');
text(ax_i,0.05,0.68,sprintf('Q4: $r\\!=\\!%.2f$',r_hi2), ...
    'Units','norm','Interpreter','latex','FontSize',8,'Color',C.red,'FontWeight','bold');
xlabel(ax_i,sprintf('$|a_{%s}|$',f_maths{pA}),'Interpreter','latex');
ylabel(ax_i,sprintf('$|a_{%s}|$',f_maths{pB}),'Interpreter','latex');
title(ax_i,'(i)','Interpreter','Latex','FontSize',9.5);
xlim(ax_i,[0 xy_max2]); ylim(ax_i,[0 xy_max2]);
grid on;

%% ---- Shared y-axis for MI panels -----------------------------------------
linkaxes([ax_d, ax_e],'y');

%% ---- Save ----------------------------------------------------------------
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Plots');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
print(fig, fullfile(out_dir,'fig_jet_combined.pdf'), '-dpdf', '-vector', '-r300');
print(fig, fullfile(out_dir,'fig_jet_combined.png'), '-dpng', '-r300');
fprintf('Saved: Plots/fig_jet_combined\n');

%% ---- Helpers -------------------------------------------------------------
function s = ternary(c,a,b), if c, s=a; else, s=b; end, end
