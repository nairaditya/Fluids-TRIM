%% FIG_KOLMOGOROV_B.M  -  Kolmogorov: regulatory structure and flow-field gating
%  Requires:
%    results_kolmogorov.mat   (run_analysis_kolmogorov.m)
%  NOTE: Theta_full_kol must be pre-computed and saved in results_kolmogorov.mat.
%  Output: ../Plots/fig_kolmogorov_B.pdf/png

clc; clear; close all;

%% ---- Load results --------------------------------------------------------
results_path = fullfile(fileparts(mfilename('fullpath')), 'results_kolmogorov.mat');
if ~exist(results_path,'file')
    error('results_kolmogorov.mat not found.');
end
load(results_path);
fprintf('Loaded: %s\n', results_path);

%% ---- Colour palette (Wong 2011, CBF-safe) --------------------------------
C.blue   = [0.000 0.447 0.698];
C.orange = [0.902 0.624 0.000];
C.red    = [0.835 0.369 0.000];
C.green  = [0.000 0.620 0.451];
C.grey   = [0.600 0.600 0.600];

div_cm = interp1(linspace(0,1,5), ...
    [C.blue;[0.3,0.6,0.9];[0.97,0.97,0.97];[0.99,0.68,0.25];C.red], linspace(0,1,256));

set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',9);

mode_tex_bare = {'a_{0,k_f}','a_{1,0}','a_{1,k_f}','a_{2,0}', ...
                 'a_{0,2k_f}','a_{1,2k_f}','a_{2,k_f}'};
mode_tex      = {'$a_{0,k_f}$','$a_{1,0}$','$a_{1,k_f}$','$a_{2,0}$', ...
                 '$a_{0,2k_f}$','$a_{1,2k_f}$','$a_{2,k_f}$'};

%% ---- Physical grid & data prep ------------------------------------------
[~,max_l51]=max(Sigma_all(:,i_large));
xi_51=edges(max_l51,1); yi_51=edges(max_l51,2);

nDFT_half=nDFT_s/2; t_vec_sp=t_vec(nDFT_half+1:end); nt_sp=length(t_vec_sp);
Z_sp=nodes_spod(:,i_large); X_sp=nodes_spod(:,xi_51); Y_sp=nodes_spod(:,yi_51);

[~,t_hi_sp]=max(movmean(Z_sp,20));
Z_mask=movmean(Z_sp,20);
Z_mask(max(1,t_hi_sp-round(100/dt_save)):min(nt_sp,t_hi_sp+round(100/dt_save)))=NaN;
[~,t_lo_sp]=min(Z_mask);
t_hi=t_hi_sp+nDFT_half; t_lo=t_lo_sp+nDFT_half;

hw_r=round(200/dt_save);
w_lo_r=max(1,t_lo_sp-hw_r):min(nt_sp,t_lo_sp+hw_r);
w_hi_r=max(1,t_hi_sp-hw_r):min(nt_sp,t_hi_sp+hw_r);
w_roll=round(25/dt_save); r_roll=NaN(nt_sp,1);
for ti=1:nt_sp
    i1=max(1,ti-w_roll); i2=min(nt_sp,ti+w_roll);
    if i2-i1>=10, r_roll(ti)=corr(X_sp(i1:i2),Y_sp(i1:i2)); end
end
r_at_lo=r_roll(t_lo_sp); r_at_hi=r_roll(t_hi_sp);
fprintf('Rolling r:  peak=%+.4f  trough=%+.4f\n',r_at_hi,r_at_lo);

Z_ylim=[0, max(Z_sp)*1.15];

%% Figure B
fig=figure('Units','centimeters','Position',[2 2 17.4 17/2],'Color','w');
tl=tiledlayout(2,4,'TileSpacing','compact','Padding','tight');

x_tick=[0,pi/2,pi,3*pi/2,2*pi];
x_tlab={'$0$','$\pi/2$','$\pi$','$3\pi/2$','$2\pi$'};

%% ---- (a) Heatmap  [rows1-2, cols2-4] -----------------------------------
ax_c=nexttile(tl,1,[2,2]);

[~,sort_l]=sort(max(Theta_full_kol,[],2,'omitnan'),'descend');
n_show=min(15,size(edges,1)); idx_show=sort_l(1:n_show);
Th_show=Theta_full_kol(idx_show,:);
Th_max_c=max(Th_show(~isnan(Th_show)));
Th_img=Th_show; Th_img(isnan(Th_img)|Th_img<0)=0;
imagesc(ax_c,Th_img); colormap(ax_c,parula); clim(ax_c,[0,Th_max_c]);
hold(ax_c,'on');
for l=1:n_show
    for zi=1:n_modes
        if isnan(Th_show(l,zi))
            fill(ax_c,[zi-0.5,zi+0.5,zi+0.5,zi-0.5,zi-0.5], ...
                      [l-0.5,l-0.5,l+0.5,l+0.5,l-0.5],[0.82 0.82 0.82],'EdgeColor','none');
            continue;
        end
        tc=ternary_local(Th_show(l,zi)>0.6*Th_max_c,'k','w');
        text(ax_c,zi,l,sprintf('%.0f',Th_show(l,zi)), ...
            'HorizontalAlignment','center','VerticalAlignment','middle', ...
            'FontSize',7,'Color',tc,'FontWeight','bold','Interpreter','none');
    end
end
rectangle(ax_c,'Position',[i_large-0.5,0.5,1,n_show], ...
    'EdgeColor','w','LineWidth',2,'FaceColor','none');
hold(ax_c,'off');
pair_lbl=cell(n_show,1);
for l=1:n_show
    pair_lbl{l}=sprintf('$(%s,\\,%s)$', ...
        mode_tex_bare{edges(idx_show(l),1)},mode_tex_bare{edges(idx_show(l),2)});
end
set(ax_c,'XTick',1:n_modes,'XTickLabel',mode_tex,'TickLabelInterpreter','latex', ...
         'YTick',1:n_show,'YTickLabel',pair_lbl,'TickLabelInterpreter','latex','FontSize',8);
xtickangle(ax_c,30);
xlabel(ax_c,'Regulator $Z$','Interpreter','latex');
ylabel(ax_c,'Coupled pair $(X,Y)$','Interpreter','latex');
title(ax_c,'(a) $\Theta_\Sigma$ per triad','Interpreter','latex','FontSize',9.5);
box(ax_c,'on');

%% ---- (b) |a(1,0)| - trough window  [row3, cols1-2] --------------------
ax_d=nexttile(tl,3,[1,2]);
plot(ax_d,t_vec(w_lo_r),Z_sp(w_lo_r),'-','Color',C.green,'LineWidth',1.8);hold on;
xline(0,'--','Color',[0.35 0.35 0.35],'LineWidth',1,'HandleVisibility','off');
title(ax_d,'(d) Trough','Interpreter','Latex','FontSize',9.5);
yline(0,'--','Color',C.grey,'LineWidth',0.9,'HandleVisibility','off');
plot(ax_d,t_vec(w_lo_r),r_roll(w_lo_r),'-','Color',C.blue,'LineWidth',2.0);
xline(t_vec(2508),':','Color',C.blue,'LineWidth',2,'HandleVisibility','off');
text(ax_d,0.57,0.12,sprintf('$r(t_c)=%+.2f$',r_at_lo), ...
    'Units','norm','Interpreter','latex','FontSize',9,'HorizontalAlignment','left', ...
    'Color',C.blue,'FontWeight','bold');
hold off;
xlabel(ax_d,'$t$','Interpreter','latex');
ylabel(ax_d,...
['\color[rgb]{0.000 0.447 0.698} r(X,Y)' ...
 ' \color[rgb]{0.000 0.620 0.451} a(1,0)']);
title(ax_d,'(b)','Interpreter','Latex','FontSize',9.5);
ylim(ax_d,[-1 1]);
xlim(ax_d,[1200 1300]);
box on; grid on;

%% ---- (c) |a(1,0)| - peak window  [row3, cols3-4] ----------------------
ax_e=nexttile(tl,7,[1,2]);
plot(ax_e,t_vec(w_hi_r),Z_sp(w_hi_r),'-','Color',C.green,'LineWidth',1.8);hold on;
xline(0,'--','Color',[0.35 0.35 0.35],'LineWidth',1,'HandleVisibility','off');
title(ax_e,'(e) Peak','Interpreter','none','FontSize',9.5);
yline(0,'--','Color',C.grey,'LineWidth',0.9,'HandleVisibility','off');
plot(ax_e,t_vec(w_hi_r),r_roll(w_hi_r),'-','Color',C.red,'LineWidth',2.0);
xline(t_vec(4992),':','Color',C.red,...
      'LineWidth',2,'HandleVisibility','off');
text(ax_e,0.97,0.12,sprintf('$r(t_c)=%+.2f$',r_roll(5000)), ...
    'Units','norm','Interpreter','latex','FontSize',9,'HorizontalAlignment','right', ...
    'Color',C.red,'FontWeight','bold');
hold off;
xlabel(ax_e,'$t$','Interpreter','latex');
ylabel(ax_e,...
['\color[rgb]{0.835 0.369 0.000} r(X,Y)' ...
 ' \color[rgb]{0.000 0.620 0.451} a(1,0)']);
title(ax_e,'(c)','Interpreter','Latex','FontSize',9.5);
ylim(ax_e,[-1 1]);
xlim(ax_e,[2450 2550]);
box on; grid on;


%% ---- Save ---------------------------------------------------------------
out_dir = fullfile(fileparts(mfilename('fullpath')), '..', 'Plots');
if ~exist(out_dir,'dir'), mkdir(out_dir); end
print(fig, fullfile(out_dir,'fig_kolmogorov_B.pdf'), '-dpdf', '-vector', '-r300');
print(fig, fullfile(out_dir,'fig_kolmogorov_B.png'), '-dpng', '-r300');
fprintf('Saved: Plots/fig_kolmogorov_B\n');

%% ---- Local helpers ------------------------------------------------------
function s=ternary_local(c,a,b), if c, s=a; else, s=b; end, end
