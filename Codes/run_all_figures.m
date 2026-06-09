%% RUN_ALL_FIGURES.M  -  Generate all JFM figures.
%  -----------------
%    fig_triadic_regulation  - synthetic INC/DEC schematic
%    fig_SL                  - Stuart-Landau TRIM validation
%    fig_kolmogorov_A        - Kolmogorov: detection + regulatory identification
%    fig_kolmogorov_B        - Kolmogorov: heatmap + time-series gating
%    fig_jet_combined        - Turbulent jet: full analysis figure

here  = fileparts(mfilename('fullpath'));
plots = fullfile(here, '..', 'Plots');
if ~exist(plots,'dir'), mkdir(plots); end

script_names = { ...
    'fig_triadic_regulation.m', ...
    'fig_stuart_landau.m', ...
    'fig_kolmogorov_A.m', ...
    'fig_kolmogorov_B.m', ...
    'fig_jet.m'};

t_total = tic;
for idx = 1:numel(script_names)
    sname = script_names{idx};
    t0 = tic;
    run_isolated(fullfile(here, sname));
    close all;
end

function run_isolated(script_path)
    run(script_path);
end
