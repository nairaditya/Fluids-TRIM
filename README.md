# TRIM: Triadic Interaction Mining in fluid flows

MATLAB figure-generation code for the TRIM framework, which detects and
quantifies amplitude-mediated triadic regulation in oscillatory and
turbulent flows. Given a coupled mode pair `(X, Y)` and a candidate
regulator mode `Z`, TRIM measures how the conditional coupling
`MI(X, Y | |Z|)` varies with the regulator amplitude and tests the
result against an amplitude-shuffled null, summarised by the regulatory
gain statistic `Theta_Sigma`.

Each script regenerates one figure from precomputed results stored in the
accompanying `.mat` files. The analysis itself (SPOD/BMD decomposition,
mutual-information estimation, null-distribution bootstrap) is run upstream;
only the precomputed results needed to reproduce the figures are included
here.

## Repository layout

```
.
|-- Codes/
|   |-- run_all_figures.m            % driver: regenerates every figure
|   |-- fig_triadic_regulation.m     % synthetic INC/DEC schematic
|   |-- fig_stuart_landau.m          % Stuart-Landau validation
|   |-- fig_kolmogorov_A.m           % Kolmogorov: detection + identification
|   |-- fig_kolmogorov_B.m           % Kolmogorov: heatmap + flow-field gating
|   |-- fig_jet.m                    % turbulent jet: combined analysis figure
|   |-- results_triadic.mat
|   |-- results_sl.mat
|   |-- results_kolmogorov.mat
|   |-- results_jet.mat
|-- Plots/                           % created on first run (figure output)
|-- README.md
```

## Requirements

- MATLAB R2022a or newer (the scripts use `tiledlayout`/`nexttile` and `clim`).
- Statistics and Machine Learning Toolbox (`corr`, `prctile`).

## Usage

From the `Codes/` directory, regenerate all figures:

```matlab
run_all_figures
```

Or run any figure script on its own, for example:

```matlab
fig_jet
```

Each script loads its `results_*.mat` file, builds the figure, and writes
both PDF (vector) and PNG (300 dpi) into a `Plots/` directory created one
level above `Codes/`.

## Scripts and outputs

| Script                     | Input data               | Output (in `Plots/`)            |
| -------------------------- | ------------------------ | ------------------------------- |
| `fig_triadic_regulation.m` | `results_triadic.mat`    | `fig_triadic_regulation.pdf/png`|
| `fig_stuart_landau.m`      | `results_sl.mat`         | `fig_SL.pdf/png`                |
| `fig_kolmogorov_A.m`       | `results_kolmogorov.mat` | `fig_kolmogorov_A.pdf/png`      |
| `fig_kolmogorov_B.m`       | `results_kolmogorov.mat` | `fig_kolmogorov_B.pdf/png`      |
| `fig_jet.m`                | `results_jet.mat`        | `fig_jet_combined.pdf/png`      |

## Test cases

- **Synthetic triadic system** (`fig_triadic_regulation`): a controlled
  example with a known regulator, used to illustrate the INC (amplifying)
  and DEC (suppressing) regulation fingerprints.
- **Stuart-Landau oscillator** (`fig_stuart_landau`): validation against a
  reference oscillator, including a null test that a genuine regulator
  rejects `H_0` while correlated and uncorrelated controls do not.
- **Kolmogorov flow** (`fig_kolmogorov_A`, `fig_kolmogorov_B`): regulatory
  gain spectrum across SPOD modes, bootstrap rank stability, the
  per-triad `Theta_Sigma` heatmap, and flow-field gating during peak and
  trough regulator windows.
- **Turbulent jet** (`fig_jet`): SPOD eigenspectrum, BMD bispectrum, the
  per-triad TRIM heatmap, and conditional-coupling profiles for the two
  leading triads.
