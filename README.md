# h2o2comp

Code and data for the H2O2-mediated competition between *Prochlorococcus* and *Synechococcus*.

## Repository layout

```
data/      raw and derived data
src/       analysis and figure-generation scripts
figures/   output figures
```

## Requirements

Python 3.10+ with the following packages:

```
numpy
pandas
scipy
matplotlib
statsmodels
openpyxl
```

Install with:

```
pip install numpy pandas scipy matplotlib statsmodels openpyxl
```

## Reproducing all figures

From the `src/` directory:

```
bash run_all.sh
```

This runs every script in order and writes outputs to `../figures/` and `../data/`.

## Running scripts individually

From the `src/` directory:

| Script | Outputs |
|---|---|
| `intro_plot.py`   | `figures/intro_plot.png` (figure 1) |
| `fit_kdam.py`     | `figures/figureS1.png`, `figureS2.png`, `figureS3.png`; tables in `data/` |
| `simulations.py`  | `figures/figure2.png` through `figure6.png` |
| `sensitivity.py`  | `figures/figureS4.png` |

`functions.py` holds the model parameters and the leak/equilibrium functions; it is imported by `simulations.py` and `sensitivity.py` rather than run directly.

## Data files

| File | Source |
|---|---|
| `Ma_fig_1_raw_data_compiled.xlsx` | Growth curves used to fit kdam |
| `data_MEGA.xlsx`                  | Field cell counts and H2O2 used by `intro_plot.py` |
| `kdam_estimates.csv`              | Full kdam fit output (mu, slope, SE, CI, R2) |
| `kdam_table_clean.csv`            | kdam point estimate, 95% CI, n, R2 |
| `mu_from_zero_h2o2_controls.csv`  | Per-strain, per-temperature growth rates from controls |
| `tidy_growth_data.csv`            | Long-form replicate-level growth data |
| `prochlorococcus_kdam_outputs.xlsx` | All of the above bundled as Excel sheets |
