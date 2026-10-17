# Data-driven GMPE weighting for Taiwan

Code and data for *Data-Driven Weighting of Ground-Motion Prediction Equations for Seismic Hazard Analysis in Taiwan* (Gao & Yang, Engineering Geology).

## Requirements

R (>= 4.2) with `dplyr`, `readr`, `reshape2`, `ggplot2`, `tidyr`, `patchwork`, `lme4`, `MASS`, `ggtext`, `MethComp`, `scales`, and `GMPEhaz` (GMPE implementations, used in steps 01 and 08).

## Layout

```
config.R                    data file, year range, spectral periods, constants
run_all.R                   runs every step
R/utils.R                   record selection, magnitude-distance cell weights, scoring
R/gmpe_models.R             input parameters passed to each GMPE
R/llh_weighted.R            weighted log-likelihood
R/edr_weighted.R            weighted Euclidean distance-based ranking
R/sp_vecnorm.R              vector normalization and SP weights
data/TSMIP_flatfile_1992_2024.csv   records and parameters used in this study
intermediate/<years>/       predict/, residual/, sigma/, calib_resid/
results/<years>/            tables
figures/<years>/            figures
```

| Script | Output | Manuscript |
|---|---|---|
| `00_prepare_flatfile.R` | reduced flat file | Section 2 |
| `01_gmpe_predictions.R` | median predictions and ln residuals | Section 3 |
| `02_score_weights.R` | cell weights, weighted LLH/EDR, SP weights | Eqs. 1-7, Tables 3-4 |
| `03_plot_scores.R` | LLH, EDR, SP and weight figures | Figures 3-5 |
| `04_residual_decomposition.R` | between-event, site-to-station and within-event residuals | Section 6.2 |
| `05_plot_residuals.R` | calibration and residual figures | Figures 8, X5-X16 |
| `06_bootstrap.R` | event- and station-level cluster bootstrap | Table 5 |
| `07_bootstrap_figures.R` | bootstrap figures | Figure X18 |
| `08_sammon_map.R` | Sammon's maps for the crustal scenarios | Figure 11a-b |

## Flat file

One row per record (64,401 records from 697 earthquakes, 1992-2024). Columns:

| Column | Description |
|---|---|
| `EQ_ID`, `Sta` | event and station identifiers |
| `EQ_type` | shallow crustal, subduction interface or subduction intraslab |
| `Mw` | moment magnitude |
| `Mag` | magnitude passed to the GMPEs (`MAG_COL` in `config.R`) |
| `depth`, `Ztor`, `Dip`, `width`, `ftype` | hypocentral depth (km), depth to top of rupture (km), dip (deg), rupture width (km), mechanism flag |
| `Rrup`, `Rjb`, `Rx`, `Ry` | distance metrics (km) |
| `Vs30`, `Vs30_class`, `Z1.0`, `Z2.5` | site parameters (m/s, flag, m, km) |
| `Dmoho`, `hwflag`, `basinflag`, `sourcetype`, `msasflag` | additional GMPE inputs |
| `T0S` ... `T10S` | RotD50 PGA and 5%-damped PSA at 24 periods (g) |

## Intermediate files

`intermediate/<years>/{predict,residual,sigma}/<model>(CRU|SUB).csv` hold one row per record of the
selected flat-file subset and one column per spectral period (25 for global GMPEs, 23 for
Taiwan-specific GMPEs). `predict` and `residual` are written by `01_gmpe_predictions.R`; `sigma` contains
the total standard deviation (ln units) of each GMPE for each record and period.

## Usage

Run from the repository root:

```r
source("run_all.R")
```

or a single step:

```r
TECTONIC <- "sub"; YEAR_MIN <- 2017; YEAR_MAX <- 2024
source("scripts/02_score_weights.R")
```

## Data sources

TSMIP records: Central Weather Administration, https://doi.org/10.7914/SN/T5.
Flat file for 1992-2014: Taiwan SSHAC Level 3 PSHA Project, https://sshac.ncree.org.tw.
