if (!exists("RAW_FLATFILE")) RAW_FLATFILE <- "data/raw/1994_2024_flatfile_filled4.csv"
if (!exists("DATA_FILE"))    DATA_FILE    <- "data/TSMIP_flatfile_1992_2024.csv"
if (!exists("YEAR_MIN"))     YEAR_MIN     <- 1992
if (!exists("YEAR_MAX"))     YEAR_MAX     <- 2024
if (!exists("TECTONIC"))     TECTONIC     <- "cru"
RUN_TAG <- paste0(YEAR_MIN, "_", YEAR_MAX)

INTER_DIR   <- file.path("intermediate", RUN_TAG)
RESULTS_DIR <- file.path("results", RUN_TAG)
FIG_DIR     <- file.path("figures", RUN_TAG)

PERIODS  <- c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12,
              0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1,
              1.5, 2, 3, 4, 5, 7.5, 10)
OBS_COLS <- paste0("T", PERIODS, "S")
N_PER    <- length(PERIODS)

DIST_BREAKS <- 10^(seq(log10(0.01), log10(1000), by = 0.1))

INCLUDE_LOG_SIGMA <- FALSE
MAG_COL <- "Mag"
LN_G    <- 6.89
