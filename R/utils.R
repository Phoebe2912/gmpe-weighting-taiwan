suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

source("R/llh_weighted.R", local = TRUE)
source("R/edr_weighted.R", local = TRUE)
source("R/sp_vecnorm.R",   local = TRUE)

FLATFILE_COLS <- c("EQ_ID", "Sta", "EQ_type", "Mw", "Mag", "depth", "Ztor", "Dip", "width", "ftype",
                   "Rrup", "Rjb", "Rx", "Ry", "Vs30", "Vs30_class", "Z1.0", "Z2.5", "Dmoho",
                   "hwflag", "basinflag", "sourcetype", "msasflag")
ANALYSIS_TYPES <- c("shallow crustal", "subduction interface", "subduction intraslab")

tectonic_setup <- function(tectonic) {
  stopifnot(tectonic %in% c("cru", "sub"))
  if (tectonic == "cru") {
    list(
      tag = "cru", suffix = "CRU", label = "Crustal",
      eq_types   = "shallow crustal",
      mag_breaks = seq(1.5, 8.5, 0.25),
      models = c(ASK14 = "ASK14", BSSA14 = "BSSA14", CB14 = "CB14", CY14 = "CY14", I14 = "I14",
                 PhungCru20 = "PhungCru20", Chao20 = "Chao20", Lin09 = "Lin09", Lin11 = "Lin11"),
      n_periods = c(ASK14 = 25, BSSA14 = 25, CB14 = 25, CY14 = 25, I14 = 25,
                    PhungCru20 = 23, Chao20 = 23, Lin09 = 23, Lin11 = 23),
      plot_order = c("Chao20", "PhungCru20", "Lin09", "Lin11",
                     "ASK14", "BSSA14", "CB14", "CY14", "I14"),
      plot_cols  = c("brown", "orangered", "gold", "yellowgreen",
                     "green4", "dodgerblue", "blue3", "purple4", "pink2")
    )
  } else {
    list(
      tag = "sub", suffix = "SUB", label = "Subduction",
      eq_types   = c("subduction interface", "subduction intraslab"),
      mag_breaks = seq(1.5, 9.5, 0.25),
      models = c(AG2020 = "AG2022", KBCG2020 = "KBCG2020", PSHAB2022 = "PSHAB2020", SI2022 = "SI2020",
                 Chao20 = "Chao20", PhungSub20 = "PhungSub20", LL08 = "LL08"),
      n_periods = c(AG2020 = 25, KBCG2020 = 25, PSHAB2022 = 25, SI2022 = 25,
                    Chao20 = 23, PhungSub20 = 23, LL08 = 23),
      plot_order = c("Chao20", "PhungSub20", "LL08",
                     "AG2020", "KBCG2020", "PSHAB2022", "SI2022"),
      plot_cols  = c("brown", "orangered", "gold",
                     "green4", "dodgerblue", "blue3", "purple4")
    )
  }
}

select_records <- function(ff, eq_types) {
  yr <- suppressWarnings(as.integer(substr(ff$EQ_ID, 1, 4)))
  ff <- ff[!is.na(yr) & yr >= YEAR_MIN & yr <= YEAR_MAX, ]
  keep <- ff$Mw >= 4
  for (v in c("Ztor", "Rjb", "Rx", "Ry", "Rrup", "Dip", "width", "ftype"))
    keep <- keep & !is.na(ff[[v]]) & ff[[v]] != -999
  ff <- ff[which(keep & ff$EQ_type %in% eq_types & !is.na(ff$T0S)), ]
  rownames(ff) <- NULL
  ff
}

load_flatfile <- function(eq_types) {
  ff <- read.csv(DATA_FILE, stringsAsFactors = FALSE, check.names = FALSE)
  stopifnot(all(c(FLATFILE_COLS, OBS_COLS) %in% names(ff)))
  ff <- select_records(ff, eq_types)
  ff$M <- ff[[MAG_COL]]
  ff
}

compute_cell_weights <- function(df, xvar = "Rrup", yvar = "Mw", xbreak, ybreak) {
  df <- df %>%
    mutate(x_bin = cut(.data[[xvar]], breaks = xbreak),
           y_bin = cut(.data[[yvar]], breaks = ybreak))
  cell_count <- df %>%
    filter(!is.na(x_bin), !is.na(y_bin)) %>%
    count(x_bin, y_bin, name = "n_b")
  B <- nrow(cell_count)
  cell_count <- cell_count %>% mutate(weight = 1 / (B * n_b))
  df %>%
    filter(!is.na(x_bin), !is.na(y_bin)) %>%
    left_join(cell_count, by = c("x_bin", "y_bin")) %>%
    mutate(B = B)
}

read_model_table <- function(kind, key, suffix) {
  f <- file.path(INTER_DIR, kind, paste0(key, c(paste0("(", suffix, ").csv"), paste0(suffix, ".csv"),
                                                 paste0("_", suffix, ".csv"))))
  f <- f[file.exists(f)]
  if (length(f) == 0) stop("Missing input file: ", file.path(INTER_DIR, kind, paste0(key, "(", suffix, ").csv")))
  as.data.frame(read_csv(f[1], show_col_types = FALSE))
}

load_scoring_data <- function(setup) {
  dat  <- load_flatfile(setup$eq_types)
  labs <- names(setup$models)
  read_all <- function(kind)
    setNames(lapply(labs, function(m) read_model_table(kind, setup$models[[m]], setup$suffix)), labs)
  resid <- read_all("residual")
  pred  <- read_all("predict")
  sigma <- read_all("sigma")

  for (m in labs) {
    for (kind in c("residual", "predict", "sigma")) {
      d <- switch(kind, residual = resid[[m]], predict = pred[[m]], sigma = sigma[[m]])
      if (nrow(d) != nrow(dat))
        stop(m, " ", kind, ": ", nrow(d), " rows, flatfile subset has ", nrow(dat), " rows.")
      if (ncol(d) < setup$n_periods[[m]])
        stop(m, " ", kind, ": fewer than ", setup$n_periods[[m]], " period columns.")
    }
  }

  valid  <- rep(TRUE, nrow(dat))
  reason <- rep("", nrow(dat))
  for (m in labs) {
    k <- seq_len(setup$n_periods[[m]])
    chk <- list(sigma = sigma[[m]][, k], pred = pred[[m]][, k],
                res = resid[[m]][, k], obs = dat[, OBS_COLS[k]])
    for (it in names(chk)) {
      x  <- suppressWarnings(matrix(as.numeric(as.matrix(chk[[it]])), nrow = nrow(dat)))
      ok <- is.finite(x)
      if (it != "res") ok <- ok & x > 0
      bad <- rowSums(!ok) > 0
      reason[bad] <- paste0(reason[bad], ifelse(reason[bad] == "", "", ";"), m, ":", it)
      valid <- valid & !bad
    }
  }

  dropped <- dat[!valid, c("EQ_ID", "Sta", "Mw", "Mag", "Vs30", "Rrup")]
  dropped$reason <- reason[!valid]
  sub_rows <- function(lst) lapply(lst, function(d) d[valid, , drop = FALSE])

  out_dat <- dat[valid, , drop = FALSE]
  rownames(out_dat) <- NULL
  message(sprintf("%s: %d records, %d dropped (%d events), %d retained.",
                  setup$label, nrow(dat), sum(!valid),
                  length(unique(dat$EQ_ID[!valid])), sum(valid)))

  list(setup = setup, dat = out_dat, n_total = nrow(dat), dropped = dropped,
       resid = sub_rows(resid), pred = sub_rows(pred), sigma = sub_rows(sigma))
}

score_matrices <- function(ds, idx = seq_len(nrow(ds$dat))) {
  setup <- ds$setup
  dat_s <- ds$dat[idx, , drop = FALSE]
  cells <- compute_cell_weights(dat_s, xbreak = DIST_BREAKS, ybreak = setup$mag_breaks)
  if (nrow(cells) != length(idx))
    stop("Some records fall outside the magnitude-distance grid.")
  w <- cells$weight

  labs <- names(setup$models)
  LLH  <- matrix(NA_real_, nrow = length(labs), ncol = N_PER,
                 dimnames = list(labs, paste0("V", seq_len(N_PER))))
  EDR  <- LLH

  for (m in labs) {
    res <- ds$resid[[m]][idx, , drop = FALSE]
    sig <- ds$sigma[[m]][idx, , drop = FALSE]
    prd <- ds$pred[[m]][idx, , drop = FALSE]
    for (j in seq_len(setup$n_periods[[m]])) {
      LLH[m, j] <- LLH_weighted(res[[j]], sigma = sig[[j]], w = w,
                                include_log_sigma = INCLUDE_LOG_SIGMA)
      obs_j  <- log(dat_s[[OBS_COLS[j]]])
      pred_j <- log(prd[[j]])
      sig_j  <- sig[[j]]
      ok <- is.finite(obs_j) & is.finite(pred_j) & is.finite(sig_j) & is.finite(w) & sig_j > 0
      EDR[m, j] <- EDR_weighted(obs_j[ok], pred_j[ok], sig_j[ok], w = w[ok])$EDR
    }
  }

  list(LLH = LLH, EDR = EDR, cells = cells, sp = SP_vecnorm_detail(LLH, EDR))
}

as_period_table <- function(mat, models) {
  df <- as.data.frame(mat)
  colnames(df) <- paste0("V", seq_len(ncol(df)))
  cbind(model = models, df, stringsAsFactors = FALSE)
}

read_csv_clean <- function(f) {
  d <- read_csv(f, show_col_types = FALSE, locale = locale(encoding = "UTF-8"))
  names(d) <- sub("^﻿", "", names(d))
  if ("model" %in% names(d)) d$model <- trimws(as.character(d$model))
  as.data.frame(d)
}
