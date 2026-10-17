# Median predictions and ln residuals of the candidate GMPEs for every record and period.

source("config.R", local = TRUE)
source("R/utils.R", local = TRUE)
source("R/gmpe_models.R", local = TRUE)
suppressPackageStartupMessages(library(GMPEhaz))

setup  <- tectonic_setup(TECTONIC)
regre  <- load_flatfile(setup$eq_types)
models <- gmpe_models(setup$tag)

pred_dir  <- file.path(INTER_DIR, "predict")
resid_dir <- file.path(INTER_DIR, "residual")
dir.create(pred_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(resid_dir, recursive = TRUE, showWarnings = FALSE)

for (lab in names(models)) {
  per  <- PERIODS[seq_len(setup$n_periods[[lab]])]
  pred <- sapply(per, function(tp)
    exp(sapply(seq_len(nrow(regre)), function(k) unlist(models[[lab]](regre, k, tp))[["lnY"]])) / exp(LN_G))
  pred <- matrix(pred, nrow = nrow(regre))
  res  <- log(as.matrix(regre[, OBS_COLS[seq_along(per)]])) - log(pred)
  colnames(pred) <- paste0("predict_", per)
  colnames(res)  <- paste0("residual_", per)

  file_name <- paste0(setup$models[[lab]], "(", setup$suffix, ").csv")
  write_csv(as.data.frame(pred), file.path(pred_dir,  file_name))
  write_csv(as.data.frame(res),  file.path(resid_dir, file_name))
  message(lab, " done")
}
