# Weighted LLH and EDR, vector normalization, and SP weights (Sections 2.1, 4; Eqs. 1-7; Tables 3-4).

source("config.R", local = TRUE)
source("R/utils.R", local = TRUE)

setup   <- tectonic_setup(TECTONIC)
tag     <- setup$tag
out_dir <- file.path(RESULTS_DIR, paste0("scores_", tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ds <- load_scoring_data(setup)
if (nrow(ds$dropped) > 0)
  write_csv(ds$dropped, file.path(out_dir, "dropped_records.csv"))

sc     <- score_matrices(ds)
models <- rownames(sc$LLH)

write_csv(sc$cells %>%
            select(EQ_ID, Sta, Mw, Rrup, x_bin, y_bin, n_b, B, weight) %>%
            mutate(across(c(x_bin, y_bin), as.character)),
          file.path(out_dir, paste0("cell_weights_", tag, ".csv")))

write_csv(data.frame(n_total = ds$n_total,
                     n_used = nrow(ds$dat),
                     n_dropped = ds$n_total - nrow(ds$dat),
                     n_events_used = length(unique(ds$dat$EQ_ID)),
                     n_cells = unique(sc$cells$B),
                     INCLUDE_LOG_SIGMA = INCLUDE_LOG_SIGMA),
          file.path(out_dir, "run_info.csv"))

write_csv(as_period_table(sc$LLH, models), file.path(out_dir, paste0("LLH_weighted_", tag, ".csv")))
write_csv(as_period_table(sc$EDR, models), file.path(out_dir, paste0("EDR_weighted_", tag, ".csv")))
write_csv(as_period_table(sc$sp$LLHs_norm, models), file.path(out_dir, paste0("LLH_vecnorm_", tag, ".csv")))
write_csv(as_period_table(sc$sp$EDRs_norm, models), file.path(out_dir, paste0("EDR_vecnorm_", tag, ".csv")))

sp_out <- data.frame(model = models, spj = sc$sp$spj[, 1], SP = sc$sp$SP[, 1])
sp_out <- sp_out[order(-sp_out$SP), ]
write_csv(sp_out, file.path(out_dir, paste0("SP_vecnorm_", tag, ".csv")))

message("Results written to ", out_dir)
print(sp_out, row.names = FALSE)
