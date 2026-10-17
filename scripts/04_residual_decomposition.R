# Mixed-effects decomposition of the residuals at PGA, 0.17 s and 1.0 s (Section 6.2).

source("config.R", local = TRUE)
source("R/utils.R", local = TRUE)
source("R/gmpe_models.R", local = TRUE)
suppressPackageStartupMessages(library(lme4))

setup <- tectonic_setup(TECTONIC)
tag   <- setup$tag

DECOMP_PERIODS <- c(0, 0.17, 1.0)
PER_TAG        <- c("PGA", "T0.17s", "T1.0s")
PER_IM         <- c("PGA", "SA(0.17s)", "SA(1.0s)")
EVENT_ID   <- "EQ_ID"
STATION_ID <- "Sta"

dir_rds    <- file.path(INTER_DIR, "calib_resid")
dir_decomp <- file.path(RESULTS_DIR, paste0("residual_decomposition_", tag))
dir.create(dir_rds,    recursive = TRUE, showWarnings = FALSE)
dir.create(dir_decomp, recursive = TRUE, showWarnings = FALSE)

regre <- load_flatfile(setup$eq_types)

decompose_residuals <- function(df) {
  df <- df[is.finite(df$residual), ]
  df[[EVENT_ID]]   <- factor(df[[EVENT_ID]])
  df[[STATION_ID]] <- factor(df[[STATION_ID]])

  m1 <- lmer(as.formula(paste0("residual ~ 1 + (1|", EVENT_ID, ")")), data = df, REML = TRUE)
  c1 <- fixef(m1)[["(Intercept)"]]
  re1 <- ranef(m1)[[EVENT_ID]]
  dB_tab <- setNames(data.frame(rownames(re1), re1[, 1]), c(EVENT_ID, "dB"))
  df$dB <- dB_tab$dB[match(df[[EVENT_ID]], dB_tab[[EVENT_ID]])]
  df$dW <- df$residual - c1 - df$dB

  m2 <- lmer(as.formula(paste0("dW ~ 1 + (1|", STATION_ID, ")")), data = df, REML = TRUE)
  c2 <- fixef(m2)[["(Intercept)"]]
  re2 <- ranef(m2)[[STATION_ID]]
  dS_tab <- setNames(data.frame(rownames(re2), re2[, 1]), c(STATION_ID, "dS2S"))
  df$dS2S <- dS_tab$dS2S[match(df[[STATION_ID]], dS_tab[[STATION_ID]])]
  df$dWS  <- df$dW - c2 - df$dS2S

  stats <- data.frame(n_records = nrow(df), n_events = nlevels(df[[EVENT_ID]]),
                      n_stations = nlevels(df[[STATION_ID]]), bias = c1 + c2,
                      tau = as.data.frame(VarCorr(m1))$sdcor[1],
                      phi_S2S = as.data.frame(VarCorr(m2))$sdcor[1],
                      phi_SS = sd(df$dWS))
  list(data = df, event_table = dB_tab, station_table = dS_tab, stats = stats)
}

summary_rows <- list()
for (lab in names(setup$models)) {
  model_name <- paste0(lab, "(", setup$suffix, ")")
  dist_col   <- if (lab %in% names(gmpe_distance)) gmpe_distance[[lab]] else "Rrup"
  pred_all   <- read_model_table("predict", setup$models[[lab]], setup$suffix)
  if (nrow(pred_all) != nrow(regre)) stop(model_name, ": prediction file does not match the flat file.")

  for (i in seq_along(DECOMP_PERIODS)) {
    obs_col  <- paste0("T", DECOMP_PERIODS[i], "S")
    file_tag <- paste0(model_name, "_", PER_TAG[i])
    regre$pred <- pred_all[[match(DECOMP_PERIODS[i], PERIODS)]]

    d <- regre[is.finite(regre[[obs_col]]) & is.finite(regre$pred) & regre[[obs_col]] > 0 & regre$pred > 0, ]
    d$residual <- log10(d[[obs_col]]) - log10(d$pred)
    decomp <- decompose_residuals(d)
    rec <- decomp$data

    vs30_sta <- rec %>%
      group_by(.data[[STATION_ID]]) %>%
      summarise(Vs30 = median(Vs30), .groups = "drop")
    station_df <- decomp$station_table %>% left_join(vs30_sta, by = STATION_ID)
    event_df   <- rec %>% distinct(.data[[EVENT_ID]], M, dB) %>% rename(Mw = M)

    saveRDS(list(model_name = model_name, per_tag = PER_TAG[i], im_label = PER_IM[i],
                 dist_col = dist_col, dist_label = dist_col,
                 calib   = data.frame(obs = rec[[obs_col]], pred = rec$pred),
                 record  = data.frame(dist = rec[[dist_col]], dWS = rec$dWS),
                 event   = as.data.frame(event_df),
                 station = as.data.frame(station_df)),
            file.path(dir_rds, paste0(file_tag, ".rds")))
    write.csv(rec[, c(EVENT_ID, STATION_ID, "Mw", "M", "Rrup", "Rjb", "Vs30", obs_col, "pred",
                      "residual", "dB", "dS2S", "dWS")],
              file.path(dir_decomp, paste0(file_tag, "_record_level_residuals.csv")), row.names = FALSE)
    write.csv(decomp$event_table, file.path(dir_decomp, paste0(file_tag, "_event_level_dB.csv")), row.names = FALSE)

    summary_rows[[length(summary_rows) + 1]] <- cbind(model = model_name, period = PER_TAG[i], decomp$stats)
    message(file_tag, " done")
  }
}

summary_tbl <- do.call(rbind, summary_rows)
write.csv(summary_tbl, file.path(dir_decomp, paste0("residual_components_", tag, ".csv")), row.names = FALSE)
print(summary_tbl, row.names = FALSE)
