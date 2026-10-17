# Event- and station-level cluster bootstrap of the GMPE weights (Section 6.3; Table 5).

source("config.R", local = TRUE)
source("R/utils.R", local = TRUE)
suppressPackageStartupMessages(library(ggplot2))

if (!exists("N_BOOT"))       N_BOOT       <- 40
if (!exists("SEED_EVENT"))   SEED_EVENT   <- 4321
if (!exists("SEED_STATION")) SEED_STATION <- 8765

setup      <- tectonic_setup(TECTONIC)
tag        <- setup$tag
scores_dir <- file.path(RESULTS_DIR, paste0("scores_", tag))
out_dir    <- file.path(RESULTS_DIR, paste0("bootstrap_", tag))
part_dir   <- file.path(out_dir, "parts")
fig_dir    <- file.path(FIG_DIR, paste0("bootstrap_", tag))
dir.create(part_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)

ds <- load_scoring_data(setup)

score_fun <- function(idx) {
  sc <- score_matrices(ds, idx)
  data.frame(model    = rownames(sc$LLH),
             LLH_norm = rowMeans(sc$sp$LLHs_norm, na.rm = TRUE),
             EDR_norm = rowMeans(sc$sp$EDRs_norm, na.rm = TRUE),
             spj      = sc$sp$spj[, 1],
             SP       = sc$sp$SP[, 1],
             row.names = NULL)
}

full <- score_fun(seq_len(nrow(ds$dat)))
ref_file <- file.path(scores_dir, paste0("SP_vecnorm_", tag, ".csv"))
if (file.exists(ref_file)) {
  ref <- read_csv_clean(ref_file)
  cmp <- merge(full[, c("model", "SP")], ref[, c("model", "SP")], by = "model")
  if (nrow(cmp) != nrow(full) || max(abs(cmp$SP.x - cmp$SP.y)) > 1e-6)
    stop("Full-sample weights differ from ", ref_file, "; rerun 02_score_weights.R first.")
}

run_cluster_bootstrap <- function(df, cluster_id, n_boot, seed, prefix) {
  clusters <- unique(df[[cluster_id]])
  idx_by_cluster <- split(seq_len(nrow(df)), df[[cluster_id]])
  results <- vector("list", n_boot)
  for (b in seq_len(n_boot)) {
    f <- file.path(part_dir, sprintf("%s_%04d.csv", prefix, b))
    if (file.exists(f)) {
      results[[b]] <- read.csv(f, stringsAsFactors = FALSE)
      next
    }
    set.seed(seed + b)
    sampled <- sample(clusters, size = length(clusters), replace = TRUE)
    idx <- unlist(idx_by_cluster[as.character(sampled)], use.names = FALSE)
    scores_b <- score_fun(idx)
    scores_b$boot_id <- b
    write.csv(scores_b, f, row.names = FALSE)
    results[[b]] <- scores_b
    message(sprintf("%s-level bootstrap: %d / %d", prefix, b, n_boot))
  }
  bind_rows(results)
}

summarize_bootstrap <- function(boot, col, margin = 0.10) {
  tbl <- boot %>%
    group_by(boot_id) %>%
    mutate(rank = rank(-.data[[col]], ties.method = "average"),
           within_margin = .data[[col]] >= (1 - margin) * max(.data[[col]], na.rm = TRUE)) %>%
    ungroup() %>%
    group_by(model) %>%
    summarise(median_weight   = median(.data[[col]], na.rm = TRUE),
              ci_lower_2.5    = quantile(.data[[col]], 0.025, na.rm = TRUE),
              ci_upper_97.5   = quantile(.data[[col]], 0.975, na.rm = TRUE),
              p_rank1         = mean(rank == 1),
              p_top3          = mean(rank <= 3),
              p_within_margin = mean(within_margin),
              .groups = "drop") %>%
    arrange(desc(median_weight))
  tbl %>% mutate(in_leading_group = ci_upper_97.5 >= tbl$ci_lower_2.5[1])
}

plot_bootstrap <- function(boot, col, title_txt, ylim) {
  lv <- boot %>% group_by(model) %>% summarise(med = median(.data[[col]], na.rm = TRUE)) %>%
    arrange(desc(med)) %>% pull(model)
  boot %>%
    mutate(model = factor(model, levels = lv)) %>%
    ggplot(aes(x = model, y = .data[[col]])) +
    geom_violin(fill = "grey85", color = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.5) +
    coord_cartesian(ylim = ylim) +
    labs(title = title_txt, x = NULL, y = col) +
    theme_test() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

boot_event <- run_cluster_bootstrap(ds$dat, "EQ_ID", N_BOOT, SEED_EVENT, "event")
write_csv(boot_event, file.path(out_dir, paste0("bootstrap_event_raw_", tag, ".csv")))
boot_station <- run_cluster_bootstrap(ds$dat, "Sta", N_BOOT, SEED_STATION, "station")
write_csv(boot_station, file.path(out_dir, paste0("bootstrap_station_raw_", tag, ".csv")))

summ <- list(event = list(), station = list())
for (col in c("LLH_norm", "EDR_norm", "spj", "SP")) {
  summ$event[[col]]   <- summarize_bootstrap(boot_event, col)
  summ$station[[col]] <- summarize_bootstrap(boot_station, col)
  write_csv(summ$event[[col]],   file.path(out_dir, paste0("bootstrap_event_",   col, "_summary_", tag, ".csv")))
  write_csv(summ$station[[col]], file.path(out_dir, paste0("bootstrap_station_", col, "_summary_", tag, ".csv")))

  rg <- range(c(boot_event[[col]], boot_station[[col]]), na.rm = TRUE)
  ylim <- rg + c(-1, 1) * diff(rg) * 0.05
  for (lvl in c("event", "station")) {
    boot <- if (lvl == "event") boot_event else boot_station
    ggsave(file.path(fig_dir, paste0("bootstrap_", lvl, "_", col, "_", tag, ".png")),
           plot_bootstrap(boot, col, sprintf("%s - %s - %s-level bootstrap (n = %d, %s)",
                                             setup$label, col, lvl, N_BOOT, RUN_TAG), ylim),
           width = 7, height = 4.5, dpi = 500)
  }
}

pick <- function(tbl, prefix)
  tbl %>% select(model, median_weight, ci_lower_2.5, ci_upper_97.5, p_rank1) %>%
  rename_with(~ paste0(prefix, .x), -model)
sp_compare <- pick(summ$event$SP, "event_") %>%
  left_join(pick(summ$station$SP, "station_"), by = "model") %>%
  left_join(select(full, model, SP_full = SP), by = "model") %>%
  arrange(desc(SP_full))
write_csv(sp_compare, file.path(out_dir, paste0("bootstrap_SP_event_vs_station_", tag, ".csv")))

write_csv(data.frame(tectonic = setup$label, RUN_TAG = RUN_TAG, n_records = nrow(ds$dat),
                     n_events = length(unique(ds$dat$EQ_ID)), n_stations = length(unique(ds$dat$Sta)),
                     n_boot = N_BOOT, seed_event = SEED_EVENT, seed_station = SEED_STATION,
                     INCLUDE_LOG_SIGMA = INCLUDE_LOG_SIGMA),
          file.path(out_dir, paste0("bootstrap_run_info_", tag, ".csv")))

print(sp_compare, row.names = FALSE)
