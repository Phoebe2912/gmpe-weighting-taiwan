# Bootstrap distributions, confidence-interval widths and rank stability (Figure X18; Table 5).

source("config.R", local = TRUE)
source("R/utils.R", local = TRUE)
suppressPackageStartupMessages({
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

if (!exists("N_BOOT")) N_BOOT <- 40

setup      <- tectonic_setup(TECTONIC)
tag        <- setup$tag
boot_dir   <- file.path(RESULTS_DIR, paste0("bootstrap_", tag))
scores_dir <- file.path(RESULTS_DIR, paste0("scores_", tag))
out_dir    <- file.path(FIG_DIR, paste0("bootstrap_", tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

score_cols  <- c("LLH_norm", "EDR_norm", "SP")
score_label <- c(LLH_norm = "LLH (normalized)", EDR_norm = "EDR (normalized)", SP = "SP (final weight)")
level_cols  <- c("Event-level" = "#4C72B0", "Station-level" = "#DD8452")

boot_all <- bind_rows(
  read_csv_clean(file.path(boot_dir, paste0("bootstrap_event_raw_", tag, ".csv"))) %>%
    mutate(level = "Event-level"),
  read_csv_clean(file.path(boot_dir, paste0("bootstrap_station_raw_", tag, ".csv"))) %>%
    mutate(level = "Station-level")
) %>% mutate(level = factor(level, levels = names(level_cols)))

rg <- range(c(boot_all$LLH_norm, boot_all$EDR_norm), na.rm = TRUE)
fixed_ylim   <- c(floor(rg[1] / 0.05) * 0.05, ceiling(rg[2] / 0.05) * 0.05)
fixed_breaks <- seq(fixed_ylim[1], fixed_ylim[2], by = 0.05)

period_mean <- function(name) {
  d <- read_csv_clean(file.path(scores_dir, paste0(name, "_", tag, ".csv")))
  data.frame(model = d$model, value = rowMeans(d[, setdiff(names(d), "model")], na.rm = TRUE))
}
orig <- read_csv_clean(file.path(scores_dir, paste0("SP_vecnorm_", tag, ".csv"))) %>%
  select(model, spj, SP) %>%
  left_join(rename(period_mean("LLH_vecnorm"), LLH_norm = value), by = "model") %>%
  left_join(rename(period_mean("EDR_vecnorm"), EDR_norm = value), by = "model")

model_order    <- orig %>% arrange(desc(SP)) %>% pull(model)
boot_all$model <- factor(boot_all$model, levels = model_order)
orig$model     <- factor(orig$model, levels = model_order)

boot_long <- boot_all %>%
  select(model, level, boot_id, all_of(score_cols)) %>%
  pivot_longer(all_of(score_cols), names_to = "score", values_to = "value")
orig_long <- orig %>%
  select(model, all_of(score_cols)) %>%
  pivot_longer(all_of(score_cols), names_to = "score", values_to = "value")

caption_txt <- "Red cross: value computed from the complete dataset."

score_plot <- function(sc, standalone = TRUE) {
  gg <- ggplot(filter(boot_long, score == sc), aes(x = model, y = value, fill = level)) +
    geom_boxplot(width = 0.6, outlier.size = 0.8, position = position_dodge(width = 0.7)) +
    geom_point(data = filter(orig_long, score == sc), aes(x = model, y = value), inherit.aes = FALSE,
               shape = 4, size = 2.2, colour = "red", stroke = 1) +
    scale_fill_manual(values = level_cols) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top") +
    labs(x = "Model", y = score_label[[sc]], fill = NULL)
  if (sc != "SP")
    gg <- gg + scale_y_continuous(limits = fixed_ylim, breaks = fixed_breaks)
  if (standalone)
    gg <- gg + labs(title = paste0(score_label[[sc]], ": event-level vs. station-level bootstrap"),
                    subtitle = paste0("Number of replicates = ", N_BOOT), caption = caption_txt) +
      theme(plot.title = element_text(face = "bold"), plot.caption = element_text(hjust = 0, size = 8))
  gg
}

for (sc in score_cols)
  ggsave(file.path(out_dir, paste0("bootstrap_dist_", sc, ".png")), score_plot(sc), width = 8.5, height = 5, dpi = 300)

p1 <- (score_plot("LLH_norm", FALSE) + theme(legend.position = "none")) /
      (score_plot("EDR_norm", FALSE) + theme(legend.position = "none")) /
      score_plot("SP", FALSE) +
  plot_annotation(title = sprintf("%s: event-level vs. station-level cluster bootstrap (n = %d each), %s",
                                  setup$label, N_BOOT, RUN_TAG),
                  caption = caption_txt)
ggsave(file.path(out_dir, "bootstrap_dist_by_level.png"), p1, width = 8, height = 10, dpi = 300)

ci_width <- boot_all %>%
  group_by(model, level) %>%
  summarise(width = unname(quantile(SP, 0.975) - quantile(SP, 0.025)), .groups = "drop")
p2 <- ggplot(ci_width, aes(x = model, y = width, fill = level)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = level_cols) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Model", y = "95% interval width of SP", fill = NULL,
       title = sprintf("%s: event- vs. station-level bootstrap (n = %d), %s", setup$label, N_BOOT, RUN_TAG))
ggsave(file.path(out_dir, "bootstrap_ci_width_compare.png"), p2, width = 8, height = 4.5, dpi = 300)

rank_probs <- boot_all %>%
  group_by(level, boot_id) %>%
  mutate(rank = rank(-SP, ties.method = "average")) %>%
  ungroup() %>%
  count(level, model, rank) %>%
  group_by(level, model) %>%
  mutate(p = n / sum(n)) %>%
  ungroup() %>%
  complete(level, model, rank = seq_along(model_order), fill = list(p = 0))
p3 <- ggplot(rank_probs, aes(x = model, y = rank, fill = p)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(p > 0, sprintf("%.2f", p), "")), size = 2.8) +
  scale_y_reverse(breaks = seq_along(model_order)) +
  scale_fill_gradient(low = "white", high = "#4C72B0", limits = c(0, 1)) +
  facet_wrap(~ level) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold")) +
  labs(x = "Model", y = "Rank (1 = highest SP)", fill = "P(rank)",
       title = sprintf("%s: rank stability of SP weights (n = %d), %s", setup$label, N_BOOT, RUN_TAG))
ggsave(file.path(out_dir, "bootstrap_rank_heatmap.png"), p3, width = 9, height = 5, dpi = 300)

lower_is_better <- c("LLH_norm", "EDR_norm")
summary_all <- bind_rows(lapply(score_cols, function(col) {
  sgn <- if (col %in% lower_is_better) 1 else -1
  stats <- boot_all %>%
    group_by(level, model) %>%
    summarise(median = median(.data[[col]]),
              ci_lower_2.5  = quantile(.data[[col]], 0.025),
              ci_upper_97.5 = quantile(.data[[col]], 0.975),
              sd = sd(.data[[col]]), cv = sd / abs(mean(.data[[col]])), .groups = "drop")
  p1st <- boot_all %>%
    group_by(level, boot_id) %>%
    mutate(rank = rank(sgn * .data[[col]], ties.method = "average")) %>%
    group_by(level, model) %>%
    summarise(p_rank1 = mean(rank == 1), .groups = "drop")
  stats %>% left_join(p1st, by = c("level", "model")) %>% mutate(score = col)
})) %>%
  left_join(orig_long %>% rename(original = value), by = c("model", "score")) %>%
  arrange(score, level, model)
write_csv(summary_all, file.path(boot_dir, paste0("bootstrap_summary_", tag, ".csv")))

message("Figures written to ", out_dir)
