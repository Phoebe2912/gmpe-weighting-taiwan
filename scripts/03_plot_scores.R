# LLH, EDR and SP by spectral period, model averages, and final logic-tree weights (Figures 3-5).

source("config.R", local = TRUE)
source("R/utils.R", local = TRUE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(reshape2)
})

setup   <- tectonic_setup(TECTONIC)
tag     <- setup$tag
in_dir  <- file.path(RESULTS_DIR, paste0("scores_", tag))
out_dir <- file.path(FIG_DIR, paste0("scores_", tag))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

period_cols <- paste0("V", seq_len(N_PER))
period_map  <- setNames(PERIODS, period_cols)
model_order <- setup$plot_order
cols        <- setup$plot_cols

to_long <- function(df, value_name) {
  long <- melt(df, id = "model", variable.name = "col", value.name = value_name)
  long$Period <- period_map[as.character(long$col)]
  long$model  <- factor(long$model, levels = model_order)
  long
}

plot_line <- function(long_df, value_col, y_lab, file_name) {
  p <- ggplot(long_df) +
    scale_x_log10(breaks = PERIODS, labels = PERIODS) +
    geom_point(aes(x = Period, y = .data[[value_col]], shape = model, col = model), size = 3) +
    geom_line(aes(x = Period, y = .data[[value_col]], col = model), linewidth = 0.5) +
    scale_color_manual(values = cols) +
    scale_shape_manual(values = seq_along(model_order)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1)) +
    labs(y = y_lab, x = "Period (s)")
  ggsave(file.path(out_dir, file_name), plot = p, width = 8, height = 2.5)
}

plot_bar <- function(wide_df, y_lab, file_name) {
  avg_df <- data.frame(model = factor(wide_df$model, levels = model_order),
                       value = round(rowMeans(wide_df[, period_cols], na.rm = TRUE), 3))
  p <- ggplot(avg_df, aes(x = reorder(model, value), y = value)) +
    geom_bar(aes(fill = model, color = model), stat = "identity") +
    geom_text(aes(label = value, col = model),
              position = position_dodge(width = 1), vjust = -0.2, size = 3) +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
    labs(y = y_lab, x = "Model")
  ggsave(file.path(out_dir, file_name), plot = p, width = 5, height = 2.8)
}

read_scores <- function(name) read_csv_clean(file.path(in_dir, paste0(name, "_", tag, ".csv")))

llh_raw  <- read_scores("LLH_weighted")
edr_raw  <- read_scores("EDR_weighted")
llh_norm <- read_scores("LLH_vecnorm")
edr_norm <- read_scores("EDR_vecnorm")
sp_norm  <- cbind(model = llh_norm$model, llh_norm[, period_cols] * edr_norm[, period_cols])

plot_line(to_long(llh_raw,  "LLH"), "LLH", "LLH (weighted)",          paste0("llh_weighted_", tag, ".png"))
plot_bar (llh_raw,                  "Average LLH (weighted)",          paste0("llh_weighted_", tag, "_bar.png"))
plot_line(to_long(llh_norm, "LLH"), "LLH", "LLH (vector-normalized)", paste0("llh_vecnorm_", tag, ".png"))
plot_bar (llh_norm,                 "Average LLH (vector-normalized)", paste0("llh_vecnorm_", tag, "_bar.png"))
plot_line(to_long(edr_raw,  "EDR"), "EDR", "EDR (weighted)",          paste0("edr_weighted_", tag, ".png"))
plot_bar (edr_raw,                  "Average EDR (weighted)",          paste0("edr_weighted_", tag, "_bar.png"))
plot_line(to_long(edr_norm, "EDR"), "EDR", "EDR (vector-normalized)", paste0("edr_vecnorm_", tag, ".png"))
plot_bar (edr_norm,                 "Average EDR (vector-normalized)", paste0("edr_vecnorm_", tag, "_bar.png"))
plot_line(to_long(sp_norm,  "SP"),  "SP",  "SP (vector-normalized)",  paste0("sp_vecnorm_", tag, ".png"))
plot_bar (sp_norm,                  "Average SP (vector-normalized)",  paste0("sp_vecnorm_", tag, "_bar.png"))

sp_final <- read_scores("SP_vecnorm")
sp_final$model <- factor(sp_final$model, levels = model_order)
sp_final$SP    <- round(sp_final$SP, 3)
p_lt <- ggplot(sp_final, aes(x = reorder(model, -SP), y = SP)) +
  geom_bar(aes(fill = model, color = model), stat = "identity") +
  geom_text(aes(label = SP, col = model),
            position = position_dodge(width = 1), vjust = -0.2, size = 3) +
  coord_cartesian(ylim = c(0, max(sp_final$SP) * 1.15)) +
  scale_fill_manual(values = cols) +
  scale_color_manual(values = cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)) +
  labs(y = "Weight", x = "Model")
ggsave(file.path(out_dir, paste0("logic_tree_weights_", tag, ".png")),
       plot = p_lt, width = 5, height = 2.8, dpi = 500)

message("Figures written to ", out_dir)
