# Calibration plots and decomposed residuals for the 16 GMPEs (Figures 8, X5-X16).

source("config.R", local = TRUE)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggtext)
  library(MethComp)
  library(patchwork)
  library(scales)
})

SHARE_ACROSS_PERIODS <- FALSE
PER_TAG <- c("PGA", "T0.17s", "T1.0s")
MODEL_ORDER <- c("Chao20(CRU)", "Lin09(CRU)", "Lin11(CRU)", "PhungCru20(CRU)",
                 "ASK14(CRU)", "BSSA14(CRU)", "CB14(CRU)", "CY14(CRU)",
                 "I14(CRU)", "Chao20(SUB)", "LL08(SUB)", "PhungSub20(SUB)",
                 "AG2020(SUB)", "KBCG2020(SUB)", "PSHAB2022(SUB)", "SI2022(SUB)")

dir_rds  <- file.path(INTER_DIR, "calib_resid")
dir_comb <- file.path(FIG_DIR, "residuals_combined")
tag_of   <- function(m) if (grepl("\\(SUB\\)", m)) "sub" else "cru"
calib_dir_of  <- function(m) file.path(FIG_DIR, paste0("calibration_", tag_of(m)))
decomp_dir_of <- function(m) file.path(FIG_DIR, paste0("residual_decomposition_", tag_of(m)))
for (d in unique(c(dir_comb, sapply(MODEL_ORDER, calib_dir_of), sapply(MODEL_ORDER, decomp_dir_of))))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

title_theme <- theme(plot.title = element_text(face = "bold", hjust = 0.5))
y_resid <- scale_y_continuous(limits = c(-2.2, 2), breaks = seq(-2, 2, 1))

stat_box <- function(xmin, xmax, xc, n_lab, n, mean_lab, m, sd_lab, s, coefs) list(
  annotate("rect", xmin = xmin, xmax = xmax, ymin = -1.1, ymax = -1.9, fill = "white",
           color = "black", linewidth = 0.3, alpha = 0.85),
  annotate("text", x = xc, y = -1.2, label = paste(n_lab, ":", n), size = 3.3),
  annotate("text", x = xc, y = -1.4, label = paste(mean_lab, ":", round(m, 3)), size = 3.3, col = "red"),
  annotate("text", x = xc, y = -1.6, label = paste(sd_lab, ":", round(s, 3)), size = 3.3),
  annotate("text", x = xc, y = -1.8, label = paste0("Regression : y = ", round(coefs[2], 3), "x + ",
                                                     round(coefs[1], 3)), size = 3.2, col = "blue"))

base_calibration <- function(o) {
  df <- o$calib
  df <- df[is.finite(df$obs) & is.finite(df$pred) & df$obs > 0 & df$pred > 0, ]
  df$X <- log(df$obs)
  df$Y <- log(df$pred)
  lim  <- c(-10, 1)
  reg  <- data.frame(t(Deming(df$X, df$Y)))
  xc   <- lim[1] + 0.5 * diff(lim)
  y_at <- function(f) lim[1] + f * diff(lim)
  leg  <- paste0("<span style='font-size:12pt;color:grey30'>- <span style='font-size:7pt;color:black'>pre/obs=1<br>",
                 "<span style='font-size:12pt;color:blue'>- <span style='font-size:7pt;color:black'>regression")
  ggplot() +
    geom_bin2d(data = df, aes(x = X, y = Y), bins = 100) +
    scale_x_continuous(limits = lim, breaks = seq(lim[1], lim[2], 1)) +
    scale_y_continuous(limits = lim, breaks = seq(lim[1], lim[2], 1)) +
    geom_abline(intercept = 0, slope = 1, col = "grey30", lwd = 0.5) +
    geom_abline(intercept = reg$Intercept, slope = reg$Slope, col = "blue", lwd = 0.5) +
    annotate("label", x = xc, y = y_at(0.086), size = 2.7, fill = "white",
             label = paste0("      Mean of Residual : ", round(mean(df$X - df$Y), 3), "     \n\n")) +
    annotate("text", x = xc, y = y_at(0.060), size = 2.7,
             label = paste0("Sigma (ln) : ", round(sd(df$X - df$Y), 3), "\n")) +
    annotate("text", x = xc, y = y_at(0.036), size = 2.7, col = "blue",
             label = paste0("Regression : y = ", round(reg$Slope, 3), "x + ", round(reg$Intercept, 3))) +
    geom_richtext(data = data.frame(x = -Inf, y = Inf, label = leg), aes(x = x, y = y, label = label),
                  hjust = 0, vjust = 1, inherit.aes = FALSE) +
    labs(title = paste0("model_", o$model_name),
         x = paste0("ln(", o$im_label, ") Observation"), y = paste0("ln(", o$im_label, ") Prediction")) +
    theme_bw() + title_theme
}

base_between <- function(o) {
  ev  <- o$event
  reg <- coef(lm(dB ~ Mw, data = ev))
  xc  <- mean(c(min(2, floor(min(ev$Mw))), max(8, ceiling(max(ev$Mw)))))
  ggplot(ev, aes(x = Mw, y = dB)) +
    geom_bin2d(bins = 100) +
    theme_bw() + title_theme +
    labs(title = paste0("model_", o$model_name), y = "Between-event Residual", x = "Magnitude (Mw)") +
    scale_x_continuous(limits = c(2, 8), breaks = 2:8) + y_resid +
    geom_hline(yintercept = 0, lwd = 0.3) +
    geom_hline(yintercept = mean(ev$dB), col = "red", lty = 2, lwd = 0.6) +
    geom_abline(intercept = reg[1], slope = reg[2], col = "blue", lwd = 0.5) +
    stat_box(3.2, 6.8, xc, "N events", nrow(ev), "Mean dB", mean(ev$dB), "Sigma(dB)", sd(ev$dB), reg)
}

base_within <- function(o) {
  df  <- o$record[is.finite(o$record$dist) & o$record$dist > 0, ]
  reg <- coef(lm(dWS ~ log10(dist), data = df))
  R_seq <- 10^seq(log10(max(0.8, min(df$dist))), log10(500), length.out = 300)
  brk <- c(1, 2, 5, 10, 20, 50, 100, 200, 500)
  ggplot(df, aes(x = dist, y = dWS)) +
    geom_bin2d(bins = 100) +
    theme_bw() + title_theme +
    labs(title = paste0("model_", o$model_name), y = "Single-station Within-event Residual",
         x = paste0(o$dist_label, " (km)")) +
    scale_x_log10(limits = c(0.8, 500), breaks = brk, labels = brk) + y_resid +
    geom_hline(yintercept = 0, lwd = 0.3) +
    geom_hline(yintercept = mean(df$dWS), col = "red", lty = 2, lwd = 0.6) +
    geom_line(data = data.frame(x = R_seq, y = reg[1] + reg[2] * log10(R_seq)), aes(x = x, y = y),
              col = "blue", lwd = 0.5, inherit.aes = FALSE) +
    stat_box(2.8, 300, 30, "N records", nrow(df), "Mean dWS", mean(df$dWS), "Sigma(dWS)", sd(df$dWS), reg)
}

base_site <- function(o) {
  st  <- o$station
  reg <- coef(lm(dS2S ~ log10(Vs30), data = st))
  V_seq <- 10^seq(log10(100), log10(2000), length.out = 300)
  brk <- c(100, 200, 300, 500, 700, 1000, 1500, 2000)
  ggplot(st, aes(x = Vs30, y = dS2S)) +
    geom_bin2d(bins = 100) +
    theme_bw() + title_theme +
    labs(title = paste0("model_", o$model_name), y = "Site/Station Term", x = "Vs30 (m/s)") +
    scale_x_log10(limits = c(100, 2000), breaks = brk, labels = brk) + y_resid +
    geom_hline(yintercept = 0, lwd = 0.3) +
    geom_hline(yintercept = mean(st$dS2S), col = "red", lty = 2, lwd = 0.6) +
    geom_line(data = data.frame(x = V_seq, y = reg[1] + reg[2] * log10(V_seq)), aes(x = x, y = y),
              col = "blue", lwd = 0.5, inherit.aes = FALSE) +
    stat_box(170, 1400, 500, "N stations", nrow(st), "Mean dS2S", mean(st$dS2S), "Sigma(dS2S)", sd(st$dS2S), reg)
}

max_count <- function(p) suppressWarnings(max(layer_data(p, 1)$count, na.rm = TRUE))
fill_scale <- function(gmax, name)
  scale_fill_continuous(type = "viridis", limits = c(1, max(gmax, 2)), oob = squish, name = name)

PLOT_TYPES <- list(
  calibration = list(base = base_calibration, legend = "Count", w = 4.2, h = 3.6, pdf = TRUE,
                     out = function(o) file.path(calib_dir_of(o$model_name),
                                                 paste0("calibration_", o$model_name, "_", o$per_tag))),
  between = list(base = base_between, legend = "Events", w = 4.8, h = 4.5, pdf = FALSE,
                 out = function(o) file.path(decomp_dir_of(o$model_name),
                                             paste0(o$model_name, "_", o$per_tag, "_between_event_vs_Mw"))),
  within = list(base = base_within, legend = "Records", w = 4.8, h = 4.5, pdf = FALSE,
                out = function(o) file.path(decomp_dir_of(o$model_name),
                                            paste0(o$model_name, "_", o$per_tag, "_within_event_vs_", o$dist_label))),
  site = list(base = base_site, legend = "Stations", w = 4.8, h = 4.5, pdf = FALSE,
              out = function(o) file.path(decomp_dir_of(o$model_name),
                                          paste0(o$model_name, "_", o$per_tag, "_site_term_vs_Vs30")))
)

read_obj <- function(m, per) {
  f <- file.path(dir_rds, paste0(m, "_", per, ".rds"))
  if (!file.exists(f)) { warning("Missing ", f); return(NULL) }
  readRDS(f)
}
objs <- lapply(setNames(PER_TAG, PER_TAG), function(per)
  lapply(setNames(MODEL_ORDER, MODEL_ORDER), function(m) read_obj(m, per)))

for (ty in names(PLOT_TYPES)) {
  cfg   <- PLOT_TYPES[[ty]]
  bases <- lapply(objs, function(pl) lapply(pl, function(o) if (is.null(o)) NULL else cfg$base(o)))
  maxes <- lapply(bases, function(pl) sapply(pl, function(p) if (is.null(p)) NA else max_count(p)))
  gmax_all <- max(unlist(maxes), na.rm = TRUE)

  for (per in PER_TAG) {
    fs <- fill_scale(if (SHARE_ACROSS_PERIODS) gmax_all else max(maxes[[per]], na.rm = TRUE), cfg$legend)
    for (m in MODEL_ORDER) {
      p <- bases[[per]][[m]]
      if (is.null(p)) next
      f <- cfg$out(objs[[per]][[m]])
      ggsave(paste0(f, ".jpg"), plot = p + fs, width = cfg$w, height = cfg$h, dpi = 300)
      if (cfg$pdf) ggsave(paste0(f, ".pdf"), plot = p + fs, width = cfg$w, height = cfg$h)
    }
    panels <- lapply(MODEL_ORDER, function(m) {
      p <- bases[[per]][[m]]
      if (is.null(p)) plot_spacer() else p + fs + theme(legend.position = "none")
    })
    f_comb <- file.path(dir_comb, paste0(ty, "_16models_", per))
    comb <- wrap_plots(panels, ncol = 4)
    ggsave(paste0(f_comb, ".jpg"), plot = comb, width = cfg$w * 4 * 0.85, height = cfg$h * 4, dpi = 300)
    if (cfg$pdf) ggsave(paste0(f_comb, ".pdf"), plot = comb, width = cfg$w * 4 * 0.85, height = cfg$h * 4)
  }
  message(ty, " figures done")
}
