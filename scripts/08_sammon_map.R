# Sammon's maps of the median predictions of the crustal GMPEs in the hazard-controlling
# scenarios (Section 6.5; Figure 11a-b).

source("config.R", local = TRUE)
suppressPackageStartupMessages({
  library(GMPEhaz)
  library(MASS)
})

if (!exists("IMT"))           IMT           <- "PGA"
if (!exists("SCEN_YEAR_MAX")) SCEN_YEAR_MAX <- 2020
stopifnot(IMT %in% c("PGA", "SA1.0"))

RP_LABEL <- "RP = 475 yr"
RP_TAG   <- "RP475"
T_IM     <- if (IMT == "PGA") 0 else 1.0
MW       <- if (IMT == "PGA") seq(5.75, 7.25, by = 0.25) else seq(6.25, 7.25, by = 0.25)
RRUPS    <- seq(10, 50, by = 10)
VS30S    <- c(360, 459.776, 557.918, 566.473, 633.091, 552.289, 565.531, 537.573,
              608.287, 232.077, 465.929, 525.994, 293.668, 437.539, 505.107, 245.398)

FAULT_TYPE <- "RV"
FTYPE <- 1
FRV   <- 1
FNM   <- 0
VS30_CLASS <- 1
SOURCETYPE <- 0
MSASFLAG   <- 0
BASINFLAG  <- 0
SITE_SIDE  <- "footwall"

out_dir <- file.path("figures", "sammon")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ff <- read.csv(DATA_FILE, stringsAsFactors = FALSE)
yr <- as.integer(substr(ff$EQ_ID, 1, 4))
rv_ev <- ff[which(yr <= SCEN_YEAR_MAX & ff$EQ_type == "shallow crustal" & ff$ftype == FTYPE & ff$Mw >= 6), ]
rv_ev <- rv_ev[!duplicated(rv_ev$EQ_ID), ]
DIP      <- median(rv_ev$Dip)
ZTOR_FIX <- median(rv_ev$Ztor)
message(sprintf("%d reverse events: Dip = %.2f, Ztor = %.2f km", nrow(rv_ev), DIP, ZTOR_FIX))

z1_from_vs30 <- function(vs30) exp(28.5 - (3.82 / 8) * log(vs30^8 + 378.7^8))
z25_from_z1  <- function(z1) (519 + 3.595 * z1) / 1000
site_tab <- data.frame(vs30 = VS30S, z1 = z1_from_vs30(VS30S))
site_tab$z25 <- z25_from_z1(site_tab$z1)

fault_plane <- function(mw, dip, fn, fr, ztor) {
  A  <- 10^(-3.49 + 0.91 * mw)
  AR <- 10^((0.01752 - 0.00472 * fn - 0.01099 * fr) * (mw - 4)^3.097)
  W  <- sqrt(A / AR)
  sd <- sin(dip * pi / 180)
  depth <- ztor + W / 2 * sd
  ztop <- depth - W / 2 * sd
  zbot <- depth + W / 2 * sd
  zt   <- max(ztop, 0)
  list(mw = mw, W = W, dip = dip, depth = depth, ztor = zt, W_eff = (zbot - zt) / sd)
}

site_distance <- function(rx, fp) {
  d <- fp$dip * pi / 180
  s <- min(max(rx * cos(d) - fp$ztor * sin(d), 0), fp$W_eff)
  c(Rrup = sqrt((rx - s * cos(d))^2 + (fp$ztor + s * sin(d))^2),
    Rjb  = max(0, -rx, rx - fp$W_eff * cos(d)),
    Rx   = rx,
    Ry   = 0)
}

geometry_for_rrup <- function(rrup, fp, side) {
  if (rrup < fp$ztor) return(NULL)
  if (side == "footwall") {
    rx <- -sqrt(rrup^2 - fp$ztor^2)
  } else {
    f  <- function(x) site_distance(x, fp)["Rrup"] - rrup
    rx <- if (abs(f(0)) < 1e-9) 0 else uniroot(f, c(0, rrup + fp$W_eff + 1), tol = 1e-10)$root
  }
  site_distance(rx, fp)
}

fps <- lapply(MW, function(m) fault_plane(m, DIP, FNM, FRV, ztor = ZTOR_FIX))

geo <- do.call(rbind, lapply(seq_along(MW), function(k)
  do.call(rbind, lapply(RRUPS, function(r) {
    g <- geometry_for_rrup(r, fps[[k]], SITE_SIDE)
    if (is.null(g)) return(NULL)
    data.frame(Mw = MW[k], t(g), Ztor = fps[[k]]$ztor, hypo = fps[[k]]$depth,
               width = 10^(-1.669 + 0.435 * MW[k]))
  }))))
geo$hwflag <- as.integer(geo$Rx > 0)
if (nrow(geo) < length(MW) * length(RRUPS))
  warning(length(MW) * length(RRUPS) - nrow(geo), " magnitude-distance combinations removed (Rrup < Ztor).")

scen <- merge(geo, site_tab, by = NULL)
scen <- scen[order(scen$Mw, scen$Rrup, scen$vs30), ]
rownames(scen) <- NULL
write.csv(scen, file.path(out_dir, sprintf("scenario_grid_%s.csv", IMT)), row.names = FALSE)

TT <- T_IM
gmpe_fns <- list(
  ASK14 = function(s) ASK14(s$Mw, s$Rrup, s$Rjb, TT, DIP, FTYPE, s$width, s$vs30, s$hwflag, s$Ztor,
                            VS30_CLASS, s$z1, s$Rx, s$Ry, regionflag = 1, msasflag = 0),
  BSSA14 = function(s) BSSA14(s$Mw, s$Rjb, TT, s$vs30, FTYPE, s$z1, regionflag = 0, BASINFLAG),
  CB14 = function(s) CB14(s$Mw, s$Rrup, s$Rjb, FTYPE, TT, s$vs30, s$Ztor, s$z25, DIP, s$hypo,
                          hwflag = 0, s$Rx, s$width, regionflag = 0),
  CY14 = function(s) CY14(s$Mw, s$Rrup, s$Rjb, TT, s$vs30, DIP, s$Ztor, FTYPE, s$z1, VS30_CLASS,
                          hwflag = 0, s$Rx, regionflag = 0),
  I14 = function(s) I14(s$Mw, s$Rjb, TT, FTYPE, s$vs30),
  Chao20 = function(s) Chao2018(s$Mw, s$Rrup, TT, FTYPE, s$vs30, VS30_CLASS, s$Ztor, s$z1,
                                SOURCETYPE, MSASFLAG),
  PhungCru20 = function(s) PhungCru18(s$Mw, s$Rrup, s$Rjb, s$Rx, TT, s$vs30, DIP, s$Ztor, FTYPE, s$z1,
                                      VS30_CLASS, regionflag = 1, s$hwflag),
  Lin09 = function(s) Lin09(s$Mw, s$Rrup, TT, s$vs30, FTYPE),
  Lin11 = function(s) Lin11(s$Mw, s$Rrup, TT, s$hwflag, s$vs30)
)

scen_list <- lapply(seq_len(nrow(scen)), function(i) as.list(scen[i, ]))
pred_matrix <- do.call(rbind, lapply(names(gmpe_fns), function(nm)
  vapply(scen_list, function(s) as.numeric(unlist(gmpe_fns[[nm]](s))["lnY"]) - LN_G, numeric(1))))
rownames(pred_matrix) <- names(gmpe_fns)
colnames(pred_matrix) <- sprintf("M%.2f_R%g_V%.1f", scen$Mw, scen$Rrup, scen$vs30)
if (anyNA(pred_matrix)) stop("A GMPE returned NA for at least one scenario.")

write.csv(cbind(scen[, c("Mw", "Rrup", "Rjb", "Rx", "vs30", "z1", "z25", "Ztor")], t(exp(pred_matrix))),
          file.path(out_dir, sprintf("median_by_model_%s.csv", IMT)), row.names = FALSE)

d <- dist(pred_matrix) / sqrt(ncol(pred_matrix))
write.csv(as.data.frame(as.matrix(d)), file.path(out_dir, sprintf("gmpe_distances_%s.csv", IMT)))

set.seed(42)
sm  <- sammon(d, k = 2, niter = 500, magic = 0.2, trace = FALSE)
pts <- sm$points
rownames(pts) <- rownames(pred_matrix)

amp_dev <- rowMeans(sweep(pred_matrix, 2, colMeans(pred_matrix), "-"))
write.csv(data.frame(model = names(amp_dev), ratio_to_mean = exp(amp_dev),
                     dim1 = pts[, 1], dim2 = pts[, 2]),
          file.path(out_dir, sprintf("sammon_coordinates_%s.csv", IMT)), row.names = FALSE)
message(sprintf("Stress = %.5f; cor(Dim 1, mean ln deviation) = %.3f", sm$stress, cor(pts[, 1], amp_dev)))

col_map <- c(ASK14 = "#E41A1C", BSSA14 = "#A65628", CB14 = "#F781BF", CY14 = "#1B9E77", I14 = "#E6AB02",
             Chao20 = "#377EB8", PhungCru20 = "#984EA3", Lin09 = "#4DAF4A", Lin11 = "#FF7F00")
TAIWAN_MODELS <- c("Chao20", "PhungCru20", "Lin09", "Lin11")
GLOBAL_MODELS <- c("ASK14", "BSSA14", "CB14", "CY14", "I14")
pch_map <- c(setNames(rep(21, 4), TAIWAN_MODELS), setNames(rep(23, 5), GLOBAL_MODELS))
cex_map <- c(setNames(rep(2.6, 4), TAIWAN_MODELS), setNames(rep(2.0, 5), GLOBAL_MODELS))

amax <- max(abs(pts)) * 1.2
step <- if (amax <= 0.5) 0.1 else if (amax <= 1.5) 0.25 else 0.5
lim  <- max(step * 3, ceiling(amax / step) * step)
BREAKS <- seq(-lim, lim, by = step)

fmt_range <- function(x, digits = 0)
  if (length(unique(x)) == 1) sprintf(paste0("%.", digits, "f"), x[1]) else
    sprintf(paste0("%.", digits, "f-%.", digits, "f"), min(x), max(x))
im_label  <- if (IMT == "PGA") "PGA" else "SA(1.0 s)"
title_txt <- paste0("Sammon's map -- ", RP_LABEL, " scenarios (reverse, ", im_label, ")\n",
                    sprintf("Mw %s (n=%d) x Rrup %s km (n=%d) x Vs30 %s m/s (n=%d)",
                            fmt_range(MW, 2), length(MW), fmt_range(RRUPS), length(RRUPS),
                            fmt_range(VS30S), length(VS30S)))

outfile <- file.path(out_dir, sprintf("sammon_%s_%s_%s.jpg", RP_TAG, FAULT_TYPE, IMT))
jpeg(outfile, width = 1100, height = 800, res = 120, quality = 100)
layout(matrix(1:2, nrow = 1), widths = c(4, 1.5))

par(mar = c(5, 5, 4, 1), pty = "s")
plot(pts, type = "n", xlim = c(-lim, lim), ylim = c(-lim, lim), asp = 1, xaxs = "i", yaxs = "i",
     axes = FALSE, xlab = "Sammon dimension 1 (ln units)", ylab = "Sammon dimension 2 (ln units)",
     main = title_txt, cex.main = 0.8)
axis(1, at = BREAKS, labels = format(BREAKS), cex.axis = 0.85)
axis(2, at = BREAKS, labels = format(BREAKS), cex.axis = 0.85, las = 1)
abline(v = BREAKS, h = BREAKS, col = "grey88")
box()
draw_order <- c(GLOBAL_MODELS, TAIWAN_MODELS)
points(pts[draw_order, 1], pts[draw_order, 2], pch = pch_map[draw_order], bg = col_map[draw_order],
       col = "white", cex = cex_map[draw_order], lwd = 1.5)
legend("bottomright", legend = sprintf("Stress = %.4f", sm$stress), bty = "n", cex = 0.85)

par(mar = c(0, 0, 0, 0), pty = "m")
plot.new()
plot.window(xlim = c(0, 1), ylim = c(0, 1))
lg <- legend(x = 0, y = 0.85, legend = TAIWAN_MODELS, pch = 21, pt.bg = col_map[TAIWAN_MODELS],
             col = col_map[TAIWAN_MODELS], pt.cex = 2.0, cex = 0.95, y.intersp = 1.6, x.intersp = 1.3,
             title = "Taiwan models", title.adj = 0, title.font = 2, bty = "n")
legend(x = 0, y = lg$rect$top - lg$rect$h - 0.05, legend = GLOBAL_MODELS, pch = 23,
       pt.bg = col_map[GLOBAL_MODELS], col = col_map[GLOBAL_MODELS], pt.cex = 2.2, cex = 0.95,
       y.intersp = 1.6, x.intersp = 1.3, title = "Global models", title.adj = 0, title.font = 2, bty = "n")
dev.off()
message("Saved ", outfile)
