run <- function(script, ...) {
  env <- new.env(parent = globalenv())
  list2env(list(...), envir = env)
  message("\n=== ", script, " ", paste(names(list(...)), unlist(list(...)), sep = "=", collapse = " "))
  source(script, local = env)
  invisible(NULL)
}

periods <- list(c(1992, 2024), c(1992, 2016), c(2017, 2024))

if (!file.exists("data/TSMIP_flatfile_1992_2024.csv")) run("scripts/00_prepare_flatfile.R")

for (tec in c("cru", "sub")) {
  for (p in periods) {
    run("scripts/01_gmpe_predictions.R", TECTONIC = tec, YEAR_MIN = p[1], YEAR_MAX = p[2])
    run("scripts/02_score_weights.R",    TECTONIC = tec, YEAR_MIN = p[1], YEAR_MAX = p[2])
    run("scripts/03_plot_scores.R",      TECTONIC = tec, YEAR_MIN = p[1], YEAR_MAX = p[2])
  }
  run("scripts/04_residual_decomposition.R", TECTONIC = tec)
  run("scripts/06_bootstrap.R",              TECTONIC = tec)
  run("scripts/07_bootstrap_figures.R",      TECTONIC = tec)
}
run("scripts/05_plot_residuals.R")

for (im in c("PGA", "SA1.0")) run("scripts/08_sammon_map.R", IMT = im)
