LLH_weighted <- function(res, sigma = rep(1, length(res)), w = rep(1, length(res)),
                         include_log_sigma = FALSE) {
  if (length(sigma) == 1) sigma <- rep(sigma, length(res))
  stopifnot(length(res) == length(w), length(res) == length(sigma))
  ok <- is.finite(res) & is.finite(sigma) & is.finite(w) & sigma > 0
  res <- res[ok]; sigma <- sigma[ok]; w <- w[ok]
  if (length(res) == 0) return(NA_real_)
  vals <- if (include_log_sigma) -log2(dnorm(res, mean = 0, sd = sigma)) else -log2(dnorm(res / sigma))
  sum(w * vals) / sum(w)
}
