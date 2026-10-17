EDR_weighted <- function(obs, est, sd, w = rep(1, length(obs))) {
  stopifnot(length(obs) == length(w))
  dd <- 0.01
  x  <- 3
  N  <- length(obs)
  sw <- sum(w)

  mua <- sum(w * obs) / sw
  muY <- sum(w * est) / sw
  B1  <- sum(w * (obs - mua) * (est - muY)) / sum(w * (obs - mua)^2)
  B0  <- muY - B1 * mua
  Yc  <- est - (B0 + B1 * obs - obs)
  Kappa <- sum(w * (obs - est)^2) / sum(w * (obs - Yc)^2)

  mu_D <- obs - est
  s_D  <- sd
  dcmax <- ceiling(max(abs(obs - (est - x * sd)), abs(obs - (est + x * sd))))
  nd <- length(seq(dd / 2, dcmax, dd))

  MDE <- 0
  for (j in seq_len(nd)) {
    di <- dd / 2 + (j - 1) * dd
    d1 <- di - dd / 2
    d2 <- di + dd / 2
    P1 <- pnorm((d1 - mu_D) / s_D) - pnorm((-d1 - mu_D) / s_D)
    P2 <- pnorm((d2 - mu_D) / s_D) - pnorm((-d2 - mu_D) / s_D)
    MDE <- MDE + (P2 - P1) * di
  }

  list(MDE_norm = sqrt(sum(w * MDE^2) / sw),
       Kappa_sq = sqrt(Kappa),
       EDR      = sqrt(Kappa * sum(w * MDE^2) / sw))
}
