normalize_vector <- function(x) {
  vec_norm <- sqrt(sum(x^2, na.rm = TRUE))
  if (is.na(vec_norm) || vec_norm == 0) return(rep(0, length(x)))
  x / vec_norm
}

SP_vecnorm_detail <- function(LLHs, EDRs) {
  LLHs_norm <- apply(LLHs, 2, normalize_vector)
  EDRs_norm <- apply(EDRs, 2, normalize_vector)
  spj <- matrix(rowSums(LLHs_norm * EDRs_norm, na.rm = TRUE), ncol = 1,
                dimnames = list(rownames(LLHs), NULL))
  upper <- max(spj) / spj
  SP <- upper / sum(upper)
  list(LLHs_norm = LLHs_norm, EDRs_norm = EDRs_norm, spj = spj, SP = SP)
}
