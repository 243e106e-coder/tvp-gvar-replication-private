gpr_girf <- function(maxlag, G, F, sig, x, horizon = 12,
                     shock_var = "US_gpr") {
  require(Matrix)
  K <- nrow(x)
  block_cov <- as.matrix(bdiag(sig))
  invG <- solve(G)
  omega <- invG %*% block_cov %*% t(invG)
  j <- match(shock_var, rownames(x))
  if (is.na(j)) stop("Shock variable not found: ", shock_var)
  impact <- omega[, j] / sqrt(omega[j, j])

  phi <- array(0, c(K, K, horizon + 1))
  phi[, , 1] <- diag(K)
  if (horizon >= 1) {
    for (h in 1:horizon) {
      acc <- matrix(0, K, K)
      for (lag in seq_len(min(maxlag, h))) {
        acc <- acc + F[[lag]] %*% phi[, , h - lag + 1]
      }
      phi[, , h + 1] <- acc
    }
  }
  out <- sapply(seq_len(horizon + 1), function(h) phi[, , h] %*% impact)
  rownames(out) <- rownames(x)
  colnames(out) <- 0:horizon
  out
}

