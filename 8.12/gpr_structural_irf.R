# Structural GPR IRF for the TVP-GVAR.
# Required ordering in the US endogenous block:
#   US_gpr, US_y, US_dp, US_r, US_de, US_deq
#
# Supports two non-US GPR specifications:
#   current_only    : GPR_t enters contemporaneously; direct GPR_{t-1} excluded.
#   current_and_lag : original unrestricted specification.
#
# In current_only mode the estimated Wexlag1 block has five columns for non-US
# countries. When constructing the global transition matrix H, a zero coefficient
# is inserted at the GPR-lag position so the local B_i matrix remains conformable
# with W_i (which still contains current global GPR as its sixth foreign variable).
#
# This version additionally reports the spectral radius of the global transition
# system for every date/posterior draw. The hard stability condition is rho < 1.

split_tvp_alpha <- function(ALPHA.draw) {
  dims <- dimnames(ALPHA.draw)[[2]]
  if (is.null(dims)) stop("ALPHA draw has no regressor dimnames.")

  get3 <- function(tag) {
    ii <- which(dims == tag)
    if (!length(ii)) stop("Missing coefficient block: ", tag)
    ALPHA.draw[, ii, , drop = FALSE]
  }

  list(
    Lambda0post = get3("Wex"),
    Lambdapost = list(get3("Wexlag1")),
    Thetapost = list(get3("Ylag1"))
  )
}

stabilize_cov <- function(S, eps = 1e-10) {
  S <- as.matrix(S)
  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  mn <- min(ev)
  if (!is.finite(mn)) stop("Non-finite covariance eigenvalue.")
  if (mn <= eps) S <- S + diag(abs(mn) + eps, nrow(S))
  S
}

# Expand the estimated lagged-foreign coefficient block so it matches W_i.
# For non-US current_only estimation:
#   estimated lag block = 5 foreign macro lags
#   W_i foreign block   = 5 foreign macro + current global GPR
# Therefore append a zero GPR-lag coefficient.
align_lag_block <- function(lag_block, W_i, k_i, country_name) {
  lag_block <- as.matrix(lag_block)
  expected_foreign <- nrow(W_i) - k_i

  if (ncol(lag_block) == expected_foreign) return(lag_block)

  if (!identical(as.character(country_name), "US") &&
      ncol(lag_block) == expected_foreign - 1L &&
      expected_foreign == 6L) {
    return(cbind(lag_block, rep(0, nrow(lag_block))))
  }

  stop(
    "Lag-block/W mismatch for ", country_name,
    ": lag columns=", ncol(lag_block),
    ", expected foreign columns=", expected_foreign,
    ", W rows=", nrow(W_i), ", own variables=", k_i
  )
}

# Spectral radius of the VAR companion matrix.
# For the current p=1 model this reduces to max(abs(eigen(F[[1]]))).
tvpgvar_spectral_radius <- function(F) {
  if (!length(F)) stop("No transition matrices supplied for stability check.")
  F <- lapply(F, as.matrix)
  K <- nrow(F[[1]])
  if (any(vapply(F, nrow, integer(1)) != K) ||
      any(vapply(F, ncol, integer(1)) != K)) {
    stop("All transition matrices must be square and have the same dimension.")
  }

  p <- length(F)
  if (p == 1L) {
    vals <- eigen(F[[1]], only.values = TRUE)$values
    rho <- max(Mod(vals))
    if (!is.finite(rho)) stop("Non-finite spectral radius.")
    return(as.numeric(rho))
  }

  top <- do.call(cbind, F)
  lower <- cbind(
    diag(K * (p - 1L)),
    matrix(0, nrow = K * (p - 1L), ncol = K)
  )
  companion <- rbind(top, lower)
  vals <- eigen(companion, only.values = TRUE)$values
  rho <- max(Mod(vals))
  if (!is.finite(rho)) stop("Non-finite companion spectral radius.")
  as.numeric(rho)
}

gpr_struct_irf <- function(maxlag, G, F, sig, x, countries,
                           horizon = 12,
                           shock_var = "US_gpr",
                           normalization = c("pct", "sd"),
                           shock_pct = 10) {
  normalization <- match.arg(normalization)
  K <- nrow(x)
  if (nrow(G) != K || ncol(G) != K) stop("Global contemporaneous G matrix is not square.")

  us_idx <- match("US", countries)
  if (is.na(us_idx)) stop("US block not found.")
  if (length(sig) != length(countries)) stop("Country covariance list and country list have different lengths.")

  S_us <- stabilize_cov(sig[[us_idx]])
  if (nrow(S_us) != 6L) stop("US covariance block must be 6x6 under GPR-first ordering.")

  # Sigma_US = L L'. With GPR first, column 1 corresponds to the recursively
  # identified GPR innovation in the US/dominant block.
  L_us <- t(chol(S_us))
  shock_us <- L_us[, 1]

  blocks <- lapply(seq_along(sig), function(i) {
    if (i == us_idx) shock_us else rep(0, nrow(sig[[i]]))
  })
  u_stack <- unlist(blocks, use.names = FALSE)
  if (length(u_stack) != K) stop("Stacked innovation dimension does not match global system.")

  impact_raw <- as.numeric(solve(G, u_stack))
  j <- match(shock_var, rownames(x))
  if (is.na(j)) stop("Shock variable not found in global vector: ", shock_var)
  if (!is.finite(impact_raw[j]) || abs(impact_raw[j]) < 1e-12) stop("Near-zero GPR impact; cannot normalize.")

  if (normalization == "pct") {
    target <- log1p(shock_pct / 100)
    impact <- impact_raw * (target / impact_raw[j])
  } else {
    impact <- if (impact_raw[j] < 0) -impact_raw else impact_raw
  }

  phi <- array(0, c(K, K, horizon + 1L))
  phi[, , 1] <- diag(K)
  if (horizon >= 1L) {
    for (h in seq_len(horizon)) {
      acc <- matrix(0, K, K)
      for (lag in seq_len(min(maxlag, h))) {
        acc <- acc + F[[lag]] %*% phi[, , h - lag + 1L]
      }
      phi[, , h + 1L] <- acc
    }
  }

  out <- sapply(seq_len(horizon + 1L), function(h) phi[, , h] %*% impact)
  rownames(out) <- rownames(x)
  colnames(out) <- 0:horizon
  out
}

get_gpr_struct_irfa_t <- function(tt, draw_i, Sig_draw_i, x, globalG, countries,
                                  horz = 12, normalization = "pct", shock_pct = 10) {
  t <- tt
  S_post <- list()

  build_local <- function(i) {
    VARi <- split_tvp_alpha(draw_i[[i]])
    k_i <- dim(VARi$Lambda0post)[[3]]
    Wi <- globalG[[i]]

    A <- cbind(diag(k_i), -t(VARi$Lambda0post[t, , ]))
    if (ncol(A) != nrow(Wi)) {
      stop("Contemporaneous A/W mismatch for ", countries[[i]],
           ": A cols=", ncol(A), ", W rows=", nrow(Wi))
    }

    H_i <- vector("list", length(VARi$Thetapost))
    for (kk in seq_along(VARi$Thetapost)) {
      theta <- t(VARi$Thetapost[[kk]][t, , ])
      lag_est <- t(VARi$Lambdapost[[kk]][t, , ])
      lag_full <- align_lag_block(lag_est, Wi, k_i, countries[[i]])
      B <- cbind(theta, lag_full)
      if (ncol(B) != nrow(Wi)) {
        stop("Lagged B/W mismatch for ", countries[[i]],
             ": B cols=", ncol(B), ", W rows=", nrow(Wi))
      }
      H_i[[kk]] <- B %*% Wi
    }

    list(G = A %*% Wi, H = H_i, S = Sig_draw_i[[i]][t, , ])
  }

  first <- build_local(1L)
  G <- first$G
  H <- first$H
  S_post[[1L]] <- first$S

  if (length(Sig_draw_i) >= 2L) {
    for (i in 2:length(Sig_draw_i)) {
      li <- build_local(i)
      G <- rbind(G, li$G)
      for (kk in seq_along(H)) H[[kk]] <- rbind(H[[kk]], li$H[[kk]])
      S_post[[i]] <- li$S
    }
  }

  F <- lapply(H, function(h) solve(G, h))
  rho <- tvpgvar_spectral_radius(F)

  list(
    IRF_post = gpr_struct_irf(
      maxlag = length(F), G = G, F = F, sig = S_post, x = x,
      countries = countries, horizon = horz,
      shock_var = "US_gpr", normalization = normalization, shock_pct = shock_pct
    ),
    max_eigen_modulus = rho
  )
}
