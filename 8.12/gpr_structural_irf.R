# Structural GPR IRF for the TVP-GVAR.
# Required ordering in the US endogenous block:
#   US_gpr, US_y, US_dp, US_r, US_de, US_deq
#
# Supports two non-US GPR specifications:
#   current_only    : GPR_t enters contemporaneously; direct GPR_{t-1} excluded.
#   current_and_lag : original unrestricted specification.
#
# Domestic dynamics may contain Ylag1 only (p=1) or Ylag1 and Ylag2 (p=2).
# Foreign dynamics remain q=1. In current_only mode the estimated Wexlag1 block
# has five columns for non-US countries, so a zero GPR-lag coefficient is inserted.
# When p=2, the entire foreign lag-2 block is zero by construction.
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

  lag_tags <- grep("^Ylag[0-9]+$", unique(dims), value = TRUE)
  if (!length(lag_tags)) stop("No domestic lag coefficient blocks were found.")
  lag_numbers <- suppressWarnings(as.integer(sub("^Ylag", "", lag_tags)))
  if (anyNA(lag_numbers)) stop("Could not parse domestic lag labels.")
  p <- max(lag_numbers)
  if (!p %in% c(1L, 2L) || !identical(sort(unique(lag_numbers)), seq_len(p))) {
    stop("Domestic lag labels must be contiguous and p must be 1 or 2.")
  }

  theta <- lapply(seq_len(p), function(kk) get3(paste0("Ylag", kk)))
  lambda <- lapply(seq_len(p), function(kk) {
    tag <- paste0("Wexlag", kk)
    if (any(dims == tag)) get3(tag) else NULL
  })

  list(
    Lambda0post = get3("Wex"),
    Lambdapost = lambda,
    Thetapost = theta,
    lag_order = p
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

# Numerical conditioning of the contemporaneous global matrix G.
# A very large condition number or a very small minimum singular value indicates
# that solving G^{-1}H can amplify otherwise moderate coefficients.
tvpgvar_G_condition <- function(G) {
  G <- as.matrix(G)
  sv <- svd(G, nu = 0, nv = 0)$d
  if (!length(sv) || any(!is.finite(sv))) {
    stop("Non-finite singular values in global contemporaneous G matrix.")
  }
  smax <- max(sv)
  smin <- min(sv)
  cond <- if (smin <= .Machine$double.eps * max(1, smax)) Inf else smax / smin
  list(
    condition_number = as.numeric(cond),
    min_singular_value = as.numeric(smin),
    max_singular_value = as.numeric(smax)
  )
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
    own_lag_i <- vector("list", length(VARi$Thetapost))

    for (kk in seq_along(VARi$Thetapost)) {
      theta <- t(VARi$Thetapost[[kk]][t, , ])
      own_lag_i[[kk]] <- theta

      if (is.null(VARi$Lambdapost[[kk]])) {
        # q=1: the lag-2 foreign block is structurally zero when p=2.
        lag_full <- matrix(
          0, nrow = k_i, ncol = nrow(Wi) - k_i
        )
      } else {
        lag_est <- t(VARi$Lambdapost[[kk]][t, , ])
        lag_full <- align_lag_block(lag_est, Wi, k_i, countries[[i]])
      }
      B <- cbind(theta, lag_full)
      if (ncol(B) != nrow(Wi)) {
        stop("Lagged B/W mismatch for ", countries[[i]],
             ": B cols=", ncol(B), ", W rows=", nrow(Wi))
      }
      H_i[[kk]] <- B %*% Wi
    }

    list(
      G = A %*% Wi,
      H = H_i,
      S = Sig_draw_i[[i]][t, , ],
      local_rho = tvpgvar_spectral_radius(own_lag_i)
    )
  }

  first <- build_local(1L)
  G <- first$G
  H <- first$H
  S_post[[1L]] <- first$S
  local_rho <- setNames(rep(NA_real_, length(countries)), countries)
  local_rho[[1L]] <- first$local_rho

  if (length(Sig_draw_i) >= 2L) {
    for (i in 2:length(Sig_draw_i)) {
      li <- build_local(i)
      G <- rbind(G, li$G)
      for (kk in seq_along(H)) H[[kk]] <- rbind(H[[kk]], li$H[[kk]])
      S_post[[i]] <- li$S
      local_rho[[i]] <- li$local_rho
    }
  }

  F <- lapply(H, function(h) solve(G, h))
  rho <- tvpgvar_spectral_radius(F)
  gcond <- tvpgvar_G_condition(G)
  max_abs_F <- max(vapply(F, function(z) max(abs(z)), numeric(1)))

  list(
    IRF_post = gpr_struct_irf(
      maxlag = length(F), G = G, F = F, sig = S_post, x = x,
      countries = countries, horizon = horz,
      shock_var = "US_gpr", normalization = normalization, shock_pct = shock_pct
    ),
    max_eigen_modulus = rho,
    G_condition_number = gcond$condition_number,
    G_min_singular_value = gcond$min_singular_value,
    G_max_singular_value = gcond$max_singular_value,
    max_abs_transition = max_abs_F,
    local_rho = local_rho
  )
}
