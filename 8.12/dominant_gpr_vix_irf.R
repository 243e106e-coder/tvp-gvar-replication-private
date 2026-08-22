# =============================================================================
# Structural GPR impulse response for a separate dominant unit [GPR,VIX].
#
# Local structural identification ONLY in dominant unit:
#   GL_gpr -> GL_vix
#
# Other country blocks are not Cholesky-ordered for the GPR shock.  The shock is
# injected into the dominant-unit innovation and propagated through the stacked
# contemporaneous GVAR system and dynamic transition matrices.
# =============================================================================

stabilize_cov <- function(S, eps = 1e-10) {
  S <- as.matrix(S); S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  mn <- min(ev)
  if (!is.finite(mn)) stop("Non-finite covariance eigenvalue.")
  if (mn <= eps) S <- S + diag(abs(mn) + eps, nrow(S))
  S
}

split_tvp_alpha_dominant <- function(A) {
  dims <- dimnames(A)[[2]]
  eq <- dimnames(A)[[3]]
  if (is.null(dims) || is.null(eq)) stop("ALPHA array lacks dimnames.")
  TT <- dim(A)[1]; k <- dim(A)[3]

  block <- function(tag, optional = FALSE) {
    ii <- which(dims == tag)
    if (!length(ii)) {
      if (!optional) stop("Missing coefficient block: ", tag)
      return(array(0, dim = c(TT, 0L, k)))
    }
    A[, ii, , drop = FALSE]
  }

  lag_tags <- grep("^Ylag[0-9]+$", unique(dims), value = TRUE)
  if (!length(lag_tags)) stop("No domestic lag blocks.")
  lag_nums <- sort(unique(as.integer(sub("^Ylag", "", lag_tags))))
  if (anyNA(lag_nums) || !identical(lag_nums, seq_len(max(lag_nums)))) {
    stop("Domestic lag labels are not contiguous.")
  }
  p <- max(lag_nums)

  list(
    Lambda0 = block("Wex", optional = TRUE),
    Lambda = lapply(seq_len(p), function(j) block(paste0("Wexlag", j), optional = TRUE)),
    Theta = lapply(seq_len(p), function(j) block(paste0("Ylag", j))),
    p = p,
    eq = eq
  )
}

slice_coef <- function(arr, tt, nrow_out) {
  if (dim(arr)[2] == 0L) return(matrix(0, nrow = nrow_out, ncol = 0L))
  z <- arr[tt, , , drop = FALSE]
  z <- matrix(z, nrow = dim(arr)[2], ncol = nrow_out)
  t(z)
}

tvpgvar_spectral_radius <- function(F) {
  if (!length(F)) stop("No transition matrices.")
  K <- nrow(F[[1]]); p <- length(F)
  if (p == 1L) return(as.numeric(max(Mod(eigen(F[[1]], only.values = TRUE)$values))))
  top <- do.call(cbind, F)
  lower <- cbind(diag(K * (p - 1L)), matrix(0, K * (p - 1L), K))
  companion <- rbind(top, lower)
  as.numeric(max(Mod(eigen(companion, only.values = TRUE)$values)))
}

tvpgvar_G_condition <- function(G) {
  sv <- svd(as.matrix(G), nu = 0, nv = 0)$d
  if (!length(sv) || any(!is.finite(sv))) stop("Non-finite G singular values.")
  smax <- max(sv); smin <- min(sv)
  cond <- if (smin <= .Machine$double.eps * max(1, smax)) Inf else smax / smin
  list(condition_number = cond, min_singular_value = smin, max_singular_value = smax)
}

dominant_gpr_struct_irf <- function(G, F, sig, x, units,
                                    horizon = 12,
                                    shock_var = "GL_gpr",
                                    shock_pct = 10) {
  K <- nrow(x)
  if (!all(dim(G) == c(K, K))) stop("G must be KxK.")

  gl_i <- match("GL", units)
  if (is.na(gl_i)) stop("Dominant unit GL not found.")
  S <- stabilize_cov(sig[[gl_i]])
  if (!all(dim(S) == c(2L, 2L))) stop("Dominant covariance must be 2x2.")

  # Recursive identification inside dominant unit only: GPR -> VIX.
  L <- t(chol(S))
  shock_gl <- L[, 1]

  blocks <- lapply(seq_along(sig), function(i) {
    if (i == gl_i) shock_gl else rep(0, nrow(sig[[i]]))
  })
  u <- unlist(blocks, use.names = FALSE)
  if (length(u) != K) stop("Innovation stack dimension mismatch.")

  impact_raw <- as.numeric(solve(G, u))
  j <- match(shock_var, rownames(x))
  if (is.na(j) || !is.finite(impact_raw[j]) || abs(impact_raw[j]) < 1e-12) {
    stop("Cannot normalize GL_gpr impact.")
  }
  target <- log1p(shock_pct / 100)
  impact <- impact_raw * (target / impact_raw[j])

  phi <- array(0, c(K, K, horizon + 1L))
  phi[, , 1L] <- diag(K)
  if (horizon >= 1L) {
    for (h in seq_len(horizon)) {
      acc <- matrix(0, K, K)
      for (lag in seq_len(min(length(F), h))) {
        acc <- acc + F[[lag]] %*% phi[, , h - lag + 1L]
      }
      phi[, , h + 1L] <- acc
    }
  }

  out <- sapply(seq_len(horizon + 1L), function(h) phi[, , h] %*% impact)
  rownames(out) <- rownames(x); colnames(out) <- 0:horizon
  out
}

get_dominant_gpr_irf_t <- function(tt, draw_i, Sig_draw_i, x, globalG, units,
                                   horizon = 12, shock_pct = 10) {
  locals <- lapply(seq_along(units), function(i) {
    V <- split_tvp_alpha_dominant(draw_i[[i]])
    k_i <- length(V$eq)
    Wi <- as.matrix(globalG[[i]])

    lam0 <- slice_coef(V$Lambda0, tt, k_i)
    A <- cbind(diag(k_i), -lam0)
    if (ncol(A) != nrow(Wi)) {
      stop("A/W mismatch for ", units[[i]], ": ", ncol(A), " vs ", nrow(Wi))
    }

    H <- vector("list", V$p)
    own_lags <- vector("list", V$p)
    for (kk in seq_len(V$p)) {
      theta <- slice_coef(V$Theta[[kk]], tt, k_i)
      lag_ex <- slice_coef(V$Lambda[[kk]], tt, k_i)
      own_lags[[kk]] <- theta
      B <- cbind(theta, lag_ex)
      if (ncol(B) != nrow(Wi)) {
        stop("B/W mismatch for ", units[[i]], " lag ", kk,
             ": ", ncol(B), " vs ", nrow(Wi))
      }
      H[[kk]] <- B %*% Wi
    }

    list(
      G = A %*% Wi,
      H = H,
      S = Sig_draw_i[[i]][tt, , ],
      local_rho = tvpgvar_spectral_radius(own_lags),
      k = k_i
    )
  })

  max_p <- max(vapply(locals, function(z) length(z$H), integer(1)))
  G <- do.call(rbind, lapply(locals, `[[`, "G"))
  H <- lapply(seq_len(max_p), function(kk) {
    do.call(rbind, lapply(locals, function(z) {
      if (kk <= length(z$H)) z$H[[kk]] else matrix(0, nrow = z$k, ncol = ncol(G))
    }))
  })
  S_post <- lapply(locals, `[[`, "S")
  local_rho <- setNames(vapply(locals, `[[`, numeric(1), "local_rho"), units)

  F <- lapply(H, function(h) solve(G, h))
  rho <- tvpgvar_spectral_radius(F)
  gc <- tvpgvar_G_condition(G)

  list(
    IRF_post = dominant_gpr_struct_irf(
      G = G, F = F, sig = S_post, x = x, units = units,
      horizon = horizon, shock_var = "GL_gpr", shock_pct = shock_pct
    ),
    max_eigen_modulus = rho,
    G_condition_number = gc$condition_number,
    G_min_singular_value = gc$min_singular_value,
    G_max_singular_value = gc$max_singular_value,
    max_abs_transition = max(vapply(F, function(z) max(abs(z)), numeric(1))),
    local_rho = local_rho
  )
}
