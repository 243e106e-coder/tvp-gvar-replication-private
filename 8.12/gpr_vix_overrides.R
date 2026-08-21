# =============================================================================
# GPR + VIX overrides for the existing structural TVP-GVAR code
#
# Identification:
#   US recursive block:
#     GPR -> VIX -> y -> dp -> r -> de -> deq
#
# Non-US current_only design:
#   contemporaneous Wex:
#     5 foreign macro + GPR + VIX
#   lagged Wex:
#     5 foreign macro + VIX
#   Direct non-US GPR(t-1) is excluded.
#
# This file MUST be sourced AFTER:
#   R/BVAR_ttvp_gprlag.r
#   R/gpr_structural_irf.R
# =============================================================================


# -----------------------------------------------------------------------------
# 1. VIX-aware Wex lag construction
# -----------------------------------------------------------------------------
tvpgvar_wex_lag <- function(
    Wex,
    lag_order = 1L,
    country_name,
    mode = c("current_only", "current_and_lag")) {

  mode <- match.arg(mode)
  lag_order <- as.integer(lag_order)

  if (
    length(lag_order) != 1L ||
      is.na(lag_order) ||
      lag_order != 1L
  ) {
    stop(
      "The GPR+VIX diagnostic requires foreign lag q=1."
    )
  }

  L <- mlag(
    Wex,
    lag_order
  )

  cc <- as.character(country_name)

  if (!identical(cc, "US")) {

    # Order:
    # 1:5 = foreign macro
    # 6   = global GPR
    # 7   = global VIX
    if (ncol(L) != 7L) {
      stop(
        "Expected seven non-US Wex lags ",
        "(5 foreign macro + global GPR + global VIX), found ",
        ncol(L),
        " for ",
        cc,
        "."
      )
    }

    if (mode == "current_only") {
      # Drop only direct lagged GPR.
      # Retain the VIX lag as a global financial-risk control.
      L <- L[
        ,
        c(1:5, 7),
        drop = FALSE
      ]
    }
  } else {

    # US GPR and VIX are endogenous domestic variables.
    # US Wex remains five trade-weighted foreign macro variables.
    if (ncol(L) != 5L) {
      stop(
        "Expected five US foreign-macro lags, found ",
        ncol(L),
        "."
      )
    }
  }

  L
}


# -----------------------------------------------------------------------------
# 2. Reconstruct omitted GPR lag as zero in the global transition system
# -----------------------------------------------------------------------------
align_lag_block <- function(
    lag_block,
    W_i,
    k_i,
    country_name) {

  lag_block <- as.matrix(lag_block)
  expected_foreign <- nrow(W_i) - k_i

  if (ncol(lag_block) == expected_foreign) {
    return(lag_block)
  }

  cc <- as.character(country_name)

  # In current_only non-US estimation:
  #
  # Estimated lag order:
  #   [foreign macro 1:5, VIX_lag]
  #
  # W_i foreign/global order:
  #   [foreign macro 1:5, GPR, VIX]
  #
  # Therefore insert a zero GPR lag coefficient before the VIX lag.
  if (
    !identical(cc, "US") &&
      expected_foreign == 7L &&
      ncol(lag_block) == 6L
  ) {

    return(
      cbind(
        lag_block[, 1:5, drop = FALSE],
        GPR_lag_zero = rep(
          0,
          nrow(lag_block)
        ),
        lag_block[, 6, drop = FALSE]
      )
    )
  }

  stop(
    "Lag-block/W mismatch for ",
    cc,
    ": lag columns=",
    ncol(lag_block),
    ", expected foreign columns=",
    expected_foreign,
    ", W rows=",
    nrow(W_i),
    ", own variables=",
    k_i
  )
}


# -----------------------------------------------------------------------------
# 3. GPR structural shock with a 7x7 US covariance block
# -----------------------------------------------------------------------------
gpr_struct_irf <- function(
    maxlag,
    G,
    F,
    sig,
    x,
    countries,
    horizon = 12,
    shock_var = "US_gpr",
    normalization = c("pct", "sd"),
    shock_pct = 10) {

  normalization <- match.arg(normalization)

  K <- nrow(x)

  if (
    nrow(G) != K ||
      ncol(G) != K
  ) {
    stop(
      "Global contemporaneous G matrix is not square."
    )
  }

  us_idx <- match(
    "US",
    countries
  )

  if (is.na(us_idx)) {
    stop("US block not found.")
  }

  if (length(sig) != length(countries)) {
    stop(
      "Country covariance list and country list have different lengths."
    )
  }

  S_us <- stabilize_cov(
    sig[[us_idx]]
  )

  if (nrow(S_us) != 7L) {
    stop(
      "US covariance block must be 7x7 under ",
      "GPR -> VIX -> y -> dp -> r -> de -> deq ordering. Found ",
      nrow(S_us),
      "x",
      ncol(S_us),
      "."
    )
  }

  # Sigma_US = L L'
  # With GPR ordered first, column 1 is the recursively identified GPR shock.
  L_us <- t(
    chol(S_us)
  )

  shock_us <- L_us[, 1]

  blocks <- lapply(
    seq_along(sig),
    function(i) {
      if (i == us_idx) {
        shock_us
      } else {
        rep(
          0,
          nrow(sig[[i]])
        )
      }
    }
  )

  u_stack <- unlist(
    blocks,
    use.names = FALSE
  )

  if (length(u_stack) != K) {
    stop(
      "Stacked innovation dimension does not match global system."
    )
  }

  impact_raw <- as.numeric(
    solve(
      G,
      u_stack
    )
  )

  j <- match(
    shock_var,
    rownames(x)
  )

  if (is.na(j)) {
    stop(
      "Shock variable not found in global vector: ",
      shock_var
    )
  }

  if (
    !is.finite(impact_raw[j]) ||
      abs(impact_raw[j]) < 1e-12
  ) {
    stop(
      "Near-zero GPR impact; cannot normalize."
    )
  }

  if (normalization == "pct") {

    target <- log1p(
      shock_pct / 100
    )

    impact <- impact_raw * (
      target / impact_raw[j]
    )

  } else {

    impact <- if (
      impact_raw[j] < 0
    ) {
      -impact_raw
    } else {
      impact_raw
    }
  }

  phi <- array(
    0,
    c(
      K,
      K,
      horizon + 1L
    )
  )

  phi[, , 1] <- diag(K)

  if (horizon >= 1L) {

    for (h in seq_len(horizon)) {

      acc <- matrix(
        0,
        K,
        K
      )

      for (
        lag in seq_len(
          min(
            maxlag,
            h
          )
        )
      ) {

        acc <- acc +
          F[[lag]] %*%
          phi[, , h - lag + 1L]
      }

      phi[, , h + 1L] <- acc
    }
  }

  out <- sapply(
    seq_len(horizon + 1L),
    function(h) {
      phi[, , h] %*% impact
    }
  )

  rownames(out) <- rownames(x)
  colnames(out) <- 0:horizon

  out
}
