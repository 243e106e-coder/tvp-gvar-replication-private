#!/usr/bin/env Rscript

# =============================================================================
# TVP-GVAR COUNTERFACTUAL SWITCH-OFF DIAGNOSTIC (NO MCMC)
#
# Uses the already-saved posterior from the dominant-unit [GPR, VIX] TVP-GVAR.
# It does NOT re-estimate the model.
#
# Main questions
# --------------
# EQUITY:
#   FULL
#   EQ_CF1_zero_current_GPR
#   EQ_CF1L_zero_lag1_GPR
#   EQ_CF1A_zero_current_and_lag1_GPR
#   EQ_CF2_zero_current_foreign_equity
#   EQ_CF3_zero_current_GPR_and_foreign_equity
#   EQ_CF3A_zero_all_GPR_and_current_foreign_equity
#   EQ_CF4_zero_current_VIX
#
# ZA REER:
#   FULL
#   ZA_CF1_zero_lag1_GPR_REER
#
# Important:
#   * Original posterior draws are unchanged.
#   * Stable draws are selected from the ORIGINAL saved model and held fixed
#     across counterfactuals.
#   * Counterfactuals modify only selected Lambda0/Lambda1 coefficients before
#     rebuilding G, H and F.
#   * The structural shock is exactly the saved model's dominant-unit recursive
#     GPR shock: GL_gpr -> GL_vix, normalized to the same +GPR percentage jump.
#   * Baseline IRFs are reconstructed and checked against the saved IRF object.
#
# Default input files:
#   prior_artifact/results/predDens_dominant_gpr_vix.rda
#   prior_artifact/results/irf_dominant_gpr_vix.rda
#
# Default output:
#   results/switch_off_diagnostic/
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

posterior_file <- Sys.getenv(
  "TVPGVAR_POSTERIOR_RDA",
  "prior_artifact/results/predDens_dominant_gpr_vix.rda"
)

irf_file <- Sys.getenv(
  "TVPGVAR_IRF_RDA",
  "prior_artifact/results/irf_dominant_gpr_vix.rda"
)

out_dir <- Sys.getenv(
  "TVPGVAR_CF_OUT",
  "results/switch_off_diagnostic"
)

selected_dates <- trimws(strsplit(
  Sys.getenv(
    "TVPGVAR_CF_DATES",
    "2003Q1,2008Q3,2014Q3,2020Q1,2022Q1,2023Q4"
  ),
  ",",
  fixed = TRUE
)[[1]])
selected_dates <- selected_dates[nzchar(selected_dates)]

shock_pct <- as.numeric(Sys.getenv("TVPGVAR_GPR_SHOCK_PCT", "10"))
horizon_req <- as.integer(Sys.getenv("TVPGVAR_CF_HORIZON", "12"))

focus_h <- as.integer(trimws(strsplit(
  Sys.getenv("TVPGVAR_CF_FOCUS_HORIZONS", "0,1,4,8,12"),
  ",",
  fixed = TRUE
)[[1]]))

check_cf_stability <- tolower(
  Sys.getenv("TVPGVAR_CF_CHECK_STABILITY", "1")
) %in% c("1", "true", "yes", "y")

# 0 = use every original stable draw.  Positive values are useful for a quick
# smoke test only.  Final diagnosis should use 0.
max_draws_per_date <- as.integer(Sys.getenv("TVPGVAR_CF_MAX_DRAWS", "0"))

if (!is.finite(shock_pct) || shock_pct <= 0) {
  stop("TVPGVAR_GPR_SHOCK_PCT must be positive.")
}
if (!is.finite(horizon_req) || horizon_req < 0L) {
  stop("TVPGVAR_CF_HORIZON must be >= 0.")
}
if (is.na(max_draws_per_date) || max_draws_per_date < 0L) {
  stop("TVPGVAR_CF_MAX_DRAWS must be >= 0.")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (f in c(posterior_file, irf_file)) {
  if (!file.exists(f)) stop("Missing required saved result: ", f)
}

# -----------------------------------------------------------------------------
# Load saved posterior + IRF
# -----------------------------------------------------------------------------

loaded_post <- load(posterior_file)
if (!all(c("predDens", "Data.setup") %in% loaded_post)) {
  stop(
    "Posterior RDA must contain predDens and Data.setup. Found: ",
    paste(loaded_post, collapse = ", ")
  )
}

loaded_irf <- load(irf_file)
if (!all(c("IRF_post", "stable_mask") %in% loaded_irf)) {
  stop(
    "IRF RDA must contain IRF_post and stable_mask. Found: ",
    paste(loaded_irf, collapse = ", ")
  )
}

xglobal <- as.matrix(Data.setup$bigx)
storage.mode(xglobal) <- "double"

units <- as.character(Data.setup$countries)
country_units <- if (!is.null(Data.setup$country_units)) {
  as.character(Data.setup$country_units)
} else {
  setdiff(units, "GL")
}

# Use the exact W matrices stored with the estimated posterior.  This mirrors the
# production IRF reconstruction and avoids any ambiguity about Data.setup copies.
globalG <- lapply(predDens, `[[`, "W")
if (any(vapply(globalG, is.null, logical(1)))) {
  if (is.null(Data.setup$gW)) {
    stop("No saved W matrices found in predDens or Data.setup.")
  }
  globalG <- Data.setup$gW
}
if (is.null(names(globalG)) || !all(units %in% names(globalG))) {
  names(globalG) <- units
}

if (length(predDens) != length(units)) {
  stop("predDens/unit count mismatch.")
}
if (length(globalG) != length(units)) {
  stop("Global W/unit count mismatch.")
}
if (!"GL" %in% units) stop("Dominant unit GL not found.")
if (!"ZA" %in% country_units) stop("ZA not found in country units.")

var_names <- colnames(xglobal)
if (is.null(var_names)) stop("Data.setup$bigx lacks column names.")
if (!identical(tail(var_names, 2L), c("GL_gpr", "GL_vix"))) {
  stop("Expected final dominant variables GL_gpr, GL_vix.")
}

irf_dates <- dimnames(IRF_post)[[1]]
irf_vars <- dimnames(IRF_post)[[2]]
irf_h_names <- dimnames(IRF_post)[[3]]

if (is.null(irf_dates) || is.null(irf_vars)) {
  stop("IRF_post lacks date/variable dimnames.")
}
if (!identical(irf_vars, var_names)) {
  stop("Saved IRF variable order differs from Data.setup$bigx.")
}

saved_h <- suppressWarnings(as.integer(irf_h_names))
if (anyNA(saved_h)) saved_h <- seq_len(dim(IRF_post)[3]) - 1L

horizon <- min(horizon_req, max(saved_h))
focus_h <- sort(unique(focus_h[focus_h >= 0L & focus_h <= horizon]))
if (!length(focus_h)) stop("No requested focus horizons are available.")

selected_dates <- intersect(selected_dates, irf_dates)
if (!length(selected_dates)) stop("No requested dates exist in the saved IRF.")

n_draws <- dim(IRF_post)[4]
if (!all(dim(stable_mask) == c(length(irf_dates), n_draws))) {
  stop("stable_mask dimensions do not match IRF_post.")
}
if (any(vapply(predDens, function(z) dim(z$ALPHA)[4], integer(1)) != n_draws)) {
  stop("Posterior ALPHA draw count differs across units / from IRF_post.")
}

# -----------------------------------------------------------------------------
# Scenario definitions
# -----------------------------------------------------------------------------

equity_scenarios <- c(
  "FULL",
  "EQ_CF1_zero_current_GPR",
  "EQ_CF1L_zero_lag1_GPR",
  "EQ_CF1A_zero_current_and_lag1_GPR",
  "EQ_CF2_zero_current_foreign_equity",
  "EQ_CF3_zero_current_GPR_and_foreign_equity",
  "EQ_CF3A_zero_all_GPR_and_current_foreign_equity",
  "EQ_CF4_zero_current_VIX"
)

za_scenarios <- c(
  "FULL",
  "ZA_CF1_zero_lag1_GPR_REER"
)

all_scenarios <- unique(c(equity_scenarios, za_scenarios))

scenario_definitions <- data.frame(
  scenario = all_scenarios,
  target = c(
    "baseline",
    "equity", "equity", "equity", "equity", "equity", "equity", "equity",
    "ZA_REER"
  ),
  restriction = c(
    "No coefficient switched off.",
    "Set Lambda0(global_gpr) = 0 only in every country equity equation.",
    "Set Lambda1(global_gpr) = 0 only in every country equity equation.",
    "Set Lambda0(global_gpr) = 0 and Lambda1(global_gpr) = 0 in every country equity equation.",
    "Set Lambda0(foreign_deq) = 0 only in every country equity equation.",
    "Set Lambda0(global_gpr) = 0 and Lambda0(foreign_deq) = 0 in every country equity equation.",
    "Set Lambda0(global_gpr) = 0, Lambda1(global_gpr) = 0, and Lambda0(foreign_deq) = 0 in every country equity equation.",
    "Set Lambda0(global_vix) = 0 only in every country equity equation.",
    "Set Lambda1(global_gpr) = 0 only in the ZA REER equation."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  scenario_definitions,
  file.path(out_dir, "00_scenario_definitions.csv"),
  row.names = FALSE
)

expected_wex <- c(
  "foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq",
  "global_gpr", "global_vix"
)

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

split_alpha <- function(A) {
  if (length(dim(A)) != 3L) {
    stop("Single-draw ALPHA must be [time, coefficient, equation].")
  }

  tags <- dimnames(A)[[2]]
  eq <- dimnames(A)[[3]]
  if (is.null(tags) || is.null(eq)) {
    stop("ALPHA array lacks coefficient/equation dimnames.")
  }

  block <- function(tag, optional = FALSE) {
    ii <- which(tags == tag)
    if (!length(ii)) {
      if (!optional) stop("Missing coefficient block: ", tag)
      return(array(0, dim = c(dim(A)[1], 0L, dim(A)[3])))
    }
    A[, ii, , drop = FALSE]
  }

  lag_tags <- grep("^Ylag[0-9]+$", unique(tags), value = TRUE)
  if (!length(lag_tags)) stop("No domestic lag blocks in ALPHA.")

  lag_nums <- sort(unique(as.integer(sub("^Ylag", "", lag_tags))))
  if (anyNA(lag_nums) || !identical(lag_nums, seq_len(max(lag_nums)))) {
    stop("Domestic lag labels are not contiguous.")
  }
  p <- max(lag_nums)

  list(
    Lambda0 = block("Wex", optional = TRUE),
    Lambda = lapply(
      seq_len(p),
      function(j) block(paste0("Wexlag", j), optional = TRUE)
    ),
    Theta = lapply(
      seq_len(p),
      function(j) block(paste0("Ylag", j))
    ),
    p = p,
    eq = eq
  )
}

slice_coef <- function(arr, tt, nrow_out) {
  if (dim(arr)[2] == 0L) {
    return(matrix(0, nrow = nrow_out, ncol = 0L))
  }
  z <- arr[tt, , , drop = FALSE]
  z <- matrix(z, nrow = dim(arr)[2], ncol = nrow_out)
  t(z)
}

get_W <- function(i, cc) {
  if (!is.null(names(globalG)) && cc %in% names(globalG)) {
    return(as.matrix(globalG[[cc]]))
  }
  as.matrix(globalG[[i]])
}

eq_index <- function(eq_names, cc, suffix) {
  out <- match(paste0(cc, "_", suffix), eq_names)
  if (is.na(out)) out <- match(suffix, eq_names)
  out
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

spectral_radius <- function(F) {
  if (!length(F)) stop("No transition matrices.")
  K <- nrow(F[[1]])
  p <- length(F)

  if (p == 1L) {
    return(as.numeric(max(Mod(eigen(F[[1]], only.values = TRUE)$values))))
  }

  top <- do.call(cbind, F)
  lower <- cbind(
    diag(K * (p - 1L)),
    matrix(0, K * (p - 1L), K)
  )
  companion <- rbind(top, lower)

  as.numeric(max(Mod(eigen(companion, only.values = TRUE)$values)))
}

qsum <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(
      low90 = NA_real_,
      low68 = NA_real_,
      median = NA_real_,
      high68 = NA_real_,
      high90 = NA_real_,
      mean = NA_real_,
      positive_share = NA_real_,
      negative_share = NA_real_,
      n = 0
    ))
  }

  qq <- quantile(x, c(.05, .16, .50, .84, .95), names = FALSE)
  c(
    low90 = qq[1],
    low68 = qq[2],
    median = qq[3],
    high68 = qq[4],
    high90 = qq[5],
    mean = mean(x),
    positive_share = mean(x > 0),
    negative_share = mean(x < 0),
    n = length(x)
  )
}

# Apply scenario restrictions to a unit's contemporaneous and lagged
# weakly-exogenous coefficient matrices.
mutate_switches <- function(
  lam0,
  lam_lags,
  cc,
  eq_names,
  wex_names,
  scenario
) {
  if (scenario == "FULL" || cc == "GL") {
    return(list(lam0 = lam0, lam_lags = lam_lags))
  }

  if (!identical(wex_names, expected_wex)) {
    stop(
      "Unexpected Wex order for ", cc, ": ",
      paste(wex_names, collapse = ", ")
    )
  }

  gpr_j <- match("global_gpr", wex_names)
  vix_j <- match("global_vix", wex_names)
  feq_j <- match("foreign_deq", wex_names)

  deq_i <- eq_index(eq_names, cc, "deq")
  de_i <- eq_index(eq_names, cc, "de")

  if (grepl("^EQ_", scenario)) {
    if (is.na(deq_i)) stop("Equity equation not found for ", cc)

    if (scenario == "EQ_CF1_zero_current_GPR") {
      lam0[deq_i, gpr_j] <- 0
    } else if (scenario == "EQ_CF1L_zero_lag1_GPR") {
      if (!length(lam_lags) || ncol(lam_lags[[1]]) == 0L) {
        stop("Lag-1 Wex block missing for ", cc)
      }
      lam_lags[[1]][deq_i, gpr_j] <- 0
    } else if (scenario == "EQ_CF1A_zero_current_and_lag1_GPR") {
      lam0[deq_i, gpr_j] <- 0
      if (!length(lam_lags) || ncol(lam_lags[[1]]) == 0L) {
        stop("Lag-1 Wex block missing for ", cc)
      }
      lam_lags[[1]][deq_i, gpr_j] <- 0
    } else if (scenario == "EQ_CF2_zero_current_foreign_equity") {
      lam0[deq_i, feq_j] <- 0
    } else if (scenario == "EQ_CF3_zero_current_GPR_and_foreign_equity") {
      lam0[deq_i, gpr_j] <- 0
      lam0[deq_i, feq_j] <- 0
    } else if (scenario == "EQ_CF3A_zero_all_GPR_and_current_foreign_equity") {
      lam0[deq_i, gpr_j] <- 0
      lam0[deq_i, feq_j] <- 0
      if (!length(lam_lags) || ncol(lam_lags[[1]]) == 0L) {
        stop("Lag-1 Wex block missing for ", cc)
      }
      lam_lags[[1]][deq_i, gpr_j] <- 0
    } else if (scenario == "EQ_CF4_zero_current_VIX") {
      lam0[deq_i, vix_j] <- 0
    } else {
      stop("Unknown equity scenario: ", scenario)
    }
  }

  if (scenario == "ZA_CF1_zero_lag1_GPR_REER" && cc == "ZA") {
    if (is.na(de_i)) stop("ZA REER equation not found.")
    if (!length(lam_lags) || ncol(lam_lags[[1]]) == 0L) {
      stop("Lag-1 Wex block missing for ZA.")
    }
    lam_lags[[1]][de_i, gpr_j] <- 0
  }

  list(lam0 = lam0, lam_lags = lam_lags)
}

# Rebuild G, H, F from the same posterior draw after applying one switch-off.
build_counterfactual_state <- function(tt, draw, scenario) {
  locals <- vector("list", length(units))

  for (i in seq_along(units)) {
    cc <- units[[i]]

    A4 <- predDens[[i]]$ALPHA
    A_draw <- A4[, , , draw, drop = FALSE]
    A_draw <- array(
      A_draw,
      dim = dim(A4)[1:3],
      dimnames = dimnames(A4)[1:3]
    )

    V <- split_alpha(A_draw)
    k_i <- length(V$eq)
    Wi <- get_W(i, cc)

    if (nrow(Wi) < k_i) {
      stop("W has fewer rows than own variables for ", cc)
    }

    expected_foreign <- nrow(Wi) - k_i
    if (cc == "GL") {
      if (expected_foreign != 0L) {
        stop("GL should have zero weakly-exogenous rows.")
      }
      wex_names <- character(0)
    } else {
      if (expected_foreign != 7L) {
        stop(
          "Expected 7 Wex rows for ", cc,
          "; found ", expected_foreign
        )
      }
      if (is.null(rownames(Wi))) {
        stop("W matrix lacks row names for ", cc)
      }
      wex_names <- rownames(Wi)[(k_i + 1L):nrow(Wi)]
      if (!identical(wex_names, expected_wex)) {
        stop(
          "Wex row order mismatch for ", cc, ". Found: ",
          paste(wex_names, collapse = ", ")
        )
      }
    }

    lam0 <- slice_coef(V$Lambda0, tt, k_i)
    lam_lags <- lapply(
      seq_len(V$p),
      function(kk) slice_coef(V$Lambda[[kk]], tt, k_i)
    )

    # For p=2, q=1, Wexlag2 is structurally zero.
    for (kk in seq_len(V$p)) {
      if (ncol(lam_lags[[kk]]) == 0L && expected_foreign > 0L) {
        if (kk <= 1L) {
          stop("Missing Wexlag1 block for ", cc)
        }
        lam_lags[[kk]] <- matrix(
          0,
          nrow = k_i,
          ncol = expected_foreign
        )
      }

      if (ncol(lam_lags[[kk]]) != expected_foreign) {
        stop(
          "Foreign-lag/W mismatch for ", cc,
          " lag ", kk,
          ": got ", ncol(lam_lags[[kk]]),
          ", expected ", expected_foreign
        )
      }
    }

    sw <- mutate_switches(
      lam0 = lam0,
      lam_lags = lam_lags,
      cc = cc,
      eq_names = V$eq,
      wex_names = wex_names,
      scenario = scenario
    )
    lam0 <- sw$lam0
    lam_lags <- sw$lam_lags

    A0 <- cbind(diag(k_i), -lam0)
    if (ncol(A0) != nrow(Wi)) {
      stop("A0/W mismatch for ", cc, " under ", scenario)
    }

    H_i <- vector("list", V$p)
    for (kk in seq_len(V$p)) {
      theta <- slice_coef(V$Theta[[kk]], tt, k_i)
      B <- cbind(theta, lam_lags[[kk]])

      if (ncol(B) != nrow(Wi)) {
        stop(
          "B/W mismatch for ", cc,
          " lag ", kk, " under ", scenario
        )
      }

      H_i[[kk]] <- B %*% Wi
    }

    S4 <- predDens[[i]]$SIGMApost
    S <- S4[tt, , , draw, drop = FALSE]
    S <- matrix(S, nrow = k_i, ncol = k_i)

    locals[[i]] <- list(
      G = A0 %*% Wi,
      H = H_i,
      S = S,
      k = k_i
    )
  }

  max_p <- max(vapply(locals, function(z) length(z$H), integer(1)))
  G <- do.call(rbind, lapply(locals, `[[`, "G"))

  H <- lapply(seq_len(max_p), function(kk) {
    do.call(
      rbind,
      lapply(locals, function(z) {
        if (kk <= length(z$H)) {
          z$H[[kk]]
        } else {
          matrix(0, nrow = z$k, ncol = ncol(G))
        }
      })
    )
  })

  F <- lapply(H, function(h) solve(G, h))
  S_post <- lapply(locals, `[[`, "S")

  rho <- if (check_cf_stability) spectral_radius(F) else NA_real_

  list(
    G = G,
    F = F,
    S = S_post,
    rho = rho
  )
}

# Same structural identification as the saved dominant-unit model.
make_structural_irf <- function(state, horizon, shock_pct) {
  G <- state$G
  F <- state$F
  S_post <- state$S

  K <- ncol(xglobal)
  if (!all(dim(G) == c(K, K))) {
    stop("Global G is not K x K.")
  }

  gl_i <- match("GL", units)
  S_gl <- stabilize_cov(S_post[[gl_i]])
  if (!all(dim(S_gl) == c(2L, 2L))) {
    stop("Dominant covariance must be 2 x 2.")
  }

  # Recursive identification in dominant block: GPR -> VIX.
  L <- t(chol(S_gl))
  shock_gl <- L[, 1]

  u <- unlist(
    lapply(seq_along(units), function(i) {
      k_i <- nrow(S_post[[i]])
      if (i == gl_i) shock_gl else rep(0, k_i)
    }),
    use.names = FALSE
  )

  if (length(u) != K) stop("Innovation stack dimension mismatch.")

  impact_raw <- as.numeric(solve(G, u))
  gpr_j <- match("GL_gpr", var_names)
  if (is.na(gpr_j) ||
      !is.finite(impact_raw[gpr_j]) ||
      abs(impact_raw[gpr_j]) < 1e-12) {
    stop("Cannot normalize GL_gpr impact.")
  }

  target <- log1p(shock_pct / 100)
  impact <- impact_raw * (target / impact_raw[gpr_j])

  out <- matrix(
    0,
    nrow = K,
    ncol = horizon + 1L,
    dimnames = list(var_names, as.character(0:horizon))
  )
  out[, 1L] <- impact

  if (horizon >= 1L) {
    for (h in seq_len(horizon)) {
      for (lag in seq_len(min(length(F), h))) {
        out[, h + 1L] <- out[, h + 1L] +
          F[[lag]] %*% out[, h - lag + 1L]
      }
    }
  }

  out
}

# -----------------------------------------------------------------------------
# Allocate compact draw arrays
# -----------------------------------------------------------------------------

nd <- length(selected_dates)
ne <- length(equity_scenarios)
nz <- length(za_scenarios)
nc <- length(country_units)
nhf <- length(focus_h)

eq_direct <- array(
  NA_real_,
  dim = c(nd, ne, nc, nhf, n_draws),
  dimnames = list(
    date = selected_dates,
    scenario = equity_scenarios,
    country = country_units,
    horizon = as.character(focus_h),
    draw = NULL
  )
)

eq_cumulative <- eq_direct

za_direct <- array(
  NA_real_,
  dim = c(nd, nz, nhf, n_draws),
  dimnames = list(
    date = selected_dates,
    scenario = za_scenarios,
    horizon = as.character(focus_h),
    draw = NULL
  )
)

za_cumulative <- za_direct

rho_arr <- array(
  NA_real_,
  dim = c(nd, length(all_scenarios), n_draws),
  dimnames = list(
    date = selected_dates,
    scenario = all_scenarios,
    draw = NULL
  )
)

solve_ok <- array(
  FALSE,
  dim = c(nd, length(all_scenarios), n_draws),
  dimnames = list(
    date = selected_dates,
    scenario = all_scenarios,
    draw = NULL
  )
)

used_draw <- matrix(
  FALSE,
  nrow = nd,
  ncol = n_draws,
  dimnames = list(selected_dates, NULL)
)

max_baseline_reproduction_error <- 0

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

cat("\n============================================================\n")
cat(" TVP-GVAR SWITCH-OFF DIAGNOSTIC — NO MCMC\n")
cat(" Posterior: ", posterior_file, "\n", sep = "")
cat(" IRF:       ", irf_file, "\n", sep = "")
cat(" Dates:     ", paste(selected_dates, collapse = ", "), "\n", sep = "")
cat(" Horizon:   0-", horizon, "\n", sep = "")
cat(" Focus h:   ", paste(focus_h, collapse = ", "), "\n", sep = "")
cat(" Shock:     +", shock_pct, "% GL_gpr\n", sep = "")
cat(" CF stability check: ", check_cf_stability, "\n", sep = "")
cat("============================================================\n\n")

for (d_i in seq_along(selected_dates)) {
  ddate <- selected_dates[[d_i]]
  tt <- match(ddate, irf_dates)

  draws <- which(stable_mask[tt, ])
  if (!length(draws)) {
    warning("No original stable draws at ", ddate, "; skipping.")
    next
  }

  if (max_draws_per_date > 0L && length(draws) > max_draws_per_date) {
    draws <- draws[seq_len(max_draws_per_date)]
  }

  used_draw[d_i, draws] <- TRUE

  cat(
    "Date ", ddate, ": using ",
    length(draws), "/", sum(stable_mask[tt, ]),
    " original stable draws\n",
    sep = ""
  )

  for (jj in seq_along(draws)) {
    draw <- draws[[jj]]

    # FULL state + exact saved-IRF reproduction check.
    base_state <- build_counterfactual_state(
      tt = tt,
      draw = draw,
      scenario = "FULL"
    )
    base_path <- make_structural_irf(
      base_state,
      horizon = horizon,
      shock_pct = shock_pct
    )

    saved <- IRF_post[tt, , seq_len(horizon + 1L), draw, drop = FALSE]
    saved <- matrix(
      saved,
      nrow = length(var_names),
      ncol = horizon + 1L,
      dimnames = list(var_names, as.character(0:horizon))
    )

    err <- max(abs(base_path - saved), na.rm = TRUE)
    max_baseline_reproduction_error <- max(
      max_baseline_reproduction_error,
      err
    )

    if (!is.finite(err) || err > 1e-7) {
      stop(
        "Baseline reconstruction mismatch at ", ddate,
        ", draw ", draw,
        ". max abs error = ", signif(err, 8)
      )
    }

    # Cache FULL once.
    states <- list(FULL = base_state)
    paths <- list(FULL = base_path)

    # Build remaining counterfactual states.
    for (scenario in setdiff(all_scenarios, "FULL")) {
      ans <- tryCatch(
        {
          st <- build_counterfactual_state(
            tt = tt,
            draw = draw,
            scenario = scenario
          )
          pa <- make_structural_irf(
            st,
            horizon = horizon,
            shock_pct = shock_pct
          )
          list(ok = TRUE, state = st, path = pa)
        },
        error = function(e) {
          warning(
            "Counterfactual failed: date=", ddate,
            ", draw=", draw,
            ", scenario=", scenario,
            " | ", conditionMessage(e)
          )
          list(ok = FALSE, state = NULL, path = NULL)
        }
      )

      if (ans$ok) {
        states[[scenario]] <- ans$state
        paths[[scenario]] <- ans$path
      }
    }

    # State diagnostics.
    for (scenario in names(paths)) {
      s_j <- match(scenario, all_scenarios)
      solve_ok[d_i, s_j, draw] <- TRUE
      rho_arr[d_i, s_j, draw] <- states[[scenario]]$rho
    }

    # Equity: direct change IRF and draw-by-draw cumulative level effect.
    for (scenario in equity_scenarios) {
      if (is.null(paths[[scenario]])) next

      s_j <- match(scenario, equity_scenarios)
      path <- paths[[scenario]]

      for (c_i in seq_along(country_units)) {
        cc <- country_units[[c_i]]
        vv <- paste0(cc, "_deq")
        v_i <- match(vv, var_names)
        if (is.na(v_i)) stop("Equity variable not found: ", vv)

        z <- as.numeric(path[v_i, ])
        cz <- cumsum(z)

        for (h_i in seq_along(focus_h)) {
          h <- focus_h[[h_i]]
          eq_direct[d_i, s_j, c_i, h_i, draw] <- z[h + 1L]
          eq_cumulative[d_i, s_j, c_i, h_i, draw] <- cz[h + 1L]
        }
      }
    }

    # ZA REER.
    za_i <- match("ZA_de", var_names)
    if (is.na(za_i)) stop("ZA_de not found in global variables.")

    for (scenario in za_scenarios) {
      if (is.null(paths[[scenario]])) next

      s_j <- match(scenario, za_scenarios)
      z <- as.numeric(paths[[scenario]][za_i, ])
      cz <- cumsum(z)

      for (h_i in seq_along(focus_h)) {
        h <- focus_h[[h_i]]
        za_direct[d_i, s_j, h_i, draw] <- z[h + 1L]
        za_cumulative[d_i, s_j, h_i, draw] <- cz[h + 1L]
      }
    }

    if (jj %% 25L == 0L || jj == length(draws)) {
      cat(
        "  ", ddate, ": draw ",
        jj, "/", length(draws), " complete\n",
        sep = ""
      )
    }
  }
}

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

summarise_equity_array <- function(arr, representation) {
  rows <- list()
  rr <- 0L

  for (d_i in seq_len(nd)) {
    for (s_i in seq_len(ne)) {
      for (c_i in seq_len(nc)) {
        for (h_i in seq_len(nhf)) {
          qs <- qsum(arr[d_i, s_i, c_i, h_i, ])

          rr <- rr + 1L
          rows[[rr]] <- data.frame(
            date = selected_dates[d_i],
            scenario = equity_scenarios[s_i],
            country = country_units[c_i],
            variable = paste0(country_units[c_i], "_deq"),
            representation = representation,
            horizon = focus_h[h_i],
            median = qs[["median"]],
            mean = qs[["mean"]],
            low68 = qs[["low68"]],
            high68 = qs[["high68"]],
            low90 = qs[["low90"]],
            high90 = qs[["high90"]],
            positive_share = qs[["positive_share"]],
            negative_share = qs[["negative_share"]],
            draws = as.integer(qs[["n"]]),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  do.call(rbind, rows)
}

summarise_za_array <- function(arr, representation) {
  rows <- list()
  rr <- 0L

  for (d_i in seq_len(nd)) {
    for (s_i in seq_len(nz)) {
      for (h_i in seq_len(nhf)) {
        qs <- qsum(arr[d_i, s_i, h_i, ])

        rr <- rr + 1L
        rows[[rr]] <- data.frame(
          date = selected_dates[d_i],
          scenario = za_scenarios[s_i],
          country = "ZA",
          variable = "ZA_de",
          representation = representation,
          horizon = focus_h[h_i],
          median = qs[["median"]],
          mean = qs[["mean"]],
          low68 = qs[["low68"]],
          high68 = qs[["high68"]],
          low90 = qs[["low90"]],
          high90 = qs[["high90"]],
          positive_share = qs[["positive_share"]],
          negative_share = qs[["negative_share"]],
          draws = as.integer(qs[["n"]]),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  do.call(rbind, rows)
}

equity_summary <- rbind(
  summarise_equity_array(eq_direct, "direct_change_IRF"),
  summarise_equity_array(eq_cumulative, "cumulative_level_effect")
)

za_summary <- rbind(
  summarise_za_array(za_direct, "direct_change_IRF"),
  summarise_za_array(za_cumulative, "cumulative_level_effect")
)

write.csv(
  equity_summary,
  file.path(out_dir, "01_equity_counterfactual_posterior_summary.csv"),
  row.names = FALSE
)

write.csv(
  za_summary,
  file.path(out_dir, "02_ZA_REER_counterfactual_posterior_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Equity sign counts and baseline comparisons
# -----------------------------------------------------------------------------

eq_cum_sum <- subset(
  equity_summary,
  representation == "cumulative_level_effect"
)

sign_rows <- list()
rr <- 0L

for (idx in split(
  seq_len(nrow(eq_cum_sum)),
  interaction(
    eq_cum_sum$date,
    eq_cum_sum$scenario,
    eq_cum_sum$horizon,
    drop = TRUE,
    lex.order = TRUE
  )
)) {
  d <- eq_cum_sum[idx, , drop = FALSE]
  rr <- rr + 1L

  sign_rows[[rr]] <- data.frame(
    date = d$date[1],
    scenario = d$scenario[1],
    horizon = d$horizon[1],
    countries = nrow(d),
    positive_median_countries = sum(d$median > 0, na.rm = TRUE),
    negative_median_countries = sum(d$median < 0, na.rm = TRUE),
    significant_positive_68 = sum(d$low68 > 0, na.rm = TRUE),
    significant_negative_68 = sum(d$high68 < 0, na.rm = TRUE),
    median_abs_response = median(abs(d$median), na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

equity_sign_counts <- do.call(rbind, sign_rows)

equity_sign_counts <- equity_sign_counts[
  order(
    match(equity_sign_counts$date, selected_dates),
    match(equity_sign_counts$scenario, equity_scenarios),
    equity_sign_counts$horizon
  ),
  ,
  drop = FALSE
]

write.csv(
  equity_sign_counts,
  file.path(out_dir, "03_equity_sign_counts_by_scenario.csv"),
  row.names = FALSE
)

baseline_eq <- eq_cum_sum[
  eq_cum_sum$scenario == "FULL",
  c("date", "country", "horizon", "median", "low68", "high68")
]
names(baseline_eq)[4:6] <- c(
  "baseline_median",
  "baseline_low68",
  "baseline_high68"
)

eq_compare <- merge(
  eq_cum_sum,
  baseline_eq,
  by = c("date", "country", "horizon"),
  all.x = TRUE,
  sort = FALSE
)

eq_compare$delta_median_vs_full <-
  eq_compare$median - eq_compare$baseline_median

eq_compare$ratio_to_full <- ifelse(
  abs(eq_compare$baseline_median) > 1e-12,
  eq_compare$median / eq_compare$baseline_median,
  NA_real_
)

eq_compare$baseline_positive <- eq_compare$baseline_median > 0
eq_compare$counterfactual_positive <- eq_compare$median > 0

eq_compare$positive_to_nonpositive_flip <-
  eq_compare$baseline_positive & !eq_compare$counterfactual_positive

eq_compare$sign_changed_vs_full <-
  sign(eq_compare$median) != sign(eq_compare$baseline_median)

eq_compare <- eq_compare[
  order(
    match(eq_compare$date, selected_dates),
    match(eq_compare$scenario, equity_scenarios),
    match(eq_compare$country, country_units),
    eq_compare$horizon
  ),
  ,
  drop = FALSE
]

write.csv(
  eq_compare,
  file.path(out_dir, "04_equity_scenario_vs_FULL.csv"),
  row.names = FALSE
)

# Compact date/scenario decision table.
decision_rows <- list()
rr <- 0L

for (idx in split(
  seq_len(nrow(eq_compare)),
  interaction(
    eq_compare$date,
    eq_compare$scenario,
    eq_compare$horizon,
    drop = TRUE,
    lex.order = TRUE
  )
)) {
  d <- eq_compare[idx, , drop = FALSE]
  nonzero_base <- is.finite(d$baseline_median) & abs(d$baseline_median) > 1e-12

  rr <- rr + 1L
  decision_rows[[rr]] <- data.frame(
    date = d$date[1],
    scenario = d$scenario[1],
    horizon = d$horizon[1],
    countries = nrow(d),
    baseline_positive_countries = sum(d$baseline_median > 0, na.rm = TRUE),
    counterfactual_positive_countries = sum(d$median > 0, na.rm = TRUE),
    positive_to_nonpositive_flips = sum(
      d$positive_to_nonpositive_flip,
      na.rm = TRUE
    ),
    any_sign_flips = sum(d$sign_changed_vs_full, na.rm = TRUE),
    median_abs_full = median(abs(d$baseline_median), na.rm = TRUE),
    median_abs_counterfactual = median(abs(d$median), na.rm = TRUE),
    median_ratio_to_full = if (any(nonzero_base)) {
      median(d$ratio_to_full[nonzero_base], na.rm = TRUE)
    } else {
      NA_real_
    },
    stringsAsFactors = FALSE
  )
}

equity_decision <- do.call(rbind, decision_rows)
equity_decision <- equity_decision[
  order(
    match(equity_decision$date, selected_dates),
    match(equity_decision$scenario, equity_scenarios),
    equity_decision$horizon
  ),
  ,
  drop = FALSE
]

write.csv(
  equity_decision,
  file.path(out_dir, "05_equity_switch_off_decision_table.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# ZA comparison
# -----------------------------------------------------------------------------

za_cum_sum <- subset(
  za_summary,
  representation == "cumulative_level_effect"
)

za_base <- za_cum_sum[
  za_cum_sum$scenario == "FULL",
  c("date", "horizon", "median", "low68", "high68")
]
names(za_base)[3:5] <- c(
  "baseline_median",
  "baseline_low68",
  "baseline_high68"
)

za_compare <- merge(
  za_cum_sum,
  za_base,
  by = c("date", "horizon"),
  all.x = TRUE,
  sort = FALSE
)

za_compare$delta_median_vs_full <-
  za_compare$median - za_compare$baseline_median

za_compare$ratio_to_full <- ifelse(
  abs(za_compare$baseline_median) > 1e-12,
  za_compare$median / za_compare$baseline_median,
  NA_real_
)

za_compare$sign_changed_vs_full <-
  sign(za_compare$median) != sign(za_compare$baseline_median)

za_compare <- za_compare[
  order(
    match(za_compare$date, selected_dates),
    match(za_compare$scenario, za_scenarios),
    za_compare$horizon
  ),
  ,
  drop = FALSE
]

write.csv(
  za_compare,
  file.path(out_dir, "06_ZA_REER_scenario_vs_FULL.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Counterfactual solve/stability diagnostics
# -----------------------------------------------------------------------------

stab_rows <- list()
rr <- 0L

for (d_i in seq_len(nd)) {
  for (s_i in seq_along(all_scenarios)) {
    used <- which(used_draw[d_i, ])
    if (!length(used)) next

    ok <- solve_ok[d_i, s_i, used]
    rho <- rho_arr[d_i, s_i, used]
    rho_ok <- rho[ok & is.finite(rho)]

    rr <- rr + 1L
    stab_rows[[rr]] <- data.frame(
      date = selected_dates[d_i],
      scenario = all_scenarios[s_i],
      original_stable_draws_used = length(used),
      counterfactual_solve_success = sum(ok),
      counterfactual_solve_success_share = mean(ok),
      stability_checked = check_cf_stability,
      counterfactual_stable_draws = if (check_cf_stability) {
        sum(rho_ok < 1)
      } else {
        NA_integer_
      },
      counterfactual_stable_share_within_success = if (
        check_cf_stability && length(rho_ok)
      ) {
        mean(rho_ok < 1)
      } else {
        NA_real_
      },
      median_counterfactual_rho = if (
        check_cf_stability && length(rho_ok)
      ) {
        median(rho_ok)
      } else {
        NA_real_
      },
      p95_counterfactual_rho = if (
        check_cf_stability && length(rho_ok)
      ) {
        as.numeric(quantile(rho_ok, .95, names = FALSE))
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  }
}

stability_summary <- do.call(rbind, stab_rows)

write.csv(
  stability_summary,
  file.path(out_dir, "07_counterfactual_solve_stability_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Focus tables
# -----------------------------------------------------------------------------

diagnostic_date <- if ("2022Q1" %in% selected_dates) {
  "2022Q1"
} else {
  tail(selected_dates, 1)
}

h_focus <- if (12L %in% focus_h) 12L else max(focus_h)

focus_equity <- eq_compare[
  eq_compare$date == diagnostic_date &
    eq_compare$horizon == h_focus,
  ,
  drop = FALSE
]

focus_za <- za_compare[
  za_compare$date == diagnostic_date &
    za_compare$horizon == h_focus,
  ,
  drop = FALSE
]

write.csv(
  focus_equity,
  file.path(
    out_dir,
    paste0("FOCUS_", diagnostic_date, "_h", h_focus, "_equity_switch_off.csv")
  ),
  row.names = FALSE
)

write.csv(
  focus_za,
  file.path(
    out_dir,
    paste0("FOCUS_", diagnostic_date, "_h", h_focus, "_ZA_REER_switch_off.csv")
  ),
  row.names = FALSE
)

# Save compact draw arrays for later inspection without giant CSV files.
save(
  eq_direct,
  eq_cumulative,
  za_direct,
  za_cumulative,
  rho_arr,
  solve_ok,
  used_draw,
  selected_dates,
  equity_scenarios,
  za_scenarios,
  focus_h,
  file = file.path(out_dir, "switch_off_draw_arrays.rda")
)

# -----------------------------------------------------------------------------
# Automated human-readable diagnosis
# -----------------------------------------------------------------------------

decision_at <- function(scenario) {
  d <- equity_decision[
    equity_decision$date == diagnostic_date &
      equity_decision$horizon == h_focus &
      equity_decision$scenario == scenario,
    ,
    drop = FALSE
  ]
  if (!nrow(d)) return(NULL)
  d[1, , drop = FALSE]
}

fmt_num <- function(x, digits = 4) {
  if (!length(x) || !is.finite(x)) return("NA")
  format(round(x, digits), nsmall = digits, trim = TRUE)
}

scenario_line <- function(scenario) {
  d <- decision_at(scenario)
  if (is.null(d)) return(paste0(scenario, ": unavailable"))

  paste0(
    scenario,
    " | positive medians=", d$counterfactual_positive_countries, "/",
    d$countries,
    " | positive->nonpositive flips=", d$positive_to_nonpositive_flips,
    " | median |response| ratio vs FULL=", fmt_num(d$median_ratio_to_full)
  )
}

za_full <- focus_za[focus_za$scenario == "FULL", , drop = FALSE]
za_off <- focus_za[
  focus_za$scenario == "ZA_CF1_zero_lag1_GPR_REER",
  ,
  drop = FALSE
]

report <- c(
  "TVP-GVAR COUNTERFACTUAL SWITCH-OFF DIAGNOSTIC — NO MCMC",
  "=======================================================",
  paste0("Posterior: ", posterior_file),
  paste0("IRF: ", irf_file),
  paste0("Dates: ", paste(selected_dates, collapse = ", ")),
  paste0("Shock: +", shock_pct, "% GL_gpr"),
  paste0("Maximum FULL-vs-saved IRF reconstruction error: ",
         format(max_baseline_reproduction_error, scientific = TRUE, digits = 6)),
  paste0("Original stable draws held fixed across counterfactuals: YES"),
  paste0("Counterfactual stability recomputed: ", check_cf_stability),
  "",
  paste0("EQUITY FOCUS: ", diagnostic_date, ", cumulative h=", h_focus),
  scenario_line("FULL"),
  scenario_line("EQ_CF1_zero_current_GPR"),
  scenario_line("EQ_CF1L_zero_lag1_GPR"),
  scenario_line("EQ_CF1A_zero_current_and_lag1_GPR"),
  scenario_line("EQ_CF2_zero_current_foreign_equity"),
  scenario_line("EQ_CF3_zero_current_GPR_and_foreign_equity"),
  scenario_line("EQ_CF3A_zero_all_GPR_and_current_foreign_equity"),
  scenario_line("EQ_CF4_zero_current_VIX"),
  "",
  "INTERPRETATION RULES",
  "1) If zeroing current GPR flips many positive equity medians, the conditional current-GPR loading is a root source.",
  "2) If zeroing lagged GPR has little effect but zeroing current GPR has a large effect, the sign problem is mainly contemporaneous/conditional rather than dynamic.",
  "3) If zeroing foreign equity mainly shrinks magnitude but leaves signs positive, foreign equity is an amplifier rather than the first mover.",
  "4) If removing current GPR + foreign equity produces negative responses while FULL is positive, the positive direct-GPR/network combination dominates the negative risk channel.",
  "5) If zeroing VIX makes equity more positive, VIX is acting as a negative counterweight, not the source of the positive IRF.",
  "",
  paste0("ZA REER FOCUS: ", diagnostic_date, ", cumulative h=", h_focus),
  if (nrow(za_full)) {
    paste0(
      "FULL median=", fmt_num(za_full$median[1]),
      " | 68% CI=[", fmt_num(za_full$low68[1]), ", ",
      fmt_num(za_full$high68[1]), "]"
    )
  } else {
    "FULL ZA result unavailable"
  },
  if (nrow(za_off)) {
    paste0(
      "Lagged-GPR-off median=", fmt_num(za_off$median[1]),
      " | 68% CI=[", fmt_num(za_off$low68[1]), ", ",
      fmt_num(za_off$high68[1]), "]",
      " | ratio vs FULL=", fmt_num(za_off$ratio_to_full[1])
    )
  } else {
    "ZA lagged-GPR-off result unavailable"
  },
  "If ZA REER collapses toward zero after Lambda1(global_gpr) is switched off, the large ZA response is generated by that lagged GPR loading rather than REER scaling.",
  "",
  "PRIMARY FILES TO READ",
  "05_equity_switch_off_decision_table.csv",
  paste0("FOCUS_", diagnostic_date, "_h", h_focus, "_equity_switch_off.csv"),
  paste0("FOCUS_", diagnostic_date, "_h", h_focus, "_ZA_REER_switch_off.csv"),
  "07_counterfactual_solve_stability_summary.csv",
  "",
  "This is a post-estimation counterfactual diagnostic, not a newly estimated structural model."
)

writeLines(
  report,
  file.path(out_dir, "README_switch_off_diagnosis.txt")
)

cat("\n")
cat(paste(report, collapse = "\n"))
cat("\n\nDone. Outputs written to: ", out_dir, "\n", sep = "")
