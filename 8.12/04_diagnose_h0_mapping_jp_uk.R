#!/usr/bin/env Rscript

# =============================================================================
# Lightweight post-estimation diagnostic for the P2-Strong GPR + VIX TVP-GVAR
#
# PURPOSE
#   A. Trace the contemporaneous (h=0) GPR -> GDP mapping:
#        US covariance -> Cholesky column 1 -> stacked innovation u
#        -> solve(G, u) -> global h=0 impact
#
#   B. Compare the current recursive GPR shock with a counterfactual
#      "GPR-only residual innovation":
#        - Full recursive: keep the complete first Cholesky column
#        - GPR-only: keep only the US_gpr element of that column, set the
#          contemporaneous VIX/GDP/inflation/rate/FX/equity loadings to zero
#        - Normalize BOTH versions to the same +10% US_gpr impact
#
#      This tells us whether the systematic positive h=0 GDP response originates
#      mainly from the US Cholesky covariance loading or from the global G^{-1}
#      contemporaneous mapping.
#
#   C. Diagnose JP/UK local instability without re-running MCMC:
#        - reconstruct domestic lag matrices Theta_1, Theta_2
#        - calculate local companion spectral radius
#        - inspect own AR(1), AR(2), AR(1)+AR(2)
#        - zero one equation row at a time as a diagnostic counterfactual
#          and measure the reduction in local spectral radius
#
# INPUTS (from the already completed GPR+VIX artifact)
#   results/predDens_gpr_vix_structural.rda
#   results/irf_gpr_structural.rda
#
# OUTPUTS
#   results/h0_diag/h0_full_vs_gpr_only_draws.csv
#   results/h0_diag/h0_gdp_full_vs_gpr_only_summary.csv
#   results/h0_diag/h0_global_impact_summary.csv
#   results/h0_diag/us_recursive_gpr_loading_draws.csv
#   results/h0_diag/us_recursive_gpr_loading_summary.csv
#   results/h0_diag/jp_uk_stability_driver_draws.csv
#   results/h0_diag/jp_uk_stability_driver_summary.csv
#   results/h0_diag/jp_uk_stability_driver_overall.csv
#   results/h0_diag/h0_mapping_diagnostic_summary.txt
#
# This script is deliberately self-contained and uses only base R.
# It does NOT estimate the TVP-GVAR and does NOT modify posterior draws.
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

results_dir <- Sys.getenv("TVPGVAR_RESULTS_DIR", "results")
out_dir <- Sys.getenv(
  "TVPGVAR_H0_DIAG_OUT",
  file.path(results_dir, "h0_diag")
)

pred_file <- Sys.getenv(
  "TVPGVAR_PRED_FILE",
  file.path(results_dir, "predDens_gpr_vix_structural.rda")
)

irf_file <- Sys.getenv(
  "TVPGVAR_IRF_FILE",
  file.path(results_dir, "irf_gpr_structural.rda")
)

selected_dates_env <- Sys.getenv(
  "TVPGVAR_DIAG_DATES",
  "2003Q1,2008Q3,2014Q3,2020Q1,2022Q1,2023Q4"
)

selected_dates <- trimws(
  strsplit(selected_dates_env, ",", fixed = TRUE)[[1]]
)
selected_dates <- selected_dates[nzchar(selected_dates)]

stable_only <- tolower(
  Sys.getenv("TVPGVAR_DIAG_STABLE_ONLY", "1")
) %in% c("1", "true", "yes", "y")

shock_pct <- as.numeric(
  Sys.getenv("TVPGVAR_GPR_SHOCK_PCT", "10")
)

driver_countries <- trimws(
  strsplit(
    Sys.getenv("TVPGVAR_DRIVER_COUNTRIES", "JP,UK"),
    ",",
    fixed = TRUE
  )[[1]]
)
driver_countries <- driver_countries[nzchar(driver_countries)]

if (!is.finite(shock_pct) || shock_pct <= 0) {
  stop("TVPGVAR_GPR_SHOCK_PCT must be positive.")
}

if (!file.exists(pred_file)) {
  stop("Missing posterior file: ", pred_file)
}

if (!file.exists(irf_file)) {
  stop("Missing IRF/stability file: ", irf_file)
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n============================================================\n")
cat(" LIGHTWEIGHT h=0 + JP/UK STABILITY DIAGNOSTIC\n")
cat(" Posterior: ", pred_file, "\n", sep = "")
cat(" IRF file:  ", irf_file, "\n", sep = "")
cat(" Stable-only selected draws: ", stable_only, "\n", sep = "")
cat(" Requested dates: ", paste(selected_dates, collapse = ", "), "\n", sep = "")
cat(" Shock normalization: +", shock_pct, "% US_gpr\n", sep = "")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

q <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, p, names = FALSE, type = 7))
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4L) return(NA_real_)
  sx <- sd(x[ok])
  sy <- sd(y[ok])
  if (!is.finite(sx) || !is.finite(sy) || sx == 0 || sy == 0) {
    return(NA_real_)
  }
  suppressWarnings(cor(x[ok], y[ok]))
}

stabilize_cov <- function(S, eps = 1e-10) {
  S <- as.matrix(S)
  S <- (S + t(S)) / 2
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  mn <- min(ev)
  if (!is.finite(mn)) {
    stop("Non-finite covariance eigenvalue.")
  }
  if (mn <= eps) {
    S <- S + diag(abs(mn) + eps, nrow(S))
  }
  S
}

spectral_radius <- function(F) {
  if (!length(F)) stop("No transition matrices supplied.")
  F <- lapply(F, as.matrix)

  K <- nrow(F[[1]])
  if (any(vapply(F, nrow, integer(1)) != K) ||
      any(vapply(F, ncol, integer(1)) != K)) {
    stop("Transition matrices must be square and conformable.")
  }

  p <- length(F)

  if (p == 1L) {
    vals <- eigen(F[[1]], only.values = TRUE)$values
    return(as.numeric(max(Mod(vals))))
  }

  top <- do.call(cbind, F)
  lower <- cbind(
    diag(K * (p - 1L)),
    matrix(0, nrow = K * (p - 1L), ncol = K)
  )
  companion <- rbind(top, lower)
  vals <- eigen(companion, only.values = TRUE)$values
  as.numeric(max(Mod(vals)))
}

condition_number <- function(G) {
  s <- svd(G, nu = 0, nv = 0)$d
  smax <- max(s)
  smin <- min(s)

  cond <- if (
    smin <= .Machine$double.eps * max(1, smax)
  ) {
    Inf
  } else {
    smax / smin
  }

  c(
    condition_number = as.numeric(cond),
    min_singular_value = as.numeric(smin),
    max_singular_value = as.numeric(smax)
  )
}

lag_order_from_alpha <- function(A) {
  rn <- dimnames(A)[[2]]
  if (is.null(rn)) {
    stop("ALPHA array has no regressor dimnames.")
  }

  tags <- unique(rn[grepl("^Ylag[0-9]+$", rn)])
  if (!length(tags)) {
    stop("No Ylag blocks found in ALPHA.")
  }

  nums <- as.integer(sub("^Ylag", "", tags))
  if (anyNA(nums)) {
    stop("Could not parse Ylag block labels.")
  }

  p <- max(nums)
  if (!p %in% c(1L, 2L)) {
    stop("This diagnostic currently expects p=1 or p=2.")
  }

  p
}

get_theta <- function(A, tt, draw, lag_no) {
  rn <- dimnames(A)[[2]]
  eq <- dimnames(A)[[3]]

  idx <- which(rn == paste0("Ylag", lag_no))
  if (length(idx) != length(eq)) {
    stop(
      "Unexpected domestic lag block width for ",
      paste(eq, collapse = ", "),
      " at Ylag", lag_no,
      ": found ", length(idx),
      ", expected ", length(eq)
    )
  }

  z <- A[tt, idx, , draw, drop = FALSE]
  z <- matrix(
    z,
    nrow = length(idx),
    ncol = length(eq)
  )

  theta <- t(z)
  rownames(theta) <- eq
  colnames(theta) <- eq
  theta
}

get_lambda0 <- function(A, tt, draw) {
  rn <- dimnames(A)[[2]]
  eq <- dimnames(A)[[3]]

  idx <- which(rn == "Wex")
  if (!length(idx)) {
    stop("No contemporaneous Wex block found.")
  }

  z <- A[tt, idx, , draw, drop = FALSE]
  z <- matrix(
    z,
    nrow = length(idx),
    ncol = length(eq)
  )

  lam <- t(z)
  rownames(lam) <- eq
  lam
}

build_global_G <- function(tt, draw, predDens, countries) {
  pieces <- vector("list", length(countries))

  for (i in seq_along(countries)) {
    A <- predDens[[i]]$ALPHA
    eq <- dimnames(A)[[3]]
    k_i <- length(eq)

    lam0 <- get_lambda0(A, tt, draw)
    Wi <- as.matrix(predDens[[i]]$W)

    A0 <- cbind(
      diag(k_i),
      -lam0
    )

    if (ncol(A0) != nrow(Wi)) {
      stop(
        "A0/W mismatch for ", countries[[i]],
        ": A0 cols=", ncol(A0),
        ", W rows=", nrow(Wi)
      )
    }

    pieces[[i]] <- A0 %*% Wi
  }

  G <- do.call(rbind, pieces)

  if (nrow(G) != ncol(G)) {
    stop("Reconstructed global G is not square.")
  }

  G
}

build_us_shock <- function(tt, draw, predDens, countries) {
  us_i <- match("US", countries)
  if (is.na(us_i)) {
    stop("US block not found.")
  }

  S <- predDens[[us_i]]$SIGMApost[tt, , , draw, drop = FALSE]
  k <- dim(predDens[[us_i]]$SIGMApost)[2]

  S <- matrix(S, nrow = k, ncol = k)
  S <- stabilize_cov(S)

  eq <- dimnames(predDens[[us_i]]$ALPHA)[[3]]
  if (length(eq) != k) {
    stop("US equation labels/covariance dimension mismatch.")
  }

  expected_us <- c(
    "US_gpr", "US_vix", "US_y", "US_dp",
    "US_r", "US_de", "US_deq"
  )

  if (!identical(eq, expected_us)) {
    stop(
      "Unexpected US recursive ordering. Found: ",
      paste(eq, collapse = ", ")
    )
  }

  L <- t(chol(S))
  shock_us <- L[, 1]
  names(shock_us) <- eq

  list(
    covariance = S,
    chol = L,
    shock_us = shock_us
  )
}

normalize_impact <- function(G, u_stack, shock_index, shock_pct) {
  raw <- as.numeric(solve(G, u_stack))

  if (!is.finite(raw[[shock_index]]) ||
      abs(raw[[shock_index]]) < 1e-12) {
    stop("Near-zero US_gpr impact; cannot normalize.")
  }

  target <- log1p(shock_pct / 100)
  scale <- target / raw[[shock_index]]

  list(
    raw = raw,
    normalized = raw * scale,
    scale = scale,
    target = target
  )
}

summarize_vector <- function(x) {
  c(
    median = safe_median(x),
    low68 = q(x, 0.16),
    high68 = q(x, 0.84),
    low90 = q(x, 0.05),
    high90 = q(x, 0.95)
  )
}

# -----------------------------------------------------------------------------
# Load completed posterior and saved IRF/stability objects
# -----------------------------------------------------------------------------

load(pred_file)

if (!exists("predDens") || !exists("Data.setup")) {
  stop(
    pred_file,
    " must contain predDens and Data.setup."
  )
}

load(irf_file)

if (!exists("IRF_post") ||
    !exists("stability_rho") ||
    !exists("stable_mask")) {
  stop(
    irf_file,
    " must contain IRF_post, stability_rho and stable_mask."
  )
}

countries <- Data.setup$countries
xglobal <- Data.setup$bigx
var_names <- colnames(xglobal)

if (is.null(var_names)) {
  stop("Data.setup$bigx has no column names.")
}

if (!all(driver_countries %in% countries)) {
  stop(
    "Requested driver countries not found: ",
    paste(
      setdiff(driver_countries, countries),
      collapse = ", "
    )
  )
}

A_list <- lapply(predDens, `[[`, "ALPHA")

p_each <- vapply(A_list, lag_order_from_alpha, integer(1))
if (length(unique(p_each)) != 1L) {
  stop("Countries do not share one domestic lag order.")
}
p <- p_each[[1]]

n_draws <- dim(A_list[[1]])[4]
if (any(vapply(A_list, function(A) dim(A)[4], integer(1)) != n_draws)) {
  stop("Countries have inconsistent retained draw counts.")
}

irf_dates <- dimnames(IRF_post)[[1]]
irf_vars <- dimnames(IRF_post)[[2]]

if (is.null(irf_dates) || is.null(irf_vars)) {
  stop("IRF_post is missing date/variable dimnames.")
}

if (!identical(irf_vars, var_names)) {
  stop("IRF variable order does not match Data.setup$bigx.")
}

if (!identical(rownames(stable_mask), irf_dates)) {
  stop("stable_mask date order does not match IRF_post.")
}

available_dates <- intersect(selected_dates, irf_dates)
missing_dates <- setdiff(selected_dates, irf_dates)

if (length(missing_dates)) {
  warning(
    "Requested dates not available and will be skipped: ",
    paste(missing_dates, collapse = ", ")
  )
}

if (!length(available_dates)) {
  stop("None of the requested diagnostic dates are available.")
}

gdp_vars <- paste0(countries, "_y")
if (!all(gdp_vars %in% var_names)) {
  stop(
    "GDP variables missing from global vector: ",
    paste(setdiff(gdp_vars, var_names), collapse = ", ")
  )
}

shock_index <- match("US_gpr", var_names)
if (is.na(shock_index)) {
  stop("US_gpr not found in global variable vector.")
}

cat(
  "Loaded posterior: p=", p,
  ", draws=", n_draws,
  ", countries=", length(countries),
  ", global variables=", length(var_names),
  "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# A. Reconstruct h=0 mapping and compare full recursive vs GPR-only shock
# -----------------------------------------------------------------------------

h0_draw_rows <- list()
loading_rows <- list()
global_draw_rows <- list()
condition_rows <- list()

counter <- 0L
max_reproduction_error <- 0

for (date in available_dates) {
  tt <- match(date, irf_dates)

  draws <- seq_len(n_draws)
  if (stable_only) {
    draws <- draws[stable_mask[tt, draws]]
  }

  if (!length(draws)) {
    warning("No usable draws for ", date)
    next
  }

  cat(
    "h=0 mapping: ", date,
    " | draws used=", length(draws),
    "\n",
    sep = ""
  )

  for (draw in draws) {
    G <- build_global_G(
      tt = tt,
      draw = draw,
      predDens = predDens,
      countries = countries
    )

    shock <- build_us_shock(
      tt = tt,
      draw = draw,
      predDens = predDens,
      countries = countries
    )

    # Full recursive shock: complete first Cholesky column in the US block.
    blocks_full <- lapply(
      seq_along(countries),
      function(i) {
        cc <- countries[[i]]
        k_i <- dim(predDens[[i]]$ALPHA)[3]

        if (cc == "US") {
          as.numeric(shock$shock_us)
        } else {
          rep(0, k_i)
        }
      }
    )

    u_full <- unlist(blocks_full, use.names = FALSE)

    if (length(u_full) != length(var_names)) {
      stop("Stacked full innovation dimension mismatch.")
    }

    full <- normalize_impact(
      G = G,
      u_stack = u_full,
      shock_index = shock_index,
      shock_pct = shock_pct
    )

    names(full$normalized) <- var_names

    # GPR-only counterfactual:
    # keep only L[US_gpr, first shock] and zero all later US loadings.
    us_only <- rep(0, length(shock$shock_us))
    us_only[[1]] <- shock$shock_us[[1]]

    blocks_only <- lapply(
      seq_along(countries),
      function(i) {
        cc <- countries[[i]]
        k_i <- dim(predDens[[i]]$ALPHA)[3]

        if (cc == "US") {
          us_only
        } else {
          rep(0, k_i)
        }
      }
    )

    u_only <- unlist(blocks_only, use.names = FALSE)

    only <- normalize_impact(
      G = G,
      u_stack = u_only,
      shock_index = shock_index,
      shock_pct = shock_pct
    )

    names(only$normalized) <- var_names

    # Exact reconstruction check against saved h=0 IRF.
    saved_h0 <- IRF_post[tt, , 1, draw]
    err <- max(
      abs(full$normalized - saved_h0),
      na.rm = TRUE
    )
    max_reproduction_error <- max(
      max_reproduction_error,
      err
    )

    if (!is.finite(err) || err > 1e-7) {
      stop(
        "Reconstructed h=0 IRF does not match saved IRF at ",
        date, ", draw ", draw,
        ". max abs diff=", signif(err, 8)
      )
    }

    cond <- condition_number(G)
    condition_rows[[length(condition_rows) + 1L]] <- data.frame(
      date = date,
      draw = draw,
      stable = stable_mask[tt, draw],
      global_rho = stability_rho[tt, draw],
      G_condition_number = unname(cond[["condition_number"]]),
      G_min_singular_value = unname(cond[["min_singular_value"]]),
      G_max_singular_value = unname(cond[["max_singular_value"]]),
      stringsAsFactors = FALSE
    )

    # US recursive first-column loadings.
    for (j in seq_along(shock$shock_us)) {
      loading_rows[[length(loading_rows) + 1L]] <- data.frame(
        date = date,
        draw = draw,
        component = names(shock$shock_us)[[j]],
        raw_cholesky_loading = shock$shock_us[[j]],
        full_normalization_scale = full$scale,
        normalized_residual_loading = shock$shock_us[[j]] * full$scale,
        stringsAsFactors = FALSE
      )
    }

    # All global h=0 variables, summarized later.
    for (j in seq_along(var_names)) {
      global_draw_rows[[length(global_draw_rows) + 1L]] <- data.frame(
        date = date,
        draw = draw,
        variable = var_names[[j]],
        full_h0 = full$normalized[[j]],
        gpr_only_h0 = only$normalized[[j]],
        difference_full_minus_gpr_only =
          full$normalized[[j]] - only$normalized[[j]],
        stringsAsFactors = FALSE
      )
    }

    # GDP-only draw-level table.
    for (cc in countries) {
      vv <- paste0(cc, "_y")
      h0_draw_rows[[length(h0_draw_rows) + 1L]] <- data.frame(
        date = date,
        draw = draw,
        country = cc,
        variable = vv,
        full_recursive_h0 = full$normalized[[vv]],
        gpr_only_h0 = only$normalized[[vv]],
        difference_full_minus_gpr_only =
          full$normalized[[vv]] - only$normalized[[vv]],
        full_positive = full$normalized[[vv]] > 0,
        gpr_only_positive = only$normalized[[vv]] > 0,
        stringsAsFactors = FALSE
      )
    }

    counter <- counter + 1L
  }
}

if (!length(h0_draw_rows)) {
  stop("No h=0 diagnostic rows were produced.")
}

h0_draws <- do.call(rbind, h0_draw_rows)
loading_draws <- do.call(rbind, loading_rows)
global_draws <- do.call(rbind, global_draw_rows)
condition_draws <- do.call(rbind, condition_rows)

write.csv(
  h0_draws,
  file.path(out_dir, "h0_full_vs_gpr_only_draws.csv"),
  row.names = FALSE
)

write.csv(
  loading_draws,
  file.path(out_dir, "us_recursive_gpr_loading_draws.csv"),
  row.names = FALSE
)

write.csv(
  condition_draws,
  file.path(out_dir, "selected_date_G_condition_draws.csv"),
  row.names = FALSE
)

# GDP summary by date-country.
gdp_split <- split(
  h0_draws,
  interaction(
    h0_draws$date,
    h0_draws$country,
    drop = TRUE,
    lex.order = TRUE
  )
)

h0_gdp_summary <- do.call(
  rbind,
  lapply(gdp_split, function(z) {
    sf <- summarize_vector(z$full_recursive_h0)
    so <- summarize_vector(z$gpr_only_h0)
    sd <- summarize_vector(z$difference_full_minus_gpr_only)

    data.frame(
      date = z$date[[1]],
      country = z$country[[1]],
      draws = nrow(z),

      full_median = sf[["median"]],
      full_low68 = sf[["low68"]],
      full_high68 = sf[["high68"]],
      full_low90 = sf[["low90"]],
      full_high90 = sf[["high90"]],
      full_positive_share = mean(z$full_positive),

      gpr_only_median = so[["median"]],
      gpr_only_low68 = so[["low68"]],
      gpr_only_high68 = so[["high68"]],
      gpr_only_low90 = so[["low90"]],
      gpr_only_high90 = so[["high90"]],
      gpr_only_positive_share = mean(z$gpr_only_positive),

      difference_median = sd[["median"]],
      difference_low68 = sd[["low68"]],
      difference_high68 = sd[["high68"]],

      median_sign_changes =
        sign(sf[["median"]]) != sign(so[["median"]]),

      full_68_excludes_zero =
        is.finite(sf[["low68"]]) &&
        is.finite(sf[["high68"]]) &&
        (sf[["low68"]] > 0 || sf[["high68"]] < 0),

      gpr_only_68_excludes_zero =
        is.finite(so[["low68"]]) &&
        is.finite(so[["high68"]]) &&
        (so[["low68"]] > 0 || so[["high68"]] < 0),

      stringsAsFactors = FALSE
    )
  })
)

h0_gdp_summary <- h0_gdp_summary[
  order(
    match(h0_gdp_summary$date, available_dates),
    match(h0_gdp_summary$country, countries)
  ),
  ,
  drop = FALSE
]

write.csv(
  h0_gdp_summary,
  file.path(out_dir, "h0_gdp_full_vs_gpr_only_summary.csv"),
  row.names = FALSE
)

# All-variable global h=0 summary.
global_split <- split(
  global_draws,
  interaction(
    global_draws$date,
    global_draws$variable,
    drop = TRUE,
    lex.order = TRUE
  )
)

h0_global_summary <- do.call(
  rbind,
  lapply(global_split, function(z) {
    sf <- summarize_vector(z$full_h0)
    so <- summarize_vector(z$gpr_only_h0)
    sd <- summarize_vector(z$difference_full_minus_gpr_only)

    data.frame(
      date = z$date[[1]],
      variable = z$variable[[1]],
      draws = nrow(z),

      full_median = sf[["median"]],
      full_low68 = sf[["low68"]],
      full_high68 = sf[["high68"]],

      gpr_only_median = so[["median"]],
      gpr_only_low68 = so[["low68"]],
      gpr_only_high68 = so[["high68"]],

      difference_median = sd[["median"]],
      difference_low68 = sd[["low68"]],
      difference_high68 = sd[["high68"]],

      stringsAsFactors = FALSE
    )
  })
)

h0_global_summary <- h0_global_summary[
  order(
    match(h0_global_summary$date, available_dates),
    match(h0_global_summary$variable, var_names)
  ),
  ,
  drop = FALSE
]

write.csv(
  h0_global_summary,
  file.path(out_dir, "h0_global_impact_summary.csv"),
  row.names = FALSE
)

# US first Cholesky-column summary.
loading_split <- split(
  loading_draws,
  interaction(
    loading_draws$date,
    loading_draws$component,
    drop = TRUE,
    lex.order = TRUE
  )
)

loading_summary <- do.call(
  rbind,
  lapply(loading_split, function(z) {
    sr <- summarize_vector(z$raw_cholesky_loading)
    sn <- summarize_vector(z$normalized_residual_loading)

    data.frame(
      date = z$date[[1]],
      component = z$component[[1]],
      draws = nrow(z),

      raw_loading_median = sr[["median"]],
      raw_loading_low68 = sr[["low68"]],
      raw_loading_high68 = sr[["high68"]],

      normalized_loading_median = sn[["median"]],
      normalized_loading_low68 = sn[["low68"]],
      normalized_loading_high68 = sn[["high68"]],

      normalized_positive_share =
        mean(z$normalized_residual_loading > 0),

      stringsAsFactors = FALSE
    )
  })
)

loading_summary <- loading_summary[
  order(
    match(loading_summary$date, available_dates),
    match(
      loading_summary$component,
      c(
        "US_gpr", "US_vix", "US_y", "US_dp",
        "US_r", "US_de", "US_deq"
      )
    )
  ),
  ,
  drop = FALSE
]

write.csv(
  loading_summary,
  file.path(out_dir, "us_recursive_gpr_loading_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# B. JP / UK local stability-driver diagnostic
# -----------------------------------------------------------------------------

driver_rows <- list()

for (cc in driver_countries) {
  i <- match(cc, countries)
  A <- predDens[[i]]$ALPHA
  eq_names <- dimnames(A)[[3]]

  if (is.null(eq_names)) {
    stop("Missing equation labels for ", cc)
  }

  cat("Local stability drivers: ", cc, "\n", sep = "")

  for (date in available_dates) {
    tt <- match(date, irf_dates)

    # For driver diagnosis we use all posterior draws by default.
    # This is intentional: unstable draws are precisely what we want to explain.
    draws <- seq_len(n_draws)

    for (draw in draws) {
      theta <- lapply(
        seq_len(p),
        function(lag_no) {
          get_theta(
            A = A,
            tt = tt,
            draw = draw,
            lag_no = lag_no
          )
        }
      )

      rho <- spectral_radius(theta)

      for (e in seq_along(eq_names)) {
        cf <- lapply(theta, function(M) {
          M2 <- M
          M2[e, ] <- 0
          M2
        })

        rho_cf <- spectral_radius(cf)

        ar1 <- theta[[1]][e, e]
        ar2 <- if (p >= 2L) theta[[2]][e, e] else 0
        arsum <- ar1 + ar2

        driver_rows[[length(driver_rows) + 1L]] <- data.frame(
          date = date,
          draw = draw,
          country = cc,
          equation = eq_names[[e]],

          local_rho = rho,
          local_stable = rho < 1,

          own_ar1 = ar1,
          own_ar2 = ar2,
          own_ar_sum = arsum,

          rho_if_equation_row_zero = rho_cf,
          rho_reduction_if_row_zero = rho - rho_cf,

          row_zero_makes_stable =
            (rho >= 1 && rho_cf < 1),

          stringsAsFactors = FALSE
        )
      }
    }
  }
}

driver_draws <- do.call(rbind, driver_rows)

write.csv(
  driver_draws,
  file.path(out_dir, "jp_uk_stability_driver_draws.csv"),
  row.names = FALSE
)

driver_split <- split(
  driver_draws,
  interaction(
    driver_draws$date,
    driver_draws$country,
    driver_draws$equation,
    drop = TRUE,
    lex.order = TRUE
  )
)

driver_summary <- do.call(
  rbind,
  lapply(driver_split, function(z) {
    unstable <- z$local_rho >= 1
    rescue <- z$row_zero_makes_stable[unstable]

    data.frame(
      date = z$date[[1]],
      country = z$country[[1]],
      equation = z$equation[[1]],
      draws = nrow(z),

      local_stable_share = mean(z$local_stable),
      median_local_rho = median(z$local_rho),
      p95_local_rho = q(z$local_rho, 0.95),

      median_own_ar1 = median(z$own_ar1),
      median_own_ar2 = median(z$own_ar2),
      median_own_ar_sum = median(z$own_ar_sum),

      corr_own_ar_sum_with_rho =
        safe_cor(z$own_ar_sum, z$local_rho),

      median_rho_reduction_if_row_zero =
        median(z$rho_reduction_if_row_zero),

      p90_rho_reduction_if_row_zero =
        q(z$rho_reduction_if_row_zero, 0.90),

      unstable_draws = sum(unstable),

      rescue_share_among_unstable =
        if (sum(unstable) > 0L) {
          mean(rescue)
        } else {
          NA_real_
        },

      stringsAsFactors = FALSE
    )
  })
)

driver_summary <- driver_summary[
  order(
    match(driver_summary$date, available_dates),
    match(driver_summary$country, driver_countries),
    driver_summary$equation
  ),
  ,
  drop = FALSE
]

write.csv(
  driver_summary,
  file.path(out_dir, "jp_uk_stability_driver_summary.csv"),
  row.names = FALSE
)

overall_split <- split(
  driver_draws,
  interaction(
    driver_draws$country,
    driver_draws$equation,
    drop = TRUE,
    lex.order = TRUE
  )
)

driver_overall <- do.call(
  rbind,
  lapply(overall_split, function(z) {
    unstable <- z$local_rho >= 1

    data.frame(
      country = z$country[[1]],
      equation = z$equation[[1]],
      observations = nrow(z),

      stable_share = mean(z$local_stable),

      median_own_ar1 = median(z$own_ar1),
      median_own_ar2 = median(z$own_ar2),
      median_own_ar_sum = median(z$own_ar_sum),

      corr_own_ar_sum_with_rho =
        safe_cor(z$own_ar_sum, z$local_rho),

      median_rho_reduction_if_row_zero =
        median(z$rho_reduction_if_row_zero),

      p90_rho_reduction_if_row_zero =
        q(z$rho_reduction_if_row_zero, 0.90),

      rescue_share_among_unstable =
        if (sum(unstable) > 0L) {
          mean(z$row_zero_makes_stable[unstable])
        } else {
          NA_real_
        },

      stringsAsFactors = FALSE
    )
  })
)

driver_overall <- driver_overall[
  order(
    match(driver_overall$country, driver_countries),
    -driver_overall$median_rho_reduction_if_row_zero
  ),
  ,
  drop = FALSE
]

write.csv(
  driver_overall,
  file.path(out_dir, "jp_uk_stability_driver_overall.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# C. Human-readable summary
# -----------------------------------------------------------------------------

summary_lines <- c(
  "P2-Strong GPR + VIX: lightweight h=0 / JP-UK diagnostic",
  paste0("Posterior file: ", pred_file),
  paste0("IRF file: ", irf_file),
  paste0("Domestic lag order p: ", p),
  paste0("Retained posterior draws: ", n_draws),
  paste0("Selected dates: ", paste(available_dates, collapse = ", ")),
  paste0("Stable-only h=0 analysis: ", stable_only),
  paste0("GPR normalization: +", shock_pct, "%"),
  paste0(
    "Maximum absolute error when reproducing saved h=0 IRF: ",
    format(max_reproduction_error, scientific = TRUE, digits = 6)
  ),
  "",
  "INTERPRETATION OF THE COUNTERFACTUAL",
  "Full recursive = complete first Cholesky column of the US covariance block.",
  "GPR-only = retain only the US_gpr element of that column and zero",
  "the immediate US_vix/US_y/US_dp/US_r/US_de/US_deq covariance loadings.",
  "Both are separately normalized to the same +GPR impact.",
  "",
  "If full GDP is positive but GPR-only GDP becomes zero/negative:",
  "the positive h=0 GDP response is mainly generated by the recursive",
  "Cholesky covariance loadings in the US block.",
  "",
  "If GPR-only GDP remains systematically positive:",
  "the positive h=0 response is already being generated by the",
  "contemporaneous global G^{-1} mapping / W matrices.",
  ""
)

# Date-level count of positive GDP medians in both specifications.
date_count_lines <- unlist(
  lapply(available_dates, function(date) {
    z <- h0_gdp_summary[h0_gdp_summary$date == date, , drop = FALSE]

    paste0(
      date,
      ": positive median GDP | full=",
      sum(z$full_median > 0, na.rm = TRUE),
      "/",
      nrow(z),
      ", GPR-only=",
      sum(z$gpr_only_median > 0, na.rm = TRUE),
      "/",
      nrow(z),
      ", median-sign changes=",
      sum(z$median_sign_changes, na.rm = TRUE),
      "/",
      nrow(z)
    )
  })
)

summary_lines <- c(
  summary_lines,
  "h=0 GDP SIGN COUNTS",
  date_count_lines,
  ""
)

# Top JP/UK equation driver by median rho reduction.
top_driver_lines <- unlist(
  lapply(driver_countries, function(cc) {
    z <- driver_overall[
      driver_overall$country == cc,
      ,
      drop = FALSE
    ]

    if (!nrow(z)) {
      return(paste0(cc, ": no driver rows"))
    }

    z <- z[
      order(
        -z$median_rho_reduction_if_row_zero
      ),
      ,
      drop = FALSE
    ]

    top <- z[1, ]

    paste0(
      cc,
      ": largest median rho reduction from zeroing equation row = ",
      top$equation,
      " (median reduction=",
      signif(top$median_rho_reduction_if_row_zero, 5),
      ", rescue share among unstable=",
      ifelse(
        is.finite(top$rescue_share_among_unstable),
        paste0(
          round(
            100 * top$rescue_share_among_unstable,
            1
          ),
          "%"
        ),
        "NA"
      ),
      ")"
    )
  })
)

summary_lines <- c(
  summary_lines,
  "JP/UK LOCAL-STABILITY DRIVER CHECK",
  top_driver_lines,
  "",
  "Important: row-zero results are diagnostic counterfactuals only.",
  "They do not alter estimation and should not be reported as causal contributions."
)

writeLines(
  summary_lines,
  file.path(out_dir, "h0_mapping_diagnostic_summary.txt")
)

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n\nCreated diagnostic files in: ", out_dir, "\n", sep = "")
