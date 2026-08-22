#!/usr/bin/env Rscript

# =============================================================================
# Post-estimation contemporaneous-channel diagnostic for P2-Strong GPR + VIX
#
# DOES NOT re-run MCMC.
#
# Uses the completed posterior from the GPR+VIX run and decomposes the h=0 GDP
# response under a pure GPR residual innovation.
#
# Scenarios:
#   A_baseline
#       Original contemporaneous G matrix.
#
#   B1_zero_nonUS_GPR_GDP
#       Set the contemporaneous global-GPR coefficient to zero ONLY in each
#       non-US GDP equation.
#
#   B2_zero_nonUS_GPR_all
#       Set the contemporaneous global-GPR coefficient to zero in ALL non-US
#       equations.
#
#   D_zero_foreign_macro
#       Keep non-US global-GPR/global-VIX coefficients, but set all five
#       contemporaneous foreign-macro coefficients to zero in every country.
#
#   C_zero_GPR_all_plus_foreign_macro
#       Combine B2 + zero all contemporaneous foreign-macro coefficients.
#       Non-US global-VIX coefficients are left unchanged.
#
# All scenarios use the SAME pure GPR residual innovation:
#   - retain only the first element of the first US Cholesky column (US_gpr)
#   - set immediate US VIX/GDP/inflation/rate/FX/equity residual loadings to 0
#   - normalize each scenario to the same +10% US_gpr impact
#
# This isolates where the positive contemporaneous GDP response is generated:
#   non-US global-GPR coefficients vs contemporaneous foreign-macro feedback.
#
# Expected model ordering confirmed by the current repository:
#   US endogenous: GPR, VIX, y, dp, r, de, deq
#   non-US Wex: foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq,
#               global_gpr, global_vix
#
# INPUTS
#   results/predDens_gpr_vix_structural.rda
#   results/irf_gpr_structural.rda
#
# OUTPUTS
#   results/channel_diag/h0_channel_scenario_draws.csv
#   results/channel_diag/h0_channel_scenario_summary.csv
#   results/channel_diag/h0_sign_counts_by_scenario.csv
#   results/channel_diag/nonUS_gdp_wex_coeff_draws.csv
#   results/channel_diag/nonUS_gdp_wex_coeff_summary.csv
#   results/channel_diag/nonUS_global_gpr_coeff_summary_all_equations.csv
#   results/channel_diag/scenario_G_condition_summary.csv
#   results/channel_diag/h0_channel_decomposition_summary.txt
# =============================================================================

options(stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

results_dir <- Sys.getenv("TVPGVAR_RESULTS_DIR", "results")
out_dir <- Sys.getenv(
  "TVPGVAR_CHANNEL_DIAG_OUT",
  file.path(results_dir, "channel_diag")
)

pred_file <- Sys.getenv(
  "TVPGVAR_PRED_FILE",
  file.path(results_dir, "predDens_gpr_vix_structural.rda")
)

irf_file <- Sys.getenv(
  "TVPGVAR_IRF_FILE",
  file.path(results_dir, "irf_gpr_structural.rda")
)

selected_dates <- trimws(
  strsplit(
    Sys.getenv(
      "TVPGVAR_DIAG_DATES",
      "2003Q1,2008Q3,2014Q3,2020Q1,2022Q1,2023Q4"
    ),
    ",",
    fixed = TRUE
  )[[1]]
)
selected_dates <- selected_dates[nzchar(selected_dates)]

stable_only <- tolower(
  Sys.getenv("TVPGVAR_DIAG_STABLE_ONLY", "1")
) %in% c("1", "true", "yes", "y")

shock_pct <- as.numeric(
  Sys.getenv("TVPGVAR_GPR_SHOCK_PCT", "10")
)

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
cat(" CONTEMPORANEOUS GPR CHANNEL DECOMPOSITION\n")
cat(" Posterior: ", pred_file, "\n", sep = "")
cat(" IRF file:  ", irf_file, "\n", sep = "")
cat(" Stable-only draw selection: ", stable_only, "\n", sep = "")
cat(" Shock normalization: +", shock_pct, "% US_gpr\n", sep = "")
cat(" Dates: ", paste(selected_dates, collapse = ", "), "\n", sep = "")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------

qv <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, p, names = FALSE, type = 7))
}

safe_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else median(x)
}

summarize_vec <- function(x) {
  c(
    median = safe_median(x),
    low68 = qv(x, 0.16),
    high68 = qv(x, 0.84),
    low90 = qv(x, 0.05),
    high90 = qv(x, 0.95)
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

get_lambda0 <- function(A, tt, draw) {
  rn <- dimnames(A)[[2]]
  eq <- dimnames(A)[[3]]
  if (is.null(rn) || is.null(eq)) {
    stop("ALPHA array is missing regressor/equation names.")
  }

  idx <- which(rn == "Wex")
  if (!length(idx)) stop("No Wex block found in ALPHA.")

  z <- A[tt, idx, , draw, drop = FALSE]
  z <- matrix(z, nrow = length(idx), ncol = length(eq))
  lam <- t(z)
  rownames(lam) <- eq
  lam
}

wex_labels_for_country <- function(cc) {
  if (identical(cc, "US")) {
    c(
      "foreign_y", "foreign_dp", "foreign_r",
      "foreign_de", "foreign_deq"
    )
  } else {
    c(
      "foreign_y", "foreign_dp", "foreign_r",
      "foreign_de", "foreign_deq",
      "global_gpr", "global_vix"
    )
  }
}

mutate_lambda0 <- function(lam, cc, scenario) {
  labels <- wex_labels_for_country(cc)
  if (ncol(lam) != length(labels)) {
    stop(
      "Unexpected Wex width for ", cc,
      ": found ", ncol(lam),
      ", expected ", length(labels)
    )
  }
  colnames(lam) <- labels

  if (scenario == "A_baseline") {
    return(lam)
  }

  if (scenario == "B1_zero_nonUS_GPR_GDP") {
    if (!identical(cc, "US")) {
      yrow <- match(paste0(cc, "_y"), rownames(lam))
      if (is.na(yrow)) stop("GDP equation not found for ", cc)
      lam[yrow, "global_gpr"] <- 0
    }
    return(lam)
  }

  if (scenario == "B2_zero_nonUS_GPR_all") {
    if (!identical(cc, "US")) {
      lam[, "global_gpr"] <- 0
    }
    return(lam)
  }

  if (scenario == "D_zero_foreign_macro") {
    lam[, c(
      "foreign_y", "foreign_dp", "foreign_r",
      "foreign_de", "foreign_deq"
    )] <- 0
    return(lam)
  }

  if (scenario == "C_zero_GPR_all_plus_foreign_macro") {
    lam[, c(
      "foreign_y", "foreign_dp", "foreign_r",
      "foreign_de", "foreign_deq"
    )] <- 0
    if (!identical(cc, "US")) {
      lam[, "global_gpr"] <- 0
    }
    return(lam)
  }

  stop("Unknown scenario: ", scenario)
}

build_global_G <- function(tt, draw, predDens, countries, scenario) {
  blocks <- vector("list", length(countries))

  for (i in seq_along(countries)) {
    cc <- countries[[i]]
    A <- predDens[[i]]$ALPHA
    eq <- dimnames(A)[[3]]
    k_i <- length(eq)

    lam0 <- get_lambda0(A, tt, draw)
    lam0 <- mutate_lambda0(lam0, cc, scenario)

    Wi <- as.matrix(predDens[[i]]$W)
    A0 <- cbind(diag(k_i), -lam0)

    if (ncol(A0) != nrow(Wi)) {
      stop(
        "A0/W mismatch for ", cc,
        " under ", scenario,
        ": A0 cols=", ncol(A0),
        ", W rows=", nrow(Wi)
      )
    }

    blocks[[i]] <- A0 %*% Wi
  }

  G <- do.call(rbind, blocks)
  if (nrow(G) != ncol(G)) {
    stop("Global G is not square under ", scenario)
  }
  G
}

build_us_shock_vectors <- function(tt, draw, predDens, countries) {
  us_i <- match("US", countries)
  if (is.na(us_i)) stop("US block not found.")

  Sarr <- predDens[[us_i]]$SIGMApost
  k <- dim(Sarr)[2]
  S <- Sarr[tt, , , draw, drop = FALSE]
  S <- matrix(S, nrow = k, ncol = k)
  S <- stabilize_cov(S)

  eq <- dimnames(predDens[[us_i]]$ALPHA)[[3]]
  expected <- c(
    "US_gpr", "US_vix", "US_y", "US_dp",
    "US_r", "US_de", "US_deq"
  )
  if (!identical(eq, expected)) {
    stop(
      "Unexpected US ordering. Found: ",
      paste(eq, collapse = ", ")
    )
  }

  L <- t(chol(S))
  full_us <- L[, 1]
  names(full_us) <- eq

  pure_us <- rep(0, length(full_us))
  pure_us[1] <- full_us[1]
  names(pure_us) <- eq

  stack_blocks <- function(us_vec) {
    unlist(
      lapply(seq_along(countries), function(i) {
        cc <- countries[[i]]
        k_i <- dim(predDens[[i]]$ALPHA)[3]
        if (identical(cc, "US")) us_vec else rep(0, k_i)
      }),
      use.names = FALSE
    )
  }

  list(
    full = stack_blocks(full_us),
    pure = stack_blocks(pure_us),
    full_us = full_us,
    pure_us = pure_us
  )
}

solve_normalized <- function(G, u, shock_index, shock_pct) {
  ans <- tryCatch(
    as.numeric(solve(G, u)),
    error = function(e) NULL
  )

  if (is.null(ans) || any(!is.finite(ans))) {
    return(list(success = FALSE, impact = rep(NA_real_, nrow(G))))
  }

  if (!is.finite(ans[shock_index]) || abs(ans[shock_index]) < 1e-12) {
    return(list(success = FALSE, impact = rep(NA_real_, nrow(G))))
  }

  target <- log1p(shock_pct / 100)
  ans <- ans * (target / ans[shock_index])

  list(success = TRUE, impact = ans)
}

G_condition <- function(G) {
  sv <- svd(G, nu = 0, nv = 0)$d
  smax <- max(sv)
  smin <- min(sv)
  cond <- if (
    smin <= .Machine$double.eps * max(1, smax)
  ) Inf else smax / smin

  c(
    condition_number = as.numeric(cond),
    min_singular_value = as.numeric(smin),
    max_singular_value = as.numeric(smax)
  )
}

# -----------------------------------------------------------------------------
# Load already-estimated objects
# -----------------------------------------------------------------------------

load(pred_file)
if (!exists("predDens") || !exists("Data.setup")) {
  stop(pred_file, " must contain predDens and Data.setup.")
}

load(irf_file)
if (!exists("IRF_post") || !exists("stability_rho") || !exists("stable_mask")) {
  stop(irf_file, " must contain IRF_post, stability_rho and stable_mask.")
}

countries <- Data.setup$countries
xglobal <- Data.setup$bigx
var_names <- colnames(xglobal)
irf_dates <- dimnames(IRF_post)[[1]]

if (is.null(var_names) || is.null(irf_dates)) {
  stop("Missing global variable names or IRF dates.")
}

if (!identical(dimnames(IRF_post)[[2]], var_names)) {
  stop("IRF variable order differs from Data.setup$bigx.")
}

available_dates <- intersect(selected_dates, irf_dates)
missing_dates <- setdiff(selected_dates, irf_dates)
if (length(missing_dates)) {
  warning(
    "Requested dates not available and skipped: ",
    paste(missing_dates, collapse = ", ")
  )
}
if (!length(available_dates)) stop("No requested dates are available.")

n_draws <- dim(predDens[[1]]$ALPHA)[4]
if (any(vapply(predDens, function(z) dim(z$ALPHA)[4], integer(1)) != n_draws)) {
  stop("Countries have inconsistent retained draw counts.")
}

shock_index <- match("US_gpr", var_names)
if (is.na(shock_index)) stop("US_gpr not found in global vector.")

gdp_vars <- paste0(countries, "_y")
if (!all(gdp_vars %in% var_names)) {
  stop("One or more GDP variables are missing from global vector.")
}

scenarios <- c(
  "A_baseline",
  "B1_zero_nonUS_GPR_GDP",
  "B2_zero_nonUS_GPR_all",
  "D_zero_foreign_macro",
  "C_zero_GPR_all_plus_foreign_macro"
)

cat(
  "Loaded posterior: countries=", length(countries),
  ", draws=", n_draws,
  ", global variables=", length(var_names),
  "\n",
  sep = ""
)

# -----------------------------------------------------------------------------
# Main loop
# -----------------------------------------------------------------------------

impact_rows <- list()
cond_rows <- list()
gdp_coeff_rows <- list()
gpr_all_eq_rows <- list()

max_baseline_reproduction_error <- 0

for (date in available_dates) {
  tt <- match(date, irf_dates)
  draws <- seq_len(n_draws)
  if (stable_only) draws <- draws[stable_mask[tt, draws]]

  if (!length(draws)) {
    warning("No usable draws for ", date)
    next
  }

  cat(date, ": draws used=", length(draws), "\n", sep = "")

  for (draw in draws) {
    shocks <- build_us_shock_vectors(tt, draw, predDens, countries)

    # Baseline full recursive reconstruction check.
    G_base <- build_global_G(
      tt, draw, predDens, countries,
      scenario = "A_baseline"
    )
    base_full <- solve_normalized(
      G_base, shocks$full, shock_index, shock_pct
    )
    if (!base_full$success) {
      stop("Baseline G solve failed at ", date, ", draw ", draw)
    }
    saved_h0 <- IRF_post[tt, , 1, draw]
    err <- max(abs(base_full$impact - saved_h0), na.rm = TRUE)
    max_baseline_reproduction_error <- max(
      max_baseline_reproduction_error,
      err
    )
    if (!is.finite(err) || err > 1e-7) {
      stop(
        "Baseline reconstruction mismatch at ", date,
        ", draw ", draw,
        ". max abs diff=", signif(err, 8)
      )
    }

    # Record baseline contemporaneous coefficients.
    for (i in seq_along(countries)) {
      cc <- countries[[i]]
      lam <- get_lambda0(predDens[[i]]$ALPHA, tt, draw)
      labels <- wex_labels_for_country(cc)
      if (ncol(lam) != length(labels)) {
        stop("Unexpected Wex width while extracting coefficients for ", cc)
      }
      colnames(lam) <- labels

      if (!identical(cc, "US")) {
        yrow <- match(paste0(cc, "_y"), rownames(lam))
        if (is.na(yrow)) stop("GDP row not found for ", cc)

        for (term in labels) {
          gdp_coeff_rows[[length(gdp_coeff_rows) + 1L]] <- data.frame(
            date = date,
            draw = draw,
            country = cc,
            equation = paste0(cc, "_y"),
            term = term,
            coefficient = lam[yrow, term],
            stringsAsFactors = FALSE
          )
        }

        for (eq in rownames(lam)) {
          gpr_all_eq_rows[[length(gpr_all_eq_rows) + 1L]] <- data.frame(
            date = date,
            draw = draw,
            country = cc,
            equation = eq,
            global_gpr_coefficient = lam[eq, "global_gpr"],
            stringsAsFactors = FALSE
          )
        }
      }
    }

    # Scenario comparison uses PURE GPR residual shock only.
    for (scenario in scenarios) {
      Gs <- if (scenario == "A_baseline") {
        G_base
      } else {
        build_global_G(tt, draw, predDens, countries, scenario)
      }

      cond <- G_condition(Gs)
      ans <- solve_normalized(
        Gs, shocks$pure, shock_index, shock_pct
      )

      cond_rows[[length(cond_rows) + 1L]] <- data.frame(
        date = date,
        draw = draw,
        scenario = scenario,
        solve_success = ans$success,
        G_condition_number = unname(cond[["condition_number"]]),
        G_min_singular_value = unname(cond[["min_singular_value"]]),
        G_max_singular_value = unname(cond[["max_singular_value"]]),
        stringsAsFactors = FALSE
      )

      if (!ans$success) next
      names(ans$impact) <- var_names

      for (cc in countries) {
        vv <- paste0(cc, "_y")
        impact_rows[[length(impact_rows) + 1L]] <- data.frame(
          date = date,
          draw = draw,
          scenario = scenario,
          country = cc,
          variable = vv,
          h0_gdp = ans$impact[[vv]],
          positive = ans$impact[[vv]] > 0,
          stringsAsFactors = FALSE
        )
      }
    }
  }
}

if (!length(impact_rows)) stop("No scenario impact rows were produced.")

impact_draws <- do.call(rbind, impact_rows)
cond_draws <- do.call(rbind, cond_rows)
gdp_coeff_draws <- do.call(rbind, gdp_coeff_rows)
gpr_all_eq_draws <- do.call(rbind, gpr_all_eq_rows)

write.csv(
  impact_draws,
  file.path(out_dir, "h0_channel_scenario_draws.csv"),
  row.names = FALSE
)

write.csv(
  gdp_coeff_draws,
  file.path(out_dir, "nonUS_gdp_wex_coeff_draws.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

impact_split <- split(
  impact_draws,
  interaction(
    impact_draws$date,
    impact_draws$scenario,
    impact_draws$country,
    drop = TRUE,
    lex.order = TRUE
  )
)

impact_summary <- do.call(
  rbind,
  lapply(impact_split, function(z) {
    s <- summarize_vec(z$h0_gdp)
    data.frame(
      date = z$date[[1]],
      scenario = z$scenario[[1]],
      country = z$country[[1]],
      draws = nrow(z),
      median = s[["median"]],
      low68 = s[["low68"]],
      high68 = s[["high68"]],
      low90 = s[["low90"]],
      high90 = s[["high90"]],
      positive_share = mean(z$positive),
      credible68_excludes_zero =
        is.finite(s[["low68"]]) && is.finite(s[["high68"]]) &&
        (s[["low68"]] > 0 || s[["high68"]] < 0),
      stringsAsFactors = FALSE
    )
  })
)

impact_summary <- impact_summary[
  order(
    match(impact_summary$date, available_dates),
    match(impact_summary$scenario, scenarios),
    match(impact_summary$country, countries)
  ),
  ,
  drop = FALSE
]

write.csv(
  impact_summary,
  file.path(out_dir, "h0_channel_scenario_summary.csv"),
  row.names = FALSE
)

# Sign counts by date/scenario.
sign_count <- do.call(
  rbind,
  lapply(
    split(
      impact_summary,
      interaction(
        impact_summary$date,
        impact_summary$scenario,
        drop = TRUE,
        lex.order = TRUE
      )
    ),
    function(z) {
      data.frame(
        date = z$date[[1]],
        scenario = z$scenario[[1]],
        countries = nrow(z),
        positive_median_count = sum(z$median > 0, na.rm = TRUE),
        negative_median_count = sum(z$median < 0, na.rm = TRUE),
        near_zero_median_count = sum(abs(z$median) <= 1e-10, na.rm = TRUE),
        positive_68_excludes_zero = sum(
          z$median > 0 & z$credible68_excludes_zero,
          na.rm = TRUE
        ),
        negative_68_excludes_zero = sum(
          z$median < 0 & z$credible68_excludes_zero,
          na.rm = TRUE
        ),
        median_abs_gdp = median(abs(z$median), na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
  )
)

sign_count <- sign_count[
  order(
    match(sign_count$date, available_dates),
    match(sign_count$scenario, scenarios)
  ),
  ,
  drop = FALSE
]

write.csv(
  sign_count,
  file.path(out_dir, "h0_sign_counts_by_scenario.csv"),
  row.names = FALSE
)

# Non-US GDP equation Wex coefficient summary.
gdp_coeff_split <- split(
  gdp_coeff_draws,
  interaction(
    gdp_coeff_draws$date,
    gdp_coeff_draws$country,
    gdp_coeff_draws$term,
    drop = TRUE,
    lex.order = TRUE
  )
)

gdp_coeff_summary <- do.call(
  rbind,
  lapply(gdp_coeff_split, function(z) {
    s <- summarize_vec(z$coefficient)
    data.frame(
      date = z$date[[1]],
      country = z$country[[1]],
      term = z$term[[1]],
      draws = nrow(z),
      median = s[["median"]],
      low68 = s[["low68"]],
      high68 = s[["high68"]],
      low90 = s[["low90"]],
      high90 = s[["high90"]],
      positive_share = mean(z$coefficient > 0),
      stringsAsFactors = FALSE
    )
  })
)

gdp_coeff_summary <- gdp_coeff_summary[
  order(
    match(gdp_coeff_summary$date, available_dates),
    match(gdp_coeff_summary$country, countries),
    match(
      gdp_coeff_summary$term,
      c(
        "foreign_y", "foreign_dp", "foreign_r",
        "foreign_de", "foreign_deq",
        "global_gpr", "global_vix"
      )
    )
  ),
  ,
  drop = FALSE
]

write.csv(
  gdp_coeff_summary,
  file.path(out_dir, "nonUS_gdp_wex_coeff_summary.csv"),
  row.names = FALSE
)

# Global-GPR coefficient summary for every non-US equation.
gpr_eq_split <- split(
  gpr_all_eq_draws,
  interaction(
    gpr_all_eq_draws$date,
    gpr_all_eq_draws$country,
    gpr_all_eq_draws$equation,
    drop = TRUE,
    lex.order = TRUE
  )
)

gpr_eq_summary <- do.call(
  rbind,
  lapply(gpr_eq_split, function(z) {
    s <- summarize_vec(z$global_gpr_coefficient)
    data.frame(
      date = z$date[[1]],
      country = z$country[[1]],
      equation = z$equation[[1]],
      draws = nrow(z),
      median = s[["median"]],
      low68 = s[["low68"]],
      high68 = s[["high68"]],
      low90 = s[["low90"]],
      high90 = s[["high90"]],
      positive_share = mean(z$global_gpr_coefficient > 0),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  gpr_eq_summary,
  file.path(out_dir, "nonUS_global_gpr_coeff_summary_all_equations.csv"),
  row.names = FALSE
)

# G-condition summary by date/scenario.
cond_split <- split(
  cond_draws,
  interaction(
    cond_draws$date,
    cond_draws$scenario,
    drop = TRUE,
    lex.order = TRUE
  )
)

cond_summary <- do.call(
  rbind,
  lapply(cond_split, function(z) {
    finite_cond <- z$G_condition_number[is.finite(z$G_condition_number)]
    data.frame(
      date = z$date[[1]],
      scenario = z$scenario[[1]],
      draws = nrow(z),
      solve_success_share = mean(z$solve_success),
      median_condition = if (length(finite_cond)) median(finite_cond) else Inf,
      p95_condition = if (length(finite_cond)) qv(finite_cond, 0.95) else Inf,
      median_min_singular = median(z$G_min_singular_value, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

cond_summary <- cond_summary[
  order(
    match(cond_summary$date, available_dates),
    match(cond_summary$scenario, scenarios)
  ),
  ,
  drop = FALSE
]

write.csv(
  cond_summary,
  file.path(out_dir, "scenario_G_condition_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Human-readable comparison against baseline pure-GPR response
# -----------------------------------------------------------------------------

baseline <- impact_summary[
  impact_summary$scenario == "A_baseline",
  c("date", "country", "median")
]
names(baseline)[3] <- "baseline_median"

comparison <- merge(
  impact_summary,
  baseline,
  by = c("date", "country"),
  all.x = TRUE,
  sort = FALSE
)

comparison$median_ratio_to_baseline <- ifelse(
  abs(comparison$baseline_median) > 1e-12,
  comparison$median / comparison$baseline_median,
  NA_real_
)
comparison$sign_changed_vs_baseline <-
  sign(comparison$median) != sign(comparison$baseline_median)

write.csv(
  comparison,
  file.path(out_dir, "h0_scenario_vs_baseline_comparison.csv"),
  row.names = FALSE
)

# Date/scenario magnitude ratios.
ratio_lines <- character(0)
for (date in available_dates) {
  base_date <- impact_summary[
    impact_summary$date == date &
      impact_summary$scenario == "A_baseline",
    ,
    drop = FALSE
  ]
  base_mag <- median(abs(base_date$median), na.rm = TRUE)

  for (scenario in scenarios) {
    z <- impact_summary[
      impact_summary$date == date &
        impact_summary$scenario == scenario,
      ,
      drop = FALSE
    ]
    zmag <- median(abs(z$median), na.rm = TRUE)
    ratio <- if (is.finite(base_mag) && base_mag > 1e-12) zmag / base_mag else NA_real_

    ratio_lines <- c(
      ratio_lines,
      paste0(
        date, " | ", scenario,
        " | positive medians=", sum(z$median > 0, na.rm = TRUE), "/", nrow(z),
        " | negative medians=", sum(z$median < 0, na.rm = TRUE), "/", nrow(z),
        " | median |GDP h0| ratio vs baseline=",
        ifelse(is.finite(ratio), format(round(ratio, 4), nsmall = 4), "NA")
      )
    )
  }
}

summary_lines <- c(
  "P2-Strong GPR + VIX: contemporaneous channel decomposition",
  paste0("Posterior file: ", pred_file),
  paste0("IRF file: ", irf_file),
  paste0("Retained posterior draws: ", n_draws),
  paste0("Stable-only draw selection: ", stable_only),
  paste0("Selected dates: ", paste(available_dates, collapse = ", ")),
  paste0("GPR normalization: +", shock_pct, "%"),
  paste0(
    "Maximum baseline full-IRF reproduction error: ",
    format(max_baseline_reproduction_error, scientific = TRUE, digits = 6)
  ),
  "",
  "SCENARIOS",
  "A_baseline: original G, pure GPR residual shock.",
  "B1_zero_nonUS_GPR_GDP: zero current global-GPR only in non-US GDP equations.",
  "B2_zero_nonUS_GPR_all: zero current global-GPR in every non-US equation.",
  "D_zero_foreign_macro: zero contemporaneous foreign macro coefficients everywhere; retain non-US global GPR/VIX.",
  "C_zero_GPR_all_plus_foreign_macro: B2 plus zero contemporaneous foreign macro coefficients everywhere.",
  "",
  "INTERPRETATION",
  "If B1 sharply removes positive GDP responses, the direct non-US GDP-equation global-GPR coefficient is the main source.",
  "If B1 does little but B2 sharply removes them, indirect contemporaneous GPR effects through non-GDP equations are important.",
  "If B2 does little but C sharply removes them, contemporaneous foreign-macro feedback is the main source.",
  "If D remains positive while B2/C collapse, direct non-US current global-GPR exposure is dominant.",
  "If B2 and C still leave systematic positive GDP responses, inspect W/index mapping and the construction of G itself.",
  "",
  "DATE / SCENARIO RESULTS",
  ratio_lines
)

writeLines(
  summary_lines,
  file.path(out_dir, "h0_channel_decomposition_summary.txt")
)

cat("\n")
cat(paste(summary_lines, collapse = "\n"))
cat("\n\nCreated diagnostic files in: ", out_dir, "\n", sep = "")
