#!/usr/bin/env Rscript

rm(list = ls())
seed <- as.integer(Sys.getenv("TVPGVAR_SEED", "20260816"))
if (length(seed) != 1L || is.na(seed)) stop("TVPGVAR_SEED must be an integer.")
RNGkind("L'Ecuyer-CMRG")
set.seed(seed)

suppressPackageStartupMessages({
  library(compiler)
  library(snowfall)
  library(Matrix)
  library(mvtnorm)
  library(threshtvp)
  library(ggplot2)
})

source("R/BVAR_ttvp_gprlag.r")
source("R/Datahandling.r")
source("R/auxilliary_functions_tvp.r")
source("R/prepare_data_structural.R")
source("R/gpr_structural_irf.R")

dir.create("results", showWarnings = FALSE, recursive = TRUE)

Data.setup <- prepare_gvar_data()
xglobal <- Data.setup$bigx
gW <- Data.setup$gW
Daten <- Data.setup$new.data
cN <- Data.setup$countries

# Confirm the US endogenous block is exactly GPR-first.
us_cols <- colnames(xglobal)[substr(colnames(xglobal), 1, 2) == "US"]
expected_us <- c("US_gpr", "US_y", "US_dp", "US_r", "US_de", "US_deq")
if (!identical(us_cols, expected_us)) {
  stop("US recursive order mismatch. Found: ", paste(us_cols, collapse = ", "))
}

CPU <- max(1L, min(4L, parallel::detectCores() - 1L))
saves <- as.integer(Sys.getenv("TVPGVAR_SAVES", "100"))
burns <- as.integer(Sys.getenv("TVPGVAR_BURNS", "100"))
thin <- as.numeric(Sys.getenv("TVPGVAR_THIN", "0.5"))
nhor <- as.integer(Sys.getenv("TVPGVAR_HORIZON", "12"))
shock_pct <- as.numeric(Sys.getenv("TVPGVAR_GPR_SHOCK_PCT", "10"))
gpr_lag_mode <- Sys.getenv("TVPGVAR_GPR_LAG_MODE", "current_only")
lag_order <- as.integer(Sys.getenv("TVPGVAR_P", "1"))
near_unit_threshold <- as.numeric(Sys.getenv("TVPGVAR_NEAR_UNIT_THRESHOLD", "0.98"))
ext.inst <- FALSE
shrink.parm <- list(B_1 = 3, B_2 = 0.03, kappa0 = -0.1/20)

if (!is.finite(saves) || saves < 10) stop("TVPGVAR_SAVES is too small.")
if (!is.finite(burns) || burns < 10) stop("TVPGVAR_BURNS is too small.")
if (!is.finite(thin) || thin <= 0 || thin > 1) stop("TVPGVAR_THIN must lie in (0,1].")
if (!is.finite(shock_pct) || shock_pct <= 0) stop("TVPGVAR_GPR_SHOCK_PCT must be positive.")
if (!is.finite(near_unit_threshold) || near_unit_threshold <= 0 || near_unit_threshold >= 1) {
  stop("TVPGVAR_NEAR_UNIT_THRESHOLD must lie in (0,1).")
}
if (!gpr_lag_mode %in% c("current_only", "current_and_lag")) {
  stop("TVPGVAR_GPR_LAG_MODE must be current_only or current_and_lag.")
}
if (length(lag_order) != 1L || is.na(lag_order) || !lag_order %in% c(1L, 2L)) {
  stop("TVPGVAR_P must be 1 or 2.")
}

# Deterministic independent random-number stream for each country.  This makes
# grid cells comparable even when snowfall schedules workers differently.
rng_streams <- vector("list", length(cN))
rng_streams[[1L]] <- .Random.seed
if (length(cN) > 1L) {
  for (ii in 2:length(cN)) {
    rng_streams[[ii]] <- parallel::nextRNGStream(rng_streams[[ii - 1L]])
  }
}

BVAR <- cmpfun(BVAR)

sfInit(parallel = TRUE, cpus = CPU)
on.exit(try(sfStop(), silent = TRUE), add = TRUE)
sfExport(list = list(
  "mlag", "BVAR", "datahandling", "xglobal", "gW", "Daten", "cN",
  "bvartvpm", "saves", "burns", "thin", "ext.inst", "shrink.parm",
  "gpr_lag_mode", "tvpgvar_wex_lag", "rng_streams"
))

predDens <- sfLapply(seq_along(cN), function(i) {
  assign(".Random.seed", rng_streams[[i]], envir = .GlobalEnv)
  BVAR(i, gW = gW, bigx = xglobal, Daten = Daten, cN = cN,
       nsave = saves, nburn = burns, thin_chain = thin,
       ext.inst = ext.inst, parms = shrink.parm)
})
sfStop()

save(predDens, Data.setup, file = "results/predDens_gpr_structural.rda")


# =============================================================================
# MODEL-VALIDATION DIAGNOSTICS
# =============================================================================
# These diagnostics follow the validation spirit of the GVAR literature, but
# the feedback test below is deliberately NOT labelled the formal cointegration-
# based weak-exogeneity test of Dées et al. because the present model uses
# transformed/stationary macro-financial variables rather than a VECM.

# 1) Residual serial correlation: Ljung-Box test at four quarterly lags.
residual_serial_diag <- do.call(rbind, lapply(seq_along(predDens), function(i) {
  cc <- cN[[i]]
  res <- as.matrix(predDens[[i]]$cc.res)
  own_cols <- colnames(xglobal)[substr(colnames(xglobal), 1, 2) == cc]

  if (ncol(res) != length(own_cols)) {
    stop("Residual-column mismatch for ", cc,
         ": residual columns=", ncol(res), ", own variables=", length(own_cols))
  }

  do.call(rbind, lapply(seq_along(own_cols), function(j) {
    z <- as.numeric(res[, j])
    z <- z[is.finite(z)]
    lag_use <- min(4L, max(1L, floor(length(z) / 5L)))
    bt <- tryCatch(
      Box.test(z, lag = lag_use, type = "Ljung-Box", fitdf = 0),
      error = function(e) NULL
    )
    data.frame(
      country = cc,
      variable = own_cols[[j]],
      n = length(z),
      test_lag = lag_use,
      statistic = if (is.null(bt)) NA_real_ else unname(bt$statistic),
      p_value = if (is.null(bt)) NA_real_ else bt$p.value,
      reject_5pct = if (is.null(bt)) NA else bt$p.value < 0.05,
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(residual_serial_diag,
          "results/residual_serial_correlation_diagnostic.csv",
          row.names = FALSE)

# 2) Pairwise cross-unit residual correlation for common macro variables.
common_suffixes <- c("y", "dp", "r", "de", "deq")
residual_cross_diag <- list()

for (vv in common_suffixes) {
  series <- list()
  for (i in seq_along(predDens)) {
    cc <- cN[[i]]
    own_cols <- colnames(xglobal)[substr(colnames(xglobal), 1, 2) == cc]
    target <- paste0(cc, "_", vv)
    j <- match(target, own_cols)
    if (!is.na(j)) series[[cc]] <- as.numeric(predDens[[i]]$cc.res[, j])
  }

  if (length(series) >= 2L) {
    pairs <- combn(names(series), 2, simplify = FALSE)
    residual_cross_diag[[vv]] <- do.call(rbind, lapply(pairs, function(pp) {
      x1 <- series[[pp[[1]]]]
      x2 <- series[[pp[[2]]]]
      ok <- is.finite(x1) & is.finite(x2)
      x1 <- x1[ok]
      x2 <- x2[ok]
      ct <- tryCatch(
        suppressWarnings(cor.test(x1, x2, method = "pearson")),
        error = function(e) NULL
      )
      data.frame(
        variable = vv,
        country_1 = pp[[1]],
        country_2 = pp[[2]],
        n = length(x1),
        correlation = if (is.null(ct)) NA_real_ else unname(ct$estimate),
        p_value = if (is.null(ct)) NA_real_ else ct$p.value,
        significant_5pct = if (is.null(ct)) NA else ct$p.value < 0.05,
        stringsAsFactors = FALSE
      )
    }))
  }
}
residual_cross_diag <- do.call(rbind, residual_cross_diag)
write.csv(residual_cross_diag,
          "results/residual_cross_unit_correlation_diagnostic.csv",
          row.names = FALSE)

residual_cross_summary <- do.call(rbind, lapply(split(residual_cross_diag,
                                                      residual_cross_diag$variable),
                                                function(z) {
  a <- abs(z$correlation[is.finite(z$correlation)])
  data.frame(
    variable = unique(z$variable)[1],
    pairs = nrow(z),
    median_abs_correlation = if (length(a)) median(a) else NA_real_,
    p95_abs_correlation = if (length(a)) as.numeric(quantile(a, 0.95)) else NA_real_,
    max_abs_correlation = if (length(a)) max(a) else NA_real_,
    significant_share_5pct = mean(z$significant_5pct, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))
write.csv(residual_cross_summary,
          "results/residual_cross_unit_correlation_summary.csv",
          row.names = FALSE)

# 3) Stationary feedback / exogeneity diagnostic.
# H0: lagged domestic variables do not add predictive power for a foreign/global
# aggregate after controlling for that aggregate's own first lag.
# This is a practical Granger-feedback diagnostic for the transformed variables,
# NOT the formal error-correction weak-exogeneity test used in levels/VECM GVARs.
feedback_diag <- do.call(rbind, lapply(seq_along(cN), function(i) {
  cc <- cN[[i]]
  End <- xglobal[, substr(colnames(xglobal), 1, 2) == cc, drop = FALSE]
  W <- predDens[[i]]$W

  weighted_all <- W %*% t(xglobal)
  Wex <- t(weighted_all[(ncol(End) + 1L):nrow(weighted_all), , drop = FALSE])

  expected_labels <- if (cc == "US") {
    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq")
  } else {
    c("foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq",
      "global_gpr")
  }

  if (ncol(Wex) != length(expected_labels)) {
    stop("Unexpected Wex dimension in feedback diagnostic for ", cc,
         ": found ", ncol(Wex), ", expected ", length(expected_labels))
  }
  colnames(Wex) <- expected_labels

  do.call(rbind, lapply(seq_len(ncol(Wex)), function(j) {
    y <- Wex[-1, j]
    wlag <- Wex[-nrow(Wex), j]
    domlag <- End[-nrow(End), , drop = FALSE]
    colnames(domlag) <- paste0("domlag_", make.names(colnames(domlag)))

    df <- data.frame(y = y, wlag = wlag, domlag, check.names = FALSE)
    df <- df[complete.cases(df), , drop = FALSE]

    ans <- tryCatch({
      fit_r <- lm(y ~ wlag, data = df)
      fit_u <- lm(y ~ ., data = df)
      aa <- anova(fit_r, fit_u)
      list(
        F = aa$F[2],
        p = aa$`Pr(>F)`[2],
        df_num = aa$Df[2],
        df_den = df.residual(fit_u)
      )
    }, error = function(e) NULL)

    data.frame(
      country = cc,
      foreign_or_global_variable = expected_labels[[j]],
      n = nrow(df),
      F_statistic = if (is.null(ans)) NA_real_ else ans$F,
      p_value = if (is.null(ans)) NA_real_ else ans$p,
      df_num = if (is.null(ans)) NA_real_ else ans$df_num,
      df_den = if (is.null(ans)) NA_real_ else ans$df_den,
      reject_no_feedback_5pct = if (is.null(ans)) NA else ans$p < 0.05,
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(feedback_diag,
          "results/stationary_feedback_exogeneity_diagnostic.csv",
          row.names = FALSE)

feedback_summary <- do.call(rbind, lapply(split(feedback_diag, feedback_diag$country),
                                          function(z) {
  data.frame(
    country = unique(z$country)[1],
    tests = nrow(z),
    rejection_share_5pct = mean(z$reject_no_feedback_5pct, na.rm = TRUE),
    min_p_value = suppressWarnings(min(z$p_value, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}))
feedback_summary$min_p_value[!is.finite(feedback_summary$min_p_value)] <- NA_real_
write.csv(feedback_summary,
          "results/stationary_feedback_exogeneity_summary.csv",
          row.names = FALSE)

writeLines(c(
  "Validation notes",
  "1. residual_serial_correlation_diagnostic.csv: Ljung-Box tests on posterior-mean residuals.",
  "2. residual_cross_unit_correlation_diagnostic.csv: pairwise Pearson residual correlations across countries.",
  "3. stationary_feedback_exogeneity_diagnostic.csv: stationary Granger-feedback diagnostic.",
  "   It is NOT the formal error-correction weak-exogeneity test used in VECM-based GVARs.",
  "4. Formal MCMC convergence diagnostics are only treated as interpretable when >=100 retained draws are available."
), "results/model_validation_notes.txt")

# Confirm that the estimated regressor structure matches the requested GPR lag mode.
reg_diag <- do.call(rbind, lapply(seq_along(predDens), function(i) {
  rn <- dimnames(predDens[[i]]$ALPHA)[[2]]
  data.frame(
    country = cN[[i]],
    Wex = sum(rn == "Wex"),
    Wexlag1 = sum(rn == "Wexlag1"),
    Wexlag2 = sum(rn == "Wexlag2"),
    Ylag1 = sum(rn == "Ylag1"),
    Ylag2 = sum(rn == "Ylag2"),
    stringsAsFactors = FALSE
  )
}))

expected_lag <- if (gpr_lag_mode == "current_only") {
  rep(5L, length(cN))
} else {
  ifelse(cN == "US", 5L, 6L)
}
expected_wex <- ifelse(cN == "US", 5L, 6L)
expected_ylag <- ifelse(cN == "US", 6L, 5L)
expected_ylag2 <- if (lag_order == 2L) expected_ylag else rep(0L, length(cN))

if (any(reg_diag$Wex != expected_wex) ||
    any(reg_diag$Wexlag1 != expected_lag) ||
    any(reg_diag$Wexlag2 != 0L) ||
    any(reg_diag$Ylag1 != expected_ylag) ||
    any(reg_diag$Ylag2 != expected_ylag2)) {
  print(reg_diag)
  stop("Estimated ALPHA regressor structure does not match requested GPR lag specification.")
}
reg_diag$gpr_lag_mode <- gpr_lag_mode
reg_diag$domestic_lag_order <- lag_order
reg_diag$foreign_lag_order <- 1L
write.csv(reg_diag, "results/gpr_lag_spec_diagnostic.csv", row.names = FALSE)

Sigma.posterior <- A.list <- globalG <- vector("list", length(predDens))
for (i in seq_along(predDens)) {
  globalG[[i]] <- predDens[[i]]$W
  Sigma.posterior[[i]] <- predDens[[i]]$SIGMApost
  A.list[[i]] <- predDens[[i]]$ALPHA
}

n_irf_draws <- dim(A.list[[1]])[4]
if (is.null(n_irf_draws) || n_irf_draws < 1L) stop("No retained posterior draws found.")
if (any(vapply(A.list, function(z) dim(z)[4], integer(1)) != n_irf_draws)) {
  stop("Countries have inconsistent retained posterior draw counts.")
}

Tirf <- nrow(xglobal) - lag_order
irf_dates <- Data.setup$quarters[-seq_len(lag_order)]
if (any(vapply(A.list, function(z) dim(z)[1], integer(1)) != Tirf)) {
  stop("Estimated coefficient-path length does not match the requested lag order.")
}

IRF_post <- array(
  NA_real_,
  c(Tirf, ncol(xglobal), nhor + 1L, n_irf_draws),
  dimnames = list(irf_dates, colnames(xglobal), 0:nhor, NULL)
)

# rho[t,draw] = maximum modulus of companion-matrix eigenvalues.
stability_rho <- matrix(
  NA_real_, nrow = Tirf, ncol = n_irf_draws,
  dimnames = list(irf_dates, paste0("draw_", seq_len(n_irf_draws)))
)

G_condition <- matrix(
  NA_real_, nrow = Tirf, ncol = n_irf_draws,
  dimnames = dimnames(stability_rho)
)
G_min_singular <- matrix(
  NA_real_, nrow = Tirf, ncol = n_irf_draws,
  dimnames = dimnames(stability_rho)
)
G_max_singular <- matrix(
  NA_real_, nrow = Tirf, ncol = n_irf_draws,
  dimnames = dimnames(stability_rho)
)
max_abs_transition <- matrix(
  NA_real_, nrow = Tirf, ncol = n_irf_draws,
  dimnames = dimnames(stability_rho)
)
local_rho_array <- array(
  NA_real_,
  dim = c(Tirf, length(cN), n_irf_draws),
  dimnames = list(irf_dates, cN, paste0("draw_", seq_len(n_irf_draws)))
)

for (irep in seq_len(n_irf_draws)) {
  A.i <- lapply(A.list, function(z) z[, , , irep])
  S.i <- lapply(Sigma.posterior, function(z) z[, , , irep])

  for (tt in seq_len(Tirf)) {
    ans <- get_gpr_struct_irfa_t(
      tt = tt, draw_i = A.i, Sig_draw_i = S.i,
      x = t(xglobal), globalG = globalG, countries = cN,
      horz = nhor, normalization = "pct", shock_pct = shock_pct
    )
    IRF_post[tt, , , irep] <- ans$IRF_post
    stability_rho[tt, irep] <- ans$max_eigen_modulus
    G_condition[tt, irep] <- ans$G_condition_number
    G_min_singular[tt, irep] <- ans$G_min_singular_value
    G_max_singular[tt, irep] <- ans$G_max_singular_value
    max_abs_transition[tt, irep] <- ans$max_abs_transition
    local_rho_array[tt, , irep] <- ans$local_rho[cN]
  }
  cat("IRF/stability posterior draw", irep, "of", n_irf_draws, "complete\n")
}

if (any(!is.finite(stability_rho))) {
  stop("Non-finite spectral radius detected in stability audit.")
}

if (any(!is.finite(G_min_singular)) || any(!is.finite(G_max_singular)) ||
    any(!is.finite(max_abs_transition)) || any(!is.finite(local_rho_array))) {
  stop("Non-finite source-diagnostic value detected in stability audit.")
}
# Infinite condition numbers are allowed and explicitly classified as ill-conditioned.

stable_mask <- stability_rho < 1
near_unit_mask <- stability_rho >= near_unit_threshold & stability_rho < 1
unstable_mask <- stability_rho >= 1

# Save all draws plus stability diagnostics.
save(IRF_post, stability_rho, stable_mask,
     file = "results/irf_gpr_structural.rda")

# Stable-only IRF object: unstable date/draw slices are set to NA.
# Near-unit-but-still-stable draws are NOT discarded; they are only flagged.
IRF_post_stable <- IRF_post
for (tt in seq_len(Tirf)) {
  bad <- which(!stable_mask[tt, ])
  if (length(bad)) {
    for (ddraw in bad) IRF_post_stable[tt, , , ddraw] <- NA_real_
  }
}
save(IRF_post_stable, stability_rho, stable_mask,
     file = "results/irf_gpr_structural_stable_only.rda")

# Long-form stability diagnostic.
stability_detail <- do.call(rbind, lapply(seq_len(n_irf_draws), function(ddraw) {
  data.frame(
    date = irf_dates,
    draw = ddraw,
    max_eigen_modulus = stability_rho[, ddraw],
    stable = stable_mask[, ddraw],
    near_unit = near_unit_mask[, ddraw],
    unstable = unstable_mask[, ddraw],
    stringsAsFactors = FALSE
  )
}))
write.csv(stability_detail, "results/stability_diagnostic.csv", row.names = FALSE)

# Date-level stability summary.
stability_summary <- do.call(rbind, lapply(seq_len(Tirf), function(tt) {
  r <- stability_rho[tt, ]
  data.frame(
    date = irf_dates[tt],
    draws = length(r),
    stable_draws = sum(r < 1),
    stable_share = mean(r < 1),
    near_unit_draws = sum(r >= near_unit_threshold & r < 1),
    near_unit_share = mean(r >= near_unit_threshold & r < 1),
    unstable_draws = sum(r >= 1),
    unstable_share = mean(r >= 1),
    median_rho = median(r),
    p95_rho = as.numeric(quantile(r, 0.95, names = FALSE)),
    max_rho = max(r),
    stringsAsFactors = FALSE
  )
}))
write.csv(stability_summary, "results/stability_summary_by_date.csv", row.names = FALSE)


# =============================================================================
# SOURCE-OF-INSTABILITY DIAGNOSTICS
# =============================================================================
G_condition_warn <- as.numeric(Sys.getenv("TVPGVAR_G_CONDITION_WARN", "10000"))
if (!is.finite(G_condition_warn) || G_condition_warn <= 1) G_condition_warn <- 10000

local_stability_detail <- do.call(rbind, lapply(seq_len(n_irf_draws), function(ddraw) {
  do.call(rbind, lapply(seq_along(cN), function(i) {
    data.frame(
      date = irf_dates,
      draw = ddraw,
      country = cN[[i]],
      local_rho = local_rho_array[, i, ddraw],
      local_stable = local_rho_array[, i, ddraw] < 1,
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(local_stability_detail,
          "results/local_stability_diagnostic.csv",
          row.names = FALSE)

local_stability_summary_country <- do.call(rbind, lapply(seq_along(cN), function(i) {
  r <- as.numeric(local_rho_array[, i, ])
  data.frame(
    country = cN[[i]],
    observations = length(r),
    stable_share = mean(r < 1),
    median_rho = median(r),
    p95_rho = as.numeric(quantile(r, 0.95)),
    max_rho = max(r),
    stringsAsFactors = FALSE
  )
}))
write.csv(local_stability_summary_country,
          "results/local_stability_summary_by_country.csv",
          row.names = FALSE)

# Coefficient-level diagnosis for the US block.  This separates the persistent
# short-rate equation from the GPR equation and records whether removing the
# short-rate row would move an otherwise unstable local system inside the unit
# circle.  The row-removal result is diagnostic only; it is not used to alter
# estimation or the reported IRFs.
us_i <- match("US", cN)
if (is.na(us_i)) stop("US block not found for coefficient-level diagnosis.")
us_A <- A.list[[us_i]]
us_reg_names <- dimnames(us_A)[[2]]
us_eq_names <- dimnames(us_A)[[3]]
if (is.null(us_reg_names) || is.null(us_eq_names)) {
  stop("US coefficient array lacks regressor or equation names.")
}
us_r_i <- match("US_r", us_eq_names)
us_gpr_i <- match("US_gpr", us_eq_names)
if (is.na(us_r_i) || is.na(us_gpr_i)) {
  stop("US_r or US_gpr was not found in the US equation block.")
}

extract_us_lag_matrices <- function(tt, ddraw) {
  lapply(seq_len(lag_order), function(kk) {
    idx <- which(us_reg_names == paste0("Ylag", kk))
    if (length(idx) != length(us_eq_names)) {
      stop("Unexpected US domestic lag-block width at lag ", kk, ".")
    }
    z <- us_A[tt, idx, , ddraw, drop = TRUE]
    z <- matrix(z, nrow = length(idx), ncol = length(us_eq_names))
    theta <- t(z)
    dimnames(theta) <- list(us_eq_names, us_eq_names)
    theta
  })
}

us_coefficient_diag <- do.call(rbind, lapply(seq_len(n_irf_draws), function(ddraw) {
  do.call(rbind, lapply(seq_len(Tirf), function(tt) {
    theta <- extract_us_lag_matrices(tt, ddraw)
    theta_no_r <- lapply(theta, function(z) {
      z[us_r_i, ] <- 0
      z
    })
    theta_no_gpr <- lapply(theta, function(z) {
      z[us_gpr_i, ] <- 0
      z
    })

    r_own <- vapply(theta, function(z) z[us_r_i, us_r_i], numeric(1))
    gpr_own <- vapply(theta, function(z) z[us_gpr_i, us_gpr_i], numeric(1))
    rho_full <- local_rho_array[tt, us_i, ddraw]
    rho_no_r <- tvpgvar_spectral_radius(theta_no_r)
    rho_no_gpr <- tvpgvar_spectral_radius(theta_no_gpr)

    data.frame(
      date = irf_dates[[tt]],
      draw = ddraw,
      lag_order = lag_order,
      US_local_rho = rho_full,
      US_local_stable = rho_full < 1,
      US_r_ar1 = r_own[[1]],
      US_r_ar2 = if (lag_order >= 2L) r_own[[2]] else 0,
      US_r_ar_sum = sum(r_own),
      US_gpr_ar1 = gpr_own[[1]],
      US_gpr_ar2 = if (lag_order >= 2L) gpr_own[[2]] else 0,
      US_gpr_ar_sum = sum(gpr_own),
      rho_without_US_r_row = rho_no_r,
      rho_without_US_gpr_row = rho_no_gpr,
      rho_reduction_without_US_r = rho_full - rho_no_r,
      rho_reduction_without_US_gpr = rho_full - rho_no_gpr,
      US_r_row_removal_stabilizes =
        rho_full >= 1 && rho_no_r < 1,
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(us_coefficient_diag,
          "results/us_local_coefficient_diagnostic.csv",
          row.names = FALSE)

us_coefficient_summary <- data.frame(
  lag_order = lag_order,
  observations = nrow(us_coefficient_diag),
  US_stable_share = mean(us_coefficient_diag$US_local_stable),
  median_US_r_ar1 = median(us_coefficient_diag$US_r_ar1),
  median_US_r_ar2 = median(us_coefficient_diag$US_r_ar2),
  median_US_r_ar_sum = median(us_coefficient_diag$US_r_ar_sum),
  p95_US_r_ar_sum = as.numeric(quantile(us_coefficient_diag$US_r_ar_sum, 0.95)),
  median_US_gpr_ar_sum = median(us_coefficient_diag$US_gpr_ar_sum),
  median_rho_reduction_without_US_r =
    median(us_coefficient_diag$rho_reduction_without_US_r),
  unstable_share_fixed_by_US_r_row =
    if (any(!us_coefficient_diag$US_local_stable)) {
      mean(us_coefficient_diag$US_r_row_removal_stabilizes[
        !us_coefficient_diag$US_local_stable
      ])
    } else {
      NA_real_
    },
  stringsAsFactors = FALSE
)
write.csv(us_coefficient_summary,
          "results/us_local_coefficient_summary.csv",
          row.names = FALSE)

source_diag <- do.call(rbind, lapply(seq_len(n_irf_draws), function(ddraw) {
  do.call(rbind, lapply(seq_len(Tirf), function(tt) {
    lr <- local_rho_array[tt, , ddraw]
    worst_i <- which.max(lr)
    global_unstable <- stability_rho[tt, ddraw] >= 1
    local_unstable <- max(lr) >= 1
    ill_G <- !is.finite(G_condition[tt, ddraw]) ||
      G_condition[tt, ddraw] >= G_condition_warn

    likely_source <- if (!global_unstable) {
      "global_stable"
    } else if (local_unstable && ill_G) {
      "local_dynamics_and_G_conditioning"
    } else if (local_unstable) {
      "local_dynamics"
    } else if (ill_G) {
      "near_singular_G"
    } else {
      "cross_country_feedback"
    }

    data.frame(
      date = irf_dates[[tt]],
      draw = ddraw,
      global_rho = stability_rho[tt, ddraw],
      G_condition_number = G_condition[tt, ddraw],
      G_min_singular_value = G_min_singular[tt, ddraw],
      G_max_singular_value = G_max_singular[tt, ddraw],
      max_abs_transition = max_abs_transition[tt, ddraw],
      max_local_rho = max(lr),
      worst_local_country = cN[[worst_i]],
      global_unstable = global_unstable,
      any_local_unstable = local_unstable,
      G_ill_conditioned = ill_G,
      likely_source = likely_source,
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(source_diag,
          "results/global_stability_source_diagnostic.csv",
          row.names = FALSE)

source_summary <- as.data.frame(table(source_diag$likely_source),
                                stringsAsFactors = FALSE)
names(source_summary) <- c("likely_source", "count")
source_summary$share <- source_summary$count / sum(source_summary$count)
source_summary <- source_summary[order(-source_summary$share), ]
write.csv(source_summary,
          "results/stability_source_summary.csv",
          row.names = FALSE)

G_condition_summary <- do.call(rbind, lapply(seq_len(Tirf), function(tt) {
  z <- G_condition[tt, ]
  data.frame(
    date = irf_dates[[tt]],
    median_condition_number = median(z[is.finite(z)]),
    p95_condition_number = if (any(is.finite(z))) as.numeric(quantile(z[is.finite(z)], 0.95)) else Inf,
    max_condition_number = max(z),
    ill_conditioned_share = mean(!is.finite(z) | z >= G_condition_warn),
    median_min_singular_value = median(G_min_singular[tt, ]),
    p05_min_singular_value = as.numeric(quantile(G_min_singular[tt, ], 0.05)),
    stringsAsFactors = FALSE
  )
}))
write.csv(G_condition_summary,
          "results/G_condition_diagnostic.csv",
          row.names = FALSE)

save(G_condition, G_min_singular, G_max_singular, max_abs_transition,
     local_rho_array, source_diag,
     file = "results/stability_source_diagnostics.rda")

# Safe posterior summaries; stable-only arrays can contain all-NA slices.
safe_median <- function(z) {
  z <- z[is.finite(z)]
  if (!length(z)) return(NA_real_)
  median(z)
}
safe_quantile <- function(z, p) {
  z <- z[is.finite(z)]
  if (!length(z)) return(NA_real_)
  as.numeric(quantile(z, probs = p, names = FALSE))
}
summarize_irf <- function(arr) {
  list(
    med = apply(arr, c(1, 2, 3), safe_median),
    lo16 = apply(arr, c(1, 2, 3), safe_quantile, p = 0.16),
    hi84 = apply(arr, c(1, 2, 3), safe_quantile, p = 0.84),
    lo05 = apply(arr, c(1, 2, 3), safe_quantile, p = 0.05),
    hi95 = apply(arr, c(1, 2, 3), safe_quantile, p = 0.95)
  )
}

sum_all <- summarize_irf(IRF_post)
sum_stable <- summarize_irf(IRF_post_stable)

# Six comparison dates:
#   GPR-focused: 2003Q1, 2014Q3, 2022Q1, 2023Q4
#   crisis-state comparisons: 2008Q3, 2020Q1
plot_dates <- c("2003Q1", "2008Q3", "2014Q3", "2020Q1", "2022Q1", "2023Q4")
plot_dates <- plot_dates[plot_dates %in% dimnames(IRF_post)[[1]]]
vars <- setdiff(colnames(xglobal), "US_gpr")


# 4) MCMC convergence check on TVP coefficient chains at the six reported dates.
# A single-chain Geweke diagnostic is only treated as interpretable here when
# at least 100 retained posterior draws are available.
if (n_irf_draws >= 100L) {
  geweke_rows <- list()
  rr <- 0L

  for (i in seq_along(A.list)) {
    arr <- A.list[[i]]
    reg_names <- dimnames(arr)[[2]]
    eq_names <- dimnames(arr)[[3]]
    if (is.null(reg_names)) reg_names <- paste0("reg_", seq_len(dim(arr)[2]))
    if (is.null(eq_names)) eq_names <- paste0(cN[[i]], "_eq_", seq_len(dim(arr)[3]))

    for (d in plot_dates) {
      tt <- match(d, irf_dates)
      for (kk in seq_len(dim(arr)[2])) {
        for (mm in seq_len(dim(arr)[3])) {
          chain <- as.numeric(arr[tt, kk, mm, ])
          if (!all(is.finite(chain)) || sd(chain) <= .Machine$double.eps) next

          gz <- tryCatch(
            unname(coda::geweke.diag(coda::mcmc(chain),
                                     frac1 = 0.1, frac2 = 0.5)$z),
            error = function(e) NA_real_
          )

          rr <- rr + 1L
          geweke_rows[[rr]] <- data.frame(
            country = cN[[i]],
            date = d,
            regressor = reg_names[[kk]],
            equation = eq_names[[mm]],
            z = gz,
            abs_z = abs(gz),
            reject_5pct = is.finite(gz) && abs(gz) > 1.96,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  mcmc_geweke <- if (length(geweke_rows)) do.call(rbind, geweke_rows) else data.frame()
  write.csv(mcmc_geweke,
            "results/mcmc_geweke_diagnostic.csv",
            row.names = FALSE)

  if (nrow(mcmc_geweke)) {
    mcmc_geweke_summary <- do.call(rbind, lapply(split(mcmc_geweke,
                                                       mcmc_geweke$country),
                                                 function(z) {
      data.frame(
        country = unique(z$country)[1],
        retained_draws = n_irf_draws,
        tested_chains = nrow(z),
        rejection_share_5pct = mean(z$reject_5pct, na.rm = TRUE),
        median_abs_z = median(z$abs_z, na.rm = TRUE),
        p95_abs_z = as.numeric(quantile(z$abs_z, 0.95, na.rm = TRUE)),
        max_abs_z = max(z$abs_z, na.rm = TRUE),
        status = "interpretable_single_chain_geweke",
        stringsAsFactors = FALSE
      )
    }))
  } else {
    mcmc_geweke_summary <- data.frame(
      country = cN,
      retained_draws = n_irf_draws,
      tested_chains = 0L,
      rejection_share_5pct = NA_real_,
      median_abs_z = NA_real_,
      p95_abs_z = NA_real_,
      max_abs_z = NA_real_,
      status = "no_valid_chains",
      stringsAsFactors = FALSE
    )
  }
} else {
  mcmc_geweke_summary <- data.frame(
    country = cN,
    retained_draws = n_irf_draws,
    tested_chains = NA_integer_,
    rejection_share_5pct = NA_real_,
    median_abs_z = NA_real_,
    p95_abs_z = NA_real_,
    max_abs_z = NA_real_,
    status = "insufficient_retained_draws_for_formal_interpretation",
    stringsAsFactors = FALSE
  )
}
write.csv(mcmc_geweke_summary,
          "results/mcmc_geweke_summary.csv",
          row.names = FALSE)

make_irf_summary <- function(S, stable = FALSE) {
  do.call(rbind, lapply(plot_dates, function(d) {
    ti <- match(d, dimnames(IRF_post)[[1]])
    nstable <- sum(stable_mask[ti, ])
    do.call(rbind, lapply(vars, function(v) data.frame(
      date = d,
      variable = v,
      horizon = 0:nhor,
      median = S$med[ti, v, ],
      low68 = S$lo16[ti, v, ],
      high68 = S$hi84[ti, v, ],
      low90 = S$lo05[ti, v, ],
      high90 = S$hi95[ti, v, ],
      stable_only = stable,
      draws_used = if (stable) nstable else n_irf_draws,
      stringsAsFactors = FALSE
    )))
  }))
}

dd <- make_irf_summary(sum_all, stable = FALSE)
dd_stable <- make_irf_summary(sum_stable, stable = TRUE)
write.csv(dd, "results/gpr_structural_irf_summary.csv", row.names = FALSE)
write.csv(dd_stable, "results/gpr_structural_irf_summary_stable_only.csv", row.names = FALSE)

# Validate shock normalization. Under the baseline each posterior draw has exactly
# log(1 + shock_pct/100) impact on log GPR at h=0.
gpr0 <- IRF_post[, "US_gpr", "0", , drop = FALSE]
validation <- data.frame(
  date = dimnames(IRF_post)[[1]],
  target_log_jump = log1p(shock_pct / 100),
  median_h0 = apply(gpr0, 1, median, na.rm = TRUE),
  min_h0 = apply(gpr0, 1, min, na.rm = TRUE),
  max_h0 = apply(gpr0, 1, max, na.rm = TRUE)
)
write.csv(validation, "results/gpr_shock_validation.csv", row.names = FALSE)

# GDP cumulative IRF: y is 100*Delta log(real GDP), so the cumulative response
# approximates the response of 100*log(real GDP level) relative to baseline.
yvars <- grep("_y$", colnames(xglobal), value = TRUE)
yidx <- match(yvars, colnames(xglobal))

make_cumulative_y <- function(arr) {
  out <- arr[, yidx, , , drop = FALSE]
  for (tt in seq_len(dim(out)[1])) {
    for (vv in seq_len(dim(out)[2])) {
      for (ddraw in seq_len(dim(out)[4])) {
        out[tt, vv, , ddraw] <- cumsum(out[tt, vv, , ddraw])
      }
    }
  }
  out
}

cumY <- make_cumulative_y(IRF_post)
cumY_stable <- make_cumulative_y(IRF_post_stable)
cum_all <- summarize_irf(cumY)
cum_stable <- summarize_irf(cumY_stable)
y_dim_names <- dimnames(IRF_post)[[2]][yidx]

make_cum_summary <- function(S, stable = FALSE) {
  do.call(rbind, lapply(plot_dates, function(d) {
    ti <- match(d, dimnames(IRF_post)[[1]])
    nstable <- sum(stable_mask[ti, ])
    do.call(rbind, lapply(seq_along(y_dim_names), function(j) data.frame(
      date = d,
      variable = y_dim_names[j],
      horizon = 0:nhor,
      median = S$med[ti, j, ],
      low68 = S$lo16[ti, j, ],
      high68 = S$hi84[ti, j, ],
      low90 = S$lo05[ti, j, ],
      high90 = S$hi95[ti, j, ],
      stable_only = stable,
      draws_used = if (stable) nstable else n_irf_draws,
      stringsAsFactors = FALSE
    )))
  }))
}

cum_summary <- make_cum_summary(cum_all, stable = FALSE)
cum_summary_stable <- make_cum_summary(cum_stable, stable = TRUE)
write.csv(cum_summary, "results/gpr_structural_cumulative_gdp.csv", row.names = FALSE)
write.csv(cum_summary_stable,
          "results/gpr_structural_cumulative_gdp_stable_only.csv",
          row.names = FALSE)

# Compact GDP sign diagnostics for both all-draw and stable-only summaries.
diag_h <- intersect(c(0, 1, 4, 8, 12), 0:nhor)

make_sign_diag <- function(z) {
  out <- z[grepl("_y$", z$variable) & z$horizon %in% diag_h, ]
  out$sign <- ifelse(out$median > 0, "positive",
                     ifelse(out$median < 0, "negative", "zero"))
  out$credible68_excludes_zero <- with(out, low68 > 0 | high68 < 0)
  out
}
sign_diag <- make_sign_diag(dd)
sign_diag_stable <- make_sign_diag(dd_stable)
write.csv(sign_diag, "results/gdp_sign_diagnostic.csv", row.names = FALSE)
write.csv(sign_diag_stable,
          "results/gdp_sign_diagnostic_stable_only.csv",
          row.names = FALSE)

# Stability summary for the six comparison dates.
stability_selected <- stability_summary[stability_summary$date %in% plot_dates, ]
write.csv(stability_selected,
          "results/stability_summary_selected_dates.csv",
          row.names = FALSE)


source_selected <- source_diag[source_diag$date %in% plot_dates, ]
source_selected_summary <- do.call(rbind, lapply(split(source_selected, source_selected$date),
                                                 function(z) {
  top_source <- names(sort(table(z$likely_source), decreasing = TRUE))[1]
  data.frame(
    date = unique(z$date)[1],
    draws = nrow(z),
    global_unstable_share = mean(z$global_unstable),
    any_local_unstable_share = mean(z$any_local_unstable),
    G_ill_conditioned_share = mean(z$G_ill_conditioned),
    median_global_rho = median(z$global_rho),
    median_max_local_rho = median(z$max_local_rho),
    median_G_condition = median(z$G_condition_number[is.finite(z$G_condition_number)]),
    dominant_diagnostic_source = top_source,
    stringsAsFactors = FALSE
  )
}))
write.csv(source_selected_summary,
          "results/stability_source_selected_dates.csv",
          row.names = FALSE)

# All-draw figures (same interpretation as previous runs).
for (d in plot_dates) {
  for (vv in c("y", "dp", "r", "de", "deq")) {
    z <- dd[dd$date == d & grepl(paste0("_", vv, "$"), dd$variable), ]
    if (!nrow(z)) next
    p <- ggplot(z, aes(horizon, median)) +
      geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
      geom_ribbon(aes(ymin = low68, ymax = high68), alpha = 0.16) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~variable, scales = "free_y", ncol = 2) +
      theme_minimal(base_size = 10) +
      labs(
        title = paste0(d, ": responses to a +", shock_pct,
                       "% structural global GPR shock (all draws, ", vv, ")"),
        x = "Quarters after shock", y = "Response"
      )
    ggsave(sprintf("results/GPR_STRUCT_%s_%s.png", d, vv),
           p, width = 10, height = 13, dpi = 220)
  }

  zc <- cum_summary[cum_summary$date == d, ]
  p2 <- ggplot(zc, aes(horizon, median)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_ribbon(aes(ymin = low68, ymax = high68), alpha = 0.16) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~variable, scales = "free_y", ncol = 2) +
    theme_minimal(base_size = 10) +
    labs(
      title = paste0(d, ": cumulative GDP response to a +", shock_pct,
                     "% structural global GPR shock (all draws)"),
      x = "Quarters after shock", y = "Cumulative response"
    )
  ggsave(sprintf("results/GPR_STRUCT_%s_CUM_GDP.png", d),
         p2, width = 10, height = 13, dpi = 220)
}

# Stable-only figures: only rho < 1 draws are used.
for (d in plot_dates) {
  ti <- match(d, irf_dates)
  nstable <- sum(stable_mask[ti, ])
  if (nstable < 5L) {
    warning("Skipping stable-only figures for ", d,
            ": only ", nstable, " stable posterior draws.")
    next
  }

  for (vv in c("y", "dp", "r", "de", "deq")) {
    z <- dd_stable[
      dd_stable$date == d & grepl(paste0("_", vv, "$"), dd_stable$variable), ]
    if (!nrow(z)) next
    p <- ggplot(z, aes(horizon, median)) +
      geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
      geom_ribbon(aes(ymin = low68, ymax = high68), alpha = 0.16) +
      geom_line(linewidth = 0.7) +
      facet_wrap(~variable, scales = "free_y", ncol = 2) +
      theme_minimal(base_size = 10) +
      labs(
        title = paste0(d, ": responses to a +", shock_pct,
                       "% structural global GPR shock (stable draws only; n=",
                       nstable, "; ", vv, ")"),
        x = "Quarters after shock", y = "Response"
      )
    ggsave(sprintf("results/GPR_STRUCT_STABLE_%s_%s.png", d, vv),
           p, width = 10, height = 13, dpi = 220)
  }

  zc <- cum_summary_stable[cum_summary_stable$date == d, ]
  p2 <- ggplot(zc, aes(horizon, median)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_ribbon(aes(ymin = low68, ymax = high68), alpha = 0.16) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~variable, scales = "free_y", ncol = 2) +
    theme_minimal(base_size = 10) +
    labs(
      title = paste0(d, ": cumulative GDP response to a +", shock_pct,
                     "% structural global GPR shock (stable draws only; n=",
                     nstable, ")"),
      x = "Quarters after shock", y = "Cumulative response"
    )
  ggsave(sprintf("results/GPR_STRUCT_STABLE_%s_CUM_GDP.png", d),
         p2, width = 10, height = 13, dpi = 220)
}

overall_stable_share <- mean(stable_mask)
overall_near_share <- mean(near_unit_mask)
overall_unstable_share <- mean(unstable_mask)

selected_stability_lines <- if (nrow(stability_selected)) {
  paste0(
    stability_selected$date, ": stable=",
    sprintf("%.1f%%", 100 * stability_selected$stable_share),
    ", near-unit=", sprintf("%.1f%%", 100 * stability_selected$near_unit_share),
    ", p95 rho=", sprintf("%.4f", stability_selected$p95_rho)
  )
} else {
  character(0)
}

run_summary <- c(
  "Structural global GPR TVP-GVAR run",
  paste0("US recursive order: ", paste(expected_us, collapse = " -> ")),
  paste0("Shock normalization: +", shock_pct, "% global GPR on impact"),
  paste0("Non-US GPR lag specification: ", gpr_lag_mode),
  paste0("Domestic lag order p=", lag_order,
         "; foreign lag order q=1"),
  paste0("saves=", saves, "; burns=", burns, "; thin=", thin),
  paste0("retained IRF draws=", n_irf_draws),
  paste0("IRF horizon=", nhor, " quarters"),
  paste0("IRF dates: ", paste(plot_dates, collapse = ", ")),
  paste0("Stability rule: spectral radius < 1; near-unit warning threshold=",
         near_unit_threshold),
  paste0("Overall stable share=", sprintf("%.1f%%", 100 * overall_stable_share),
         "; near-unit share=", sprintf("%.1f%%", 100 * overall_near_share),
         "; unstable share=", sprintf("%.1f%%", 100 * overall_unstable_share)),
  paste0("Dominant instability diagnostic source: ",
         source_summary$likely_source[[1]], " (",
         sprintf("%.1f%%", 100 * source_summary$share[[1]]), " of date-draw cases)"),
  paste0("Residual Ljung-Box rejection share (5%): ",
         sprintf("%.1f%%", 100 * mean(residual_serial_diag$reject_5pct, na.rm = TRUE))),
  paste0("US local stable share: ",
         sprintf("%.1f%%", 100 * us_coefficient_summary$US_stable_share)),
  paste0("Median US short-rate AR sum: ",
         sprintf("%.4f", us_coefficient_summary$median_US_r_ar_sum)),
  paste0("Share of unstable US cases stabilized when the US_r row is removed (diagnostic only): ",
         ifelse(
           is.finite(us_coefficient_summary$unstable_share_fixed_by_US_r_row),
           sprintf("%.1f%%",
                   100 * us_coefficient_summary$unstable_share_fixed_by_US_r_row),
           "NA"
         )),
  paste0("Stationary feedback-test rejection share (5%): ",
         sprintf("%.1f%%", 100 * mean(feedback_diag$reject_no_feedback_5pct, na.rm = TRUE))),
  paste0("MCMC Geweke status: ",
         paste(unique(mcmc_geweke_summary$status), collapse = ", ")),
  "Selected-date stability:",
  selected_stability_lines,
  "All-draw and stable-only IRF summaries are both saved; near-unit stable draws are flagged but not discarded.",
  "Interpret y as GDP-growth response; use cumulative GDP IRF for the implied GDP-level path."
)
writeLines(run_summary, "results/run_structural_summary.txt")
cat(paste(run_summary, collapse = "\n"), "\n")
