#!/usr/bin/env Rscript

# =============================================================================
# Reverse-order identification robustness:
#   Main:    GPR -> VIX
#   Reverse: VIX -> GPR
#
# IMPORTANT
# - Reuses the SAME saved TVP-GVAR posterior.
# - Does NOT rerun MCMC.
# - Does NOT change data, weights, lag order, priors, or coefficient draws.
# - Only changes the recursive contemporaneous ordering inside the dominant
#   2x2 block [GL_gpr, GL_vix].
# - Primary inference uses stable-only draws (global spectral radius < 1).
# - Only evaluates six crisis dates and y / r / deq responses.
#
# Expected repository path:
#   8.12/09_reverse_order_irf_6crises_from_saved_posterior.R
# =============================================================================

rm(list = ls())

posterior_file <- Sys.getenv(
  "TVPGVAR_DOMINANT_POSTERIOR",
  "prior_artifact/results/predDens_dominant_gpr_vix.rda"
)
out_dir <- Sys.getenv("TVPGVAR_REVERSE_RESULTS", "results_reverse_order")
nhor <- as.integer(Sys.getenv("TVPGVAR_HORIZON", "12"))
shock_pct <- as.numeric(Sys.getenv("TVPGVAR_GPR_SHOCK_PCT", "10"))

crisis_dates_requested <- c("2003Q1", "2008Q3", "2014Q3", "2020Q1", "2022Q1", "2023Q4")
responses_requested <- c("y", "r", "deq")
comparison_horizons_requested <- c(0L, 1L, 4L, 8L, 12L)

if (!file.exists(posterior_file)) {
  stop("Missing saved posterior: ", posterior_file)
}
if (!is.finite(nhor) || nhor < 1L) stop("Invalid TVPGVAR_HORIZON.")
if (!is.finite(shock_pct) || shock_pct <= 0) stop("Invalid TVPGVAR_GPR_SHOCK_PCT.")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Reuse the current p=2,q=1 reconstruction/stability machinery.
source("8.12/dominant_gpr_vix_irf.R")

# Keep a reference to the original implementation for provenance/debugging.
dominant_gpr_struct_irf_main_original <- dominant_gpr_struct_irf

# -----------------------------------------------------------------------------
# Identification override
# -----------------------------------------------------------------------------
# get_dominant_gpr_irf_t() calls dominant_gpr_struct_irf() by name.
# We replace ONLY that structural-shock constructor. Everything else in the
# existing GVAR reconstruction remains unchanged.
#
# Correct reverse ordering:
#   1) Reorder S from [GPR,VIX] to [VIX,GPR]
#   2) Cholesky in the reordered system
#   3) Take the GPR structural shock (column 2 in [VIX,GPR])
#   4) Map the shock vector back to [GPR,VIX]
#
# Directly taking column 2 of chol(S) in the original [GPR,VIX] ordering is NOT
# equivalent to a VIX -> GPR recursive ordering.
# -----------------------------------------------------------------------------

.tvpgvar_identification_order <- "GPR_VIX"

dominant_gpr_struct_irf <- function(
    G, F, sig, x, units,
    horizon = 12,
    shock_var = "GL_gpr",
    shock_pct = 10) {

  K <- nrow(x)
  if (!all(dim(G) == c(K, K))) stop("G must be KxK.")

  gl_i <- match("GL", units)
  if (is.na(gl_i)) stop("Dominant unit GL not found.")

  S <- stabilize_cov(sig[[gl_i]])
  if (!all(dim(S) == c(2L, 2L))) {
    stop("Dominant covariance must be 2x2.")
  }

  base_names <- c("GL_gpr", "GL_vix")
  ordering <- .tvpgvar_identification_order

  if (identical(ordering, "GPR_VIX")) {
    ordered_names <- c("GL_gpr", "GL_vix")
  } else if (identical(ordering, "VIX_GPR")) {
    ordered_names <- c("GL_vix", "GL_gpr")
  } else {
    stop("Unknown identification ordering: ", ordering)
  }

  # Reorder covariance matrix to the requested recursive ordering.
  ord_idx <- match(ordered_names, base_names)
  if (anyNA(ord_idx)) stop("Internal dominant-variable ordering error.")

  S_ord <- stabilize_cov(S[ord_idx, ord_idx, drop = FALSE])
  L_ord <- t(chol(S_ord))

  # Structural GPR shock is whichever column corresponds to GL_gpr in the
  # requested recursive order.
  gpr_col <- match("GL_gpr", ordered_names)
  shock_ord <- L_ord[, gpr_col]

  # Map innovation vector back into the model's stored [GL_gpr, GL_vix] order.
  shock_gl <- numeric(2L)
  shock_gl[ord_idx] <- shock_ord
  names(shock_gl) <- base_names

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

  # Same normalization in both orderings: +shock_pct% GPR impact at h=0.
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

  out <- sapply(
    seq_len(horizon + 1L),
    function(h) phi[, , h] %*% impact
  )
  rownames(out) <- rownames(x)
  colnames(out) <- 0:horizon
  out
}

# -----------------------------------------------------------------------------
# Load the exact same estimated posterior
# -----------------------------------------------------------------------------
load(posterior_file)

if (!exists("predDens") || !exists("Data.setup")) {
  stop("Saved posterior must contain predDens and Data.setup.")
}

xglobal <- Data.setup$bigx
globalG <- lapply(predDens, `[[`, "W")
cN <- Data.setup$countries
country_units <- Data.setup$country_units

if (is.null(xglobal) || is.null(cN) || is.null(country_units)) {
  stop("Data.setup is incomplete.")
}

A.list <- lapply(predDens, `[[`, "ALPHA")
S.list <- lapply(predDens, `[[`, "SIGMApost")

# Infer domestic lag p from the posterior labels; do not impose a new p.
rn0 <- dimnames(predDens[[1]]$ALPHA)[[2]]
lag_tags0 <- grep("^Ylag[0-9]+$", unique(rn0), value = TRUE)
lag_nums0 <- suppressWarnings(as.integer(sub("^Ylag", "", lag_tags0)))
if (!length(lag_nums0) || anyNA(lag_nums0)) {
  stop("Cannot infer domestic lag order from saved posterior.")
}
lag_order <- max(lag_nums0)

n_irf_draws <- dim(A.list[[1]])[4]
if (any(vapply(A.list, function(z) dim(z)[4], integer(1)) != n_irf_draws)) {
  stop("Inconsistent retained posterior draw counts.")
}

Tirf <- nrow(xglobal) - lag_order
irf_dates <- Data.setup$quarters[-seq_len(lag_order)]

if (length(irf_dates) != Tirf) stop("IRF date alignment mismatch.")
if (any(vapply(A.list, function(z) dim(z)[1], integer(1)) != Tirf)) {
  stop("Coefficient paths are not aligned across units.")
}

selected_dates <- crisis_dates_requested[crisis_dates_requested %in% irf_dates]
missing_dates <- setdiff(crisis_dates_requested, selected_dates)
if (length(missing_dates)) {
  warning("Requested crisis dates not available: ", paste(missing_dates, collapse = ", "))
}
if (!length(selected_dates)) stop("None of the requested crisis dates are available.")

selected_tt <- match(selected_dates, irf_dates)

# Only y/r/deq country variables.
all_vars <- colnames(xglobal)
response_suffix <- sub("^.*_", "", all_vars)
target_vars <- all_vars[
  response_suffix %in% responses_requested &
    !grepl("^GL_", all_vars)
]
if (!length(target_vars)) {
  stop("No country y/r/deq variables found in xglobal.")
}

comparison_horizons <- comparison_horizons_requested[
  comparison_horizons_requested <= nhor
]

cat("============================================================\n")
cat(" REVERSE-ORDER IDENTIFICATION ROBUSTNESS — NO MCMC\n")
cat(" Posterior: ", posterior_file, "\n", sep = "")
cat(" Inferred p: ", lag_order, "\n", sep = "")
cat(" Foreign/global q: unchanged from saved specification\n")
cat(" Retained posterior draws: ", n_irf_draws, "\n", sep = "")
cat(" Dates: ", paste(selected_dates, collapse = ", "), "\n", sep = "")
cat(" Responses: y, r, deq\n")
cat(" Horizon: 0-", nhor, "\n", sep = "")
cat(" GPR normalization: +", shock_pct, "% at h=0\n", sep = "")
cat(" Main ordering: GPR -> VIX\n")
cat(" Reverse ordering: VIX -> GPR\n")
cat("============================================================\n\n")

# -----------------------------------------------------------------------------
# Reconstruct only the 4 requested dates from the same posterior.
# Stability is identification-order invariant because G and F are unchanged.
# -----------------------------------------------------------------------------

nd <- length(selected_dates)
K <- ncol(xglobal)

IRF_main <- array(
  NA_real_,
  dim = c(nd, K, nhor + 1L, n_irf_draws),
  dimnames = list(selected_dates, all_vars, 0:nhor, NULL)
)
IRF_reverse <- IRF_main

stability_rho <- matrix(
  NA_real_,
  nrow = nd,
  ncol = n_irf_draws,
  dimnames = list(selected_dates, NULL)
)

rho_ordering_diff <- matrix(
  NA_real_,
  nrow = nd,
  ncol = n_irf_draws,
  dimnames = list(selected_dates, NULL)
)

for (dd in seq_len(n_irf_draws)) {

  A.i <- lapply(A.list, function(z) z[, , , dd])
  S.i <- lapply(S.list, function(z) z[, , , dd])

  for (jj in seq_along(selected_tt)) {

    tt <- selected_tt[jj]

    .tvpgvar_identification_order <- "GPR_VIX"
    ans_main <- get_dominant_gpr_irf_t(
      tt = tt,
      draw_i = A.i,
      Sig_draw_i = S.i,
      x = t(xglobal),
      globalG = globalG,
      units = cN,
      horizon = nhor,
      shock_pct = shock_pct
    )

    .tvpgvar_identification_order <- "VIX_GPR"
    ans_reverse <- get_dominant_gpr_irf_t(
      tt = tt,
      draw_i = A.i,
      Sig_draw_i = S.i,
      x = t(xglobal),
      globalG = globalG,
      units = cN,
      horizon = nhor,
      shock_pct = shock_pct
    )

    IRF_main[jj, , , dd] <- ans_main$IRF_post
    IRF_reverse[jj, , , dd] <- ans_reverse$IRF_post

    stability_rho[jj, dd] <- ans_main$max_eigen_modulus
    rho_ordering_diff[jj, dd] <- abs(
      ans_main$max_eigen_modulus - ans_reverse$max_eigen_modulus
    )
  }

  cat("Posterior draw ", dd, "/", n_irf_draws, " complete\n", sep = "")
}

# The ordering must not alter model stability.
max_rho_diff <- max(rho_ordering_diff, na.rm = TRUE)
if (!is.finite(max_rho_diff) || max_rho_diff > 1e-10) {
  stop(
    "Unexpected stability difference across orderings. max |rho_main-rho_reverse| = ",
    format(max_rho_diff, scientific = TRUE)
  )
}

stable_mask <- stability_rho < 1

# -----------------------------------------------------------------------------
# Stable-only posterior summaries
# -----------------------------------------------------------------------------

summarize_order <- function(IRF_array, ordering_label) {
  rows <- list()

  for (jj in seq_along(selected_dates)) {
    date <- selected_dates[jj]
    draws <- which(stable_mask[jj, ])

    if (!length(draws)) {
      warning("No stable draws at ", date, "; stable-only summary omitted.")
      next
    }

    for (v in target_vars) {
      vi <- match(v, all_vars)

      for (h in 0:nhor) {
        z <- IRF_array[jj, vi, h + 1L, draws]
        z <- z[is.finite(z)]
        if (!length(z)) next

        qs <- quantile(z, c(.05, .16, .50, .84, .95), names = FALSE)

        rows[[length(rows) + 1L]] <- data.frame(
          ordering = ordering_label,
          date = date,
          country = sub("_(y|r|deq)$", "", v),
          response = sub("^.*_", "", v),
          variable = v,
          horizon = h,
          median = qs[3],
          low68 = qs[2],
          high68 = qs[4],
          low90 = qs[1],
          high90 = qs[5],
          p_positive = mean(z > 0),
          p_negative = mean(z < 0),
          sig68 = (qs[2] > 0 || qs[4] < 0),
          sig90 = (qs[1] > 0 || qs[5] < 0),
          stable_draws_used = length(z),
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (!length(rows)) return(data.frame())
  do.call(rbind, rows)
}

summary_main <- summarize_order(IRF_main, "GPR_to_VIX")
summary_reverse <- summarize_order(IRF_reverse, "VIX_to_GPR")

write.csv(
  summary_main,
  file.path(out_dir, "01_GPR_to_VIX_irf_stable_only.csv"),
  row.names = FALSE
)
write.csv(
  summary_reverse,
  file.path(out_dir, "02_VIX_to_GPR_irf_stable_only.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Direct comparison at h = 0,1,4,8,12 (when available)
# -----------------------------------------------------------------------------
main_cmp <- summary_main[
  summary_main$horizon %in% comparison_horizons,
  ,
  drop = FALSE
]
rev_cmp <- summary_reverse[
  summary_reverse$horizon %in% comparison_horizons,
  ,
  drop = FALSE
]

keep_main <- c(
  "date", "country", "response", "variable", "horizon",
  "median", "low68", "high68", "low90", "high90",
  "p_positive", "p_negative", "sig68", "sig90", "stable_draws_used"
)
keep_rev <- keep_main

main_cmp <- main_cmp[, keep_main, drop = FALSE]
rev_cmp <- rev_cmp[, keep_rev, drop = FALSE]

names(main_cmp)[6:ncol(main_cmp)] <- paste0(
  names(main_cmp)[6:ncol(main_cmp)], "_main"
)
names(rev_cmp)[6:ncol(rev_cmp)] <- paste0(
  names(rev_cmp)[6:ncol(rev_cmp)], "_reverse"
)

comparison <- merge(
  main_cmp,
  rev_cmp,
  by = c("date", "country", "response", "variable", "horizon"),
  all = TRUE,
  sort = TRUE
)

sign_eps <- 1e-12
sign_class <- function(z) {
  ifelse(!is.finite(z), NA_integer_,
         ifelse(abs(z) <= sign_eps, 0L, ifelse(z > 0, 1L, -1L)))
}

comparison$sign_main <- sign_class(comparison$median_main)
comparison$sign_reverse <- sign_class(comparison$median_reverse)
comparison$same_sign <- (
  comparison$sign_main == comparison$sign_reverse
)
comparison$sign_flip <- (
  comparison$sign_main * comparison$sign_reverse == -1L
)
comparison$median_difference_reverse_minus_main <- (
  comparison$median_reverse - comparison$median_main
)

write.csv(
  comparison,
  file.path(out_dir, "03_ordering_comparison_selected_horizons.csv"),
  row.names = FALSE
)

# No arbitrary pass/fail threshold is imposed here.
# Report the share of countries with the same sign / flipped sign.
agreement_rows <- list()
group_keys <- unique(
  comparison[, c("date", "response", "horizon"), drop = FALSE]
)

for (ii in seq_len(nrow(group_keys))) {
  key <- group_keys[ii, ]
  z <- comparison[
    comparison$date == key$date &
      comparison$response == key$response &
      comparison$horizon == key$horizon,
    ,
    drop = FALSE
  ]

  valid <- is.finite(z$sign_main) & is.finite(z$sign_reverse)
  z <- z[valid, , drop = FALSE]
  if (!nrow(z)) next

  agreement_rows[[length(agreement_rows) + 1L]] <- data.frame(
    date = key$date,
    response = key$response,
    horizon = key$horizon,
    countries_compared = nrow(z),
    same_sign_countries = sum(z$same_sign),
    same_sign_share = mean(z$same_sign),
    sign_flip_countries = sum(z$sign_flip),
    sign_flip_share = mean(z$sign_flip),
    stringsAsFactors = FALSE
  )
}

agreement_summary <- if (length(agreement_rows)) {
  do.call(rbind, agreement_rows)
} else {
  data.frame()
}

write.csv(
  agreement_summary,
  file.path(out_dir, "04_sign_agreement_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Stability summary for the six crisis dates
# -----------------------------------------------------------------------------
stability_summary <- do.call(
  rbind,
  lapply(seq_along(selected_dates), function(jj) {
    r <- stability_rho[jj, ]
    data.frame(
      date = selected_dates[jj],
      total_draws = length(r),
      stable_draws = sum(r < 1, na.rm = TRUE),
      stable_share = mean(r < 1, na.rm = TRUE),
      median_rho = median(r, na.rm = TRUE),
      p95_rho = as.numeric(quantile(r, .95, na.rm = TRUE)),
      max_rho = max(r, na.rm = TRUE),
      max_ordering_rho_difference = max(rho_ordering_diff[jj, ], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

write.csv(
  stability_summary,
  file.path(out_dir, "05_stability_selected_dates.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Identification validation
# -----------------------------------------------------------------------------
target_log_jump <- log1p(shock_pct / 100)
gpr_i <- match("GL_gpr", all_vars)
vix_i <- match("GL_vix", all_vars)

if (is.na(gpr_i) || is.na(vix_i)) {
  stop("GL_gpr / GL_vix not found in xglobal colnames.")
}

validation_rows <- list()

for (jj in seq_along(selected_dates)) {
  stable_draws <- which(stable_mask[jj, ])
  if (!length(stable_draws)) next

  for (ordering in c("GPR_to_VIX", "VIX_to_GPR")) {
    X <- if (ordering == "GPR_to_VIX") IRF_main else IRF_reverse

    gpr_h0 <- X[jj, gpr_i, 1L, stable_draws]
    vix_h0 <- X[jj, vix_i, 1L, stable_draws]

    validation_rows[[length(validation_rows) + 1L]] <- data.frame(
      date = selected_dates[jj],
      ordering = ordering,
      target_gpr_h0 = target_log_jump,
      median_gpr_h0 = median(gpr_h0),
      max_abs_gpr_normalization_error = max(abs(gpr_h0 - target_log_jump)),
      median_vix_h0 = median(vix_h0),
      max_abs_vix_h0 = max(abs(vix_h0)),
      stable_draws_used = length(stable_draws),
      stringsAsFactors = FALSE
    )
  }
}

identification_validation <- do.call(rbind, validation_rows)

write.csv(
  identification_validation,
  file.path(out_dir, "06_identification_validation.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Compact plotting PDF: 18 pages = 6 crisis dates x 3 response variables.
# Each page shows countries as small multiples.
# Solid = GPR -> VIX; dashed = VIX -> GPR.
# Outer dotted/dot-dash lines = 90% intervals.
# -----------------------------------------------------------------------------

plot_file <- file.path(out_dir, "07_reverse_order_irf_compare_stable_only.pdf")
pdf(plot_file, width = 11, height = 8.5, onefile = TRUE)

for (date in selected_dates) {
  for (resp in responses_requested) {

    vars_resp <- target_vars[sub("^.*_", "", target_vars) == resp]
    if (!length(vars_resp)) next

    n_pan <- length(vars_resp)
    ncol_pan <- 4L
    nrow_pan <- ceiling(n_pan / ncol_pan)

    old_par <- par(
      mfrow = c(nrow_pan, ncol_pan),
      mar = c(2.6, 2.8, 2.0, 0.8),
      oma = c(2.2, 2.2, 3.2, 0.5)
    )

    for (v in vars_resp) {
      m <- summary_main[
        summary_main$date == date & summary_main$variable == v,
        ,
        drop = FALSE
      ]
      r <- summary_reverse[
        summary_reverse$date == date & summary_reverse$variable == v,
        ,
        drop = FALSE
      ]

      if (!nrow(m) || !nrow(r)) {
        plot.new()
        title(main = v)
        text(.5, .5, "No stable draws")
        next
      }

      yr <- range(
        c(m$low90, m$high90, r$low90, r$high90, 0),
        finite = TRUE
      )
      if (!all(is.finite(yr)) || diff(yr) == 0) yr <- c(-1, 1)

      plot(
        m$horizon, m$median,
        type = "l",
        lty = 1,
        lwd = 2,
        ylim = yr,
        xlab = "",
        ylab = "",
        main = sub("_(y|r|deq)$", "", v)
      )
      abline(h = 0, lty = 3)

      lines(m$horizon, m$low90, lty = 3)
      lines(m$horizon, m$high90, lty = 3)

      lines(r$horizon, r$median, lty = 2, lwd = 2)
      lines(r$horizon, r$low90, lty = 4)
      lines(r$horizon, r$high90, lty = 4)
    }

    # Fill unused panels.
    unused <- nrow_pan * ncol_pan - n_pan
    if (unused > 0) {
      for (ii in seq_len(unused)) plot.new()
    }

    mtext(
      paste0(
        date, " | response=", resp,
        " | stable-only | solid: GPR->VIX, dashed: VIX->GPR"
      ),
      side = 3, outer = TRUE, line = 1.0, cex = 1.0
    )
    mtext("Horizon (quarters)", side = 1, outer = TRUE, line = 0.6)
    mtext("IRF", side = 2, outer = TRUE, line = 0.6)

    par(old_par)
  }
}

dev.off()

# -----------------------------------------------------------------------------
# Save draw-level arrays for later diagnostics without another reconstruction.
# -----------------------------------------------------------------------------
save(
  IRF_main,
  IRF_reverse,
  stability_rho,
  stable_mask,
  rho_ordering_diff,
  selected_dates,
  target_vars,
  lag_order,
  nhor,
  shock_pct,
  file = file.path(out_dir, "reverse_order_irf_draw_arrays.rda")
)

# -----------------------------------------------------------------------------
# README / provenance
# -----------------------------------------------------------------------------
writeLines(
  c(
    "TVP-GVAR reverse-order identification robustness",
    "================================================",
    paste0("source_posterior=", posterior_file),
    paste0("inferred_domestic_p=", lag_order),
    "foreign_global_q=UNCHANGED_FROM_SAVED_MODEL",
    paste0("horizon=", nhor),
    paste0("shock_pct=", shock_pct),
    paste0("crisis_dates=", paste(selected_dates, collapse = ",")),
    "responses=y,r,deq",
    "main_ordering=GPR->VIX",
    "reverse_ordering=VIX->GPR",
    "normalization=same +GPR impact in both orderings",
    "stable_rule=global spectral radius < 1",
    "primary_inference=stable-only draws",
    "mcmc_rerun=FALSE",
    "data_changed=FALSE",
    "weights_changed=FALSE",
    "lag_order_changed=FALSE",
    "posterior_changed=FALSE",
    paste0("max_stability_difference_across_orderings=", format(max_rho_diff, scientific = TRUE)),
    "",
    "Interpretation:",
    "1) Same direction after reversal -> main sign is not driven by this recursive ordering.",
    "2) Same direction but different magnitude -> shared contemporaneous GPR/VIX component.",
    "3) Broad sign reversal -> avoid interpreting the baseline as a pure GPR shock; treat identification uncertainty as material."
  ),
  file.path(out_dir, "README_reverse_order_identification.txt")
)

cat("\nReverse-order identification robustness complete.\n")
cat("Output directory: ", out_dir, "\n", sep = "")
cat("Stable shares:\n")
print(stability_summary[, c("date", "stable_draws", "stable_share")])
cat("\nSign-agreement summary:\n")
print(agreement_summary)
