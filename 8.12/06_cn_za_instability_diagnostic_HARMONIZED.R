#!/usr/bin/env Rscript

# ============================================================
# Harmonized local-stability diagnostic for GPR + Brent GVAR
#
# This is the successor to 06_cn_za_instability_diagnostic.R.
# It deliberately reproduces the local VARX design used by
# 05_country_specific_lag_selection.R:
#   (1) domestic lags L1..Lp;
#   (2) foreign-star variables at t and t-1;
#   (3) current and L1 GPR and Brent controls;
#   (4) a common p=2-trimmed estimation window for p=1 and p=2.
#
# Its companion roots are therefore comparable to the lag-selection step.
# ============================================================

source("8.12/06_cn_za_instability_diagnostic.R", local = FALSE)

# Do not overwrite the previous diagnostic artifact.
OUT_DIR <- "8.12/cn_za_instability_diagnostic_harmonized"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 05 starts every candidate model after the largest candidate domestic lag.
# Retaining this rule makes p=1 and p=2 diagnostics sample-comparable.
P_CANDIDATES <- P_SET

fit_local_varx <- function(panel, cc, p) {
  i <- match(cc, COUNTRIES)
  if (is.na(i)) stopf("Country not found: %s", cc)

  Y  <- panel$X[, i, , drop = FALSE][, 1, ]
  Xs <- panel$Xstar[, i, , drop = FALSE][, 1, ]
  Y  <- as.matrix(Y)
  Xs <- as.matrix(Xs)
  colnames(Y) <- VARS
  colnames(Xs) <- VARS

  Tn <- nrow(Y)
  first_row <- max(P_CANDIDATES) + 1L
  rows <- seq.int(first_row, Tn)

  D <- data.frame(const = rep(1, length(rows)))

  # Domestic lags L1..Lp.
  for (L in seq_len(p)) {
    for (v in seq_along(VARS)) {
      D[[paste0(VARS[v], "_L", L)]] <- lag_vec(Y[, v], L)[rows]
    }
  }

  # Exact 05 foreign-star convention: contemporaneous x*_t plus x*_{t-1}.
  for (v in seq_along(VARS)) {
    D[[paste0(VARS[v], "_star_0")]] <- Xs[rows, v]
    D[[paste0(VARS[v], "_star_L1")]] <- lag_vec(Xs[, v], 1L)[rows]
  }

  # Exact 05 global convention: current and L1 controls.
  D$gpr_0  <- panel$gpr[rows]
  D$gpr_L1 <- lag_vec(panel$gpr, 1L)[rows]
  D$oil_0  <- panel$oil[rows]
  D$oil_L1 <- lag_vec(panel$oil, 1L)[rows]

  YY <- Y[rows, , drop = FALSE]
  ok <- complete.cases(D) & complete.cases(YY)
  D <- D[ok, , drop = FALSE]
  YY <- YY[ok, , drop = FALSE]
  qid <- panel$qid[rows][ok]

  Xmat <- as.matrix(D)
  B <- matrix(NA_real_, nrow = ncol(Xmat), ncol = ncol(YY),
              dimnames = list(colnames(Xmat), VARS))
  res <- matrix(NA_real_, nrow = nrow(YY), ncol = ncol(YY),
                dimnames = list(NULL, VARS))
  for (j in seq_along(VARS)) {
    fit <- lm.fit(Xmat, YY[, j])
    B[, j] <- fit$coefficients
    res[, j] <- fit$residuals
  }

  A <- vector("list", p)
  for (L in seq_len(p)) {
    Am <- matrix(0, length(VARS), length(VARS), dimnames = list(VARS, VARS))
    for (eq in VARS) for (v in VARS) {
      Am[eq, v] <- B[paste0(v, "_L", L), eq]
    }
    A[[L]] <- Am
  }

  list(B = B, residuals = res, A = A, n = nrow(YY), Y = YY, X = Xmat,
       rows = rows[ok], qid = qid)
}

# Run the original main block a second time, now with the harmonized function
# and a separate output directory.
stability_rows <- list(); eig_rows <- list(); coef_rows <- list()
eq_rows <- list(); univ_rows <- list(); cf_rows <- list(); design_rows <- list()

for (cc in TARGETS) for (p in P_SET) {
  msg("Fitting harmonized %s p=%d ...", cc, p)
  fit <- fit_local_varx(panel, cc, p)
  C <- companion_matrix(fit$A)
  ev <- eigen(C, only.values = TRUE)$values
  ord <- order(Mod(ev), decreasing = TRUE)
  stability_rows[[length(stability_rows) + 1L]] <- data.frame(
    country = cc, p = p, nobs = fit$n,
    first_quarter = quarter_label(min(fit$qid)),
    last_quarter = quarter_label(max(fit$qid)),
    spectral_radius = max(Mod(ev)), stable = max(Mod(ev)) < 1,
    borderline = max(Mod(ev)) >= 0.98 & max(Mod(ev)) < 1)
  eig_rows[[length(eig_rows) + 1L]] <- data.frame(
    country = cc, p = p, rank = seq_along(ord), real = Re(ev[ord]),
    imaginary = Im(ev[ord]), modulus = Mod(ev[ord]))
  for (L in seq_len(p)) {
    A_l <- fit$A[[L]]
    for (v in VARS) coef_rows[[length(coef_rows) + 1L]] <- data.frame(
      country = cc, p = p, lag = L, lagged_variable = v,
      column_l1 = sum(abs(A_l[, v])),
      column_l2 = sqrt(sum(A_l[, v]^2)),
      column_maxabs = max(abs(A_l[, v])))
    for (eq in VARS) eq_rows[[length(eq_rows) + 1L]] <- data.frame(
      country = cc, p = p, lag = L, equation = eq,
      own_lag_coef = A_l[eq, eq], row_l1 = sum(abs(A_l[eq, ])),
      row_l2 = sqrt(sum(A_l[eq, ]^2)), row_maxabs = max(abs(A_l[eq, ])))
  }
  cf_rows[[length(cf_rows) + 1L]] <- zero_block_counterfactuals(fit, cc, p)
  design_rows[[length(design_rows) + 1L]] <- data.frame(
    country = cc, p = p, qid = fit$qid, quarter = quarter_label(fit$qid),
    n_regressors = ncol(fit$X), stringsAsFactors = FALSE)
}

# These are data diagnostics, so calculate each country once from the common
# underlying panel rather than duplicating them for each candidate p.
for (cc in TARGETS) {
  i <- match(cc, COUNTRIES)
  for (v in VARS) {
    x <- panel$X[, i, v]
    ar <- safe_ar1(x)
    univ_rows[[length(univ_rows) + 1L]] <- data.frame(
      country = cc, variable = v, ar1 = ar["phi"], n = ar["n"],
      mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE),
      min = min(x, na.rm = TRUE), max = max(x, na.rm = TRUE))
  }
}

stability <- do.call(rbind, stability_rows); eigs <- do.call(rbind, eig_rows)
coefnorm <- do.call(rbind, coef_rows); eqsum <- do.call(rbind, eq_rows)
univ <- do.call(rbind, univ_rows); cf <- do.call(rbind, cf_rows)
sus <- cf[cf$block_type == "variable_column", ]
sus <- sus[order(sus$country, sus$p, -sus$improvement), ]
sus$rank_within_model <- ave(-sus$improvement, interaction(sus$country, sus$p),
                             FUN = function(x) rank(x, ties.method = "first"))
sus <- sus[sus$rank_within_model <= 5, ]

write.csv(stability, file.path(OUT_DIR, "01_local_stability_summary.csv"), row.names = FALSE)
write.csv(eigs, file.path(OUT_DIR, "02_companion_eigenvalues.csv"), row.names = FALSE)
write.csv(coefnorm, file.path(OUT_DIR, "03_domestic_lag_coefficient_norms.csv"), row.names = FALSE)
write.csv(eqsum, file.path(OUT_DIR, "04_equation_persistence_summary.csv"), row.names = FALSE)
write.csv(univ, file.path(OUT_DIR, "05_variable_univariate_persistence.csv"), row.names = FALSE)
write.csv(cf, file.path(OUT_DIR, "06_counterfactual_zero_block_stability.csv"), row.names = FALSE)
write.csv(sus, file.path(OUT_DIR, "07_suspect_ranking.csv"), row.names = FALSE)
write.csv(do.call(rbind, design_rows), file.path(OUT_DIR, "08_estimation_sample_audit.csv"), row.names = FALSE)

writeLines(c(
  "HARMONIZED CN + ZA LOCAL INSTABILITY DIAGNOSTIC",
  "", "This diagnostic exactly follows 05_country_specific_lag_selection.R.",
  "Foreign-star regressors: current x*_t and lagged x*_{t-1}.",
  "Global regressors: current and L1 GPR and Brent.",
  sprintf("All p models start at max(P_CANDIDATES)+1 = %d.", max(P_CANDIDATES) + 1L),
  "Use these results, rather than the prior 06 outputs, to assess local stability.",
  "08_estimation_sample_audit.csv records the exact retained quarters."),
  file.path(OUT_DIR, "README_cn_za_instability_diagnostic_harmonized.txt"))

msg("=== Harmonized stability summary ===")
print(stability, row.names = FALSE)
