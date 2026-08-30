#!/usr/bin/env Rscript

# ============================================================
# Global GVAR stability diagnostic, harmonized specification
#
# Tests the stacked GVAR transition system implied by:
#   x_i,t = Phi_i(L)x_i,t-1 + Lambda_i,0 x*_i,t
#         + Lambda_i,1 x*_i,t-1 + controls + e_i,t
# where x*_i,t = sum_j w_ij x_j,t.
#
# The current x*_i,t term is moved to the left-hand side to form
# G0 x_t = G1 x_{t-1} + ...; global roots are then computed from
# the companion matrix of F_l = solve(G0, G_l).
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

# This source supplies the verified data readers, panel constructor, weights,
# and constants. It also runs the legacy local diagnostic; its output is not
# used below. The global design is rebuilt explicitly in this file.
source("8.12/06_cn_za_instability_diagnostic.R", local = FALSE)

OUT_DIR <- "8.12/global_gvar_stability_diagnostic"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Initial global checks: the same p in every country. This isolates whether
# the intended contemporaneous foreign-star GVAR architecture is globally
# stable. After lag selection is finalized, replace these with one named
# country-specific vector and run again.
P_GRIDS <- list(p1_all = setNames(rep(1L, length(COUNTRIES)), COUNTRIES),
                p2_all = setNames(rep(2L, length(COUNTRIES)), COUNTRIES))

fit_harmonized_local <- function(panel, cc, p) {
  i <- match(cc, COUNTRIES)
  Y <- as.matrix(panel$X[, i, , drop = FALSE][, 1, ])
  Xs <- as.matrix(panel$Xstar[, i, , drop = FALSE][, 1, ])
  colnames(Y) <- colnames(Xs) <- VARS

  # Match the 05 selection convention: all candidate p models use the same
  # p=2-trimmed start, which is 2000Q4 in the current panel.
  rows <- seq.int(max(c(1L, 2L)) + 1L, nrow(Y))
  D <- data.frame(const = rep(1, length(rows)))
  for (L in seq_len(p)) for (v in seq_along(VARS))
    D[[paste0(VARS[v], "_L", L)]] <- lag_vec(Y[, v], L)[rows]
  for (v in seq_along(VARS)) {
    D[[paste0(VARS[v], "_star_0")]] <- Xs[rows, v]
    D[[paste0(VARS[v], "_star_L1")]] <- lag_vec(Xs[, v], 1L)[rows]
  }
  D$gpr_0 <- panel$gpr[rows]
  D$gpr_L1 <- lag_vec(panel$gpr, 1L)[rows]
  D$oil_0 <- panel$oil[rows]
  D$oil_L1 <- lag_vec(panel$oil, 1L)[rows]

  YY <- Y[rows, , drop = FALSE]
  ok <- complete.cases(D) & complete.cases(YY)
  X <- as.matrix(D[ok, , drop = FALSE]); YY <- YY[ok, , drop = FALSE]
  B <- matrix(NA_real_, nrow(X), ncol(YY)) # placeholder only to set dimensions
  B <- matrix(NA_real_, nrow = ncol(X), ncol = length(VARS),
              dimnames = list(colnames(X), VARS))
  for (eq in seq_along(VARS)) B[, eq] <- lm.fit(X, YY[, eq])$coefficients

  extract_block <- function(prefix) {
    A <- matrix(0, length(VARS), length(VARS), dimnames = list(VARS, VARS))
    for (eq in VARS) for (v in VARS) A[eq, v] <- B[paste0(v, prefix), eq]
    A
  }
  domestic <- lapply(seq_len(p), function(L) extract_block(paste0("_L", L)))
  list(B = B, domestic = domestic, lambda0 = extract_block("_star_0"),
       lambda1 = extract_block("_star_L1"), nobs = nrow(YY),
       first_q = panel$qid[rows][ok][1], last_q = tail(panel$qid[rows][ok], 1))
}

make_global_system <- function(panel, W, p_by_country) {
  k <- length(VARS); N <- length(COUNTRIES); nk <- N * k
  max_p <- max(p_by_country)
  fits <- setNames(vector("list", N), COUNTRIES)
  for (cc in COUNTRIES) fits[[cc]] <- fit_harmonized_local(panel, cc, p_by_country[[cc]])

  G0 <- diag(nk)
  G <- lapply(seq_len(max_p), function(x) matrix(0, nk, nk))
  for (i in seq_along(COUNTRIES)) {
    cc <- COUNTRIES[i]; ii <- ((i - 1L) * k + 1L):(i * k); fit <- fits[[cc]]
    for (j in seq_along(COUNTRIES)) if (j != i) {
      jj <- ((j - 1L) * k + 1L):(j * k); wij <- W[cc, COUNTRIES[j]]
      G0[ii, jj] <- G0[ii, jj] - wij * fit$lambda0
      G[[1L]][ii, jj] <- G[[1L]][ii, jj] + wij * fit$lambda1
    }
    for (L in seq_along(fit$domestic)) G[[L]][ii, ii] <- G[[L]][ii, ii] + fit$domestic[[L]]
  }
  list(G0 = G0, G = G, fits = fits, max_p = max_p)
}

global_companion <- function(F) {
  k <- nrow(F[[1]]); p <- length(F)
  if (p == 1L) return(F[[1]])
  C <- matrix(0, k * p, k * p)
  C[seq_len(k), seq_len(k * p)] <- do.call(cbind, F)
  C[(k + 1L):(k * p), seq_len(k * (p - 1L))] <- diag(k * (p - 1L))
  C
}

macro <- read_macro(); W <- read_weights()
gpr <- read_global_series(GPR_PATH, exact = GPR_COL_EXACT, label = "GPR")
oil <- read_global_series(OIL_PATH, label = "Brent")
panel <- make_panel(macro, W, gpr, oil)

summary_rows <- list(); eigen_rows <- list(); g0_rows <- list(); lag_rows <- list()
for (name in names(P_GRIDS)) {
  p_by_country <- P_GRIDS[[name]]
  sys <- make_global_system(panel, W, p_by_country)
  rcond_g0 <- rcond(sys$G0); kappa_g0 <- kappa(sys$G0, exact = FALSE)
  if (!is.finite(rcond_g0) || rcond_g0 < 1e-10)
    stopf("G0 is singular or ill-conditioned for %s (rcond=%.3e)", name, rcond_g0)
  F <- lapply(sys$G, function(G_l) solve(sys$G0, G_l))
  C <- global_companion(F); ev <- eigen(C, only.values = TRUE)$values
  ord <- order(Mod(ev), decreasing = TRUE)
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    specification = name, max_p = sys$max_p,
    n_countries = length(COUNTRIES), n_variables = length(VARS),
    nobs = unique(vapply(sys$fits, `[[`, numeric(1), "nobs")),
    first_quarter = quarter_label(min(vapply(sys$fits, `[[`, integer(1), "first_q"))),
    last_quarter = quarter_label(max(vapply(sys$fits, `[[`, integer(1), "last_q"))),
    g0_rcond = rcond_g0, g0_kappa = kappa_g0,
    spectral_radius = max(Mod(ev)), stable = max(Mod(ev)) < 1,
    borderline = max(Mod(ev)) >= .98 & max(Mod(ev)) < 1)
  eigen_rows[[length(eigen_rows) + 1L]] <- data.frame(
    specification = name, rank = seq_along(ord), real = Re(ev[ord]),
    imaginary = Im(ev[ord]), modulus = Mod(ev[ord]), angle = Arg(ev[ord]))
  g0_rows[[length(g0_rows) + 1L]] <- data.frame(
    specification = name, rcond = rcond_g0, kappa = kappa_g0,
    dimension = nrow(sys$G0), determinant_sign = sign(det(sys$G0)))
  lag_rows[[length(lag_rows) + 1L]] <- data.frame(specification = name,
    country = COUNTRIES, p = as.integer(p_by_country[COUNTRIES]))
}

write.csv(do.call(rbind, summary_rows), file.path(OUT_DIR, "01_global_stability_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, eigen_rows), file.path(OUT_DIR, "02_global_companion_eigenvalues.csv"), row.names = FALSE)
write.csv(do.call(rbind, g0_rows), file.path(OUT_DIR, "03_G0_condition_diagnostic.csv"), row.names = FALSE)
write.csv(do.call(rbind, lag_rows), file.path(OUT_DIR, "04_lag_vector_audit.csv"), row.names = FALSE)
writeLines(c("GLOBAL GVAR STABILITY DIAGNOSTIC", "",
  "Equation: G0 x_t = sum_l G_l x_{t-l} + exogenous controls + e_t.",
  "Stability is determined from the stacked global companion matrix of solve(G0, G_l).",
  "A model is stable only when its spectral radius is strictly below one.",
  "This run reports uniform p=1 and p=2 systems; use the same code with the final",
  "country-specific lag vector before final TVP-GVAR estimation."),
  file.path(OUT_DIR, "README_global_gvar_stability.txt"))

print(do.call(rbind, summary_rows), row.names = FALSE)
