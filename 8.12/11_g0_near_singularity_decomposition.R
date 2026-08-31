#!/usr/bin/env Rscript
# Diagnostic only: identifies the source of near-singularity in G0 for p=1.
# It does not alter estimation, weights, data, or the formal model.

source("8.12/09_global_gvar_stability_diagnostic.R", local = FALSE)
OUT_DIR <- "8.12/g0_near_singularity_decomposition"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

p1 <- setNames(rep(1L, length(COUNTRIES)), COUNTRIES)
sys <- make_global_system(panel, W, p1)
k <- length(VARS)

score_system <- function(G0, G) {
  rc <- rcond(G0); kp <- kappa(G0, exact = FALSE)
  F <- lapply(G, function(z) solve(G0, z))
  rho <- max(Mod(eigen(global_companion(F), only.values = TRUE)$values))
  c(rcond = rc, kappa = kp, spectral_radius = rho, stable = as.numeric(rho < 1))
}

base <- score_system(sys$G0, sys$G)
s <- svd(sys$G0, nu = nrow(sys$G0), nv = ncol(sys$G0))
sv <- data.frame(rank = seq_along(s$d), singular_value = s$d,
                 relative_to_largest = s$d / max(s$d),
                 condition_number_from_svd = max(s$d) / s$d)
write.csv(sv, file.path(OUT_DIR, "01_G0_singular_values.csv"), row.names = FALSE)

# Signed vector loadings identify the direction of near-singularity; absolute
# loadings permit country-family aggregation without cancellation.
labels <- as.vector(outer(COUNTRIES, VARS, paste, sep = "_"))
tail_n <- min(10L, ncol(s$u))
loads <- do.call(rbind, lapply(seq_len(tail_n), function(q) {
  ix <- ncol(s$u) - q + 1L
  data.frame(singular_rank = ix, singular_value = s$d[ix],
             country = rep(COUNTRIES, each = k), variable = rep(VARS, length(COUNTRIES)),
             label = labels, left_loading = s$u[, ix], right_loading = s$v[, ix],
             abs_left_loading = abs(s$u[, ix]), abs_right_loading = abs(s$v[, ix]))
}))
write.csv(loads, file.path(OUT_DIR, "02_small_singular_vector_loadings.csv"), row.names = FALSE)

aggregate_loadings <- function(x, by) {
  z <- aggregate(x[, c("abs_left_loading", "abs_right_loading")], x[by], sum)
  z[order(-pmax(z$abs_left_loading, z$abs_right_loading)), ]
}
write.csv(aggregate_loadings(loads, c("singular_rank", "singular_value", "country")),
          file.path(OUT_DIR, "03_small_singular_vector_country_loadings.csv"), row.names = FALSE)
write.csv(aggregate_loadings(loads, c("singular_rank", "singular_value", "variable")),
          file.path(OUT_DIR, "04_small_singular_vector_variable_loadings.csv"), row.names = FALSE)

# A family switch-off removes that response-variable column from every country's
# contemporaneous foreign-star block. Country switch-off removes its full Lambda_i,0 block.
rows <- list(data.frame(scenario = "baseline", target = "baseline", target_type = "none",
                        base_rcond = base["rcond"], base_kappa = base["kappa"],
                        base_spectral_radius = base["spectral_radius"],
                        rcond = base["rcond"], kappa = base["kappa"],
                        spectral_radius = base["spectral_radius"], stable = base["stable"]))
record <- function(G0cf, target, target_type) {
  z <- score_system(G0cf, sys$G)
  rows[[length(rows) + 1L]] <<- data.frame(scenario = "zero_contemporaneous_foreign_star",
    target = target, target_type = target_type, base_rcond = base["rcond"], base_kappa = base["kappa"],
    base_spectral_radius = base["spectral_radius"], rcond = z["rcond"], kappa = z["kappa"],
    spectral_radius = z["spectral_radius"], stable = z["stable"])
}
for (v in seq_along(VARS)) {
  G0cf <- sys$G0
  for (i in seq_along(COUNTRIES)) {
    ii <- ((i - 1L) * k + 1L):(i * k)
    for (j in seq_along(COUNTRIES)) if (j != i) {
      jj <- ((j - 1L) * k + 1L):(j * k)
      G0cf[ii, jj[v]] <- 0
    }
  }
  record(G0cf, VARS[v], "variable_family")
}
for (i in seq_along(COUNTRIES)) {
  G0cf <- sys$G0; ii <- ((i - 1L) * k + 1L):(i * k)
  for (j in seq_along(COUNTRIES)) if (j != i) {
    jj <- ((j - 1L) * k + 1L):(j * k)
    G0cf[ii, jj] <- 0
  }
  record(G0cf, COUNTRIES[i], "country_block")
}
out <- do.call(rbind, rows)
out$delta_rcond <- out$rcond - out$base_rcond
out$delta_kappa <- out$kappa - out$base_kappa
out$delta_spectral_radius <- out$spectral_radius - out$base_spectral_radius
out$improves_stability <- out$delta_spectral_radius < 0
out <- out[order(out$target_type, out$delta_spectral_radius, -out$delta_rcond, na.last = TRUE), ]
write.csv(out, file.path(OUT_DIR, "05_G0_family_and_country_switch_off.csv"), row.names = FALSE)
write.csv(out[out$target_type == "variable_family", ], file.path(OUT_DIR, "06_variable_family_summary.csv"), row.names = FALSE)
write.csv(out[out$target_type == "country_block", ], file.path(OUT_DIR, "07_country_block_summary.csv"), row.names = FALSE)
writeLines(c("G0 NEAR-SINGULARITY DECOMPOSITION (P=1)", "",
  "This is diagnostic-only. No formal-model coefficient, lag, weight, data, or identification change is made.",
  "Files 01--04 locate the small-SVD directions; files 05--07 test whole foreign-star families and whole country Lambda_i,0 blocks.",
  "Interpret a large improvement in rcond / fall in kappa as evidence for the source of G0 ill-conditioning.",
  "Do not use any switch-off as a final specification without a separate economic and estimation justification."),
  file.path(OUT_DIR, "README_g0_near_singularity_decomposition.txt"))
print(out, row.names = FALSE)
