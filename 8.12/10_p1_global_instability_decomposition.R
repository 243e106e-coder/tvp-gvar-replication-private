#!/usr/bin/env Rscript
# P=1 global GVAR instability decomposition
# Zeros one contemporaneous foreign-star response-variable column at a time.

source("8.12/09_global_gvar_stability_diagnostic.R", local = FALSE)
OUT_DIR <- "8.12/p1_global_instability_decomposition"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

p1 <- setNames(rep(1L, length(COUNTRIES)), COUNTRIES)
sys <- make_global_system(panel, W, p1)

score_system <- function(G0, G) {
  rc <- rcond(G0); kp <- kappa(G0, exact = FALSE)
  F <- lapply(G, function(z) solve(G0, z))
  ev <- eigen(global_companion(F), only.values = TRUE)$values
  c(rcond = rc, kappa = kp, spectral_radius = max(Mod(ev)),
    stable = as.numeric(max(Mod(ev)) < 1))
}

base <- score_system(sys$G0, sys$G)
rows <- list(data.frame(
  scenario = "baseline", country = NA_character_, variable = NA_character_,
  base_rcond = base["rcond"], base_kappa = base["kappa"],
  base_spectral_radius = base["spectral_radius"],
  rcond = base["rcond"], kappa = base["kappa"],
  spectral_radius = base["spectral_radius"],
  delta_rcond = 0, delta_kappa = 0, delta_spectral_radius = 0,
  stable = base["stable"]
))

k <- length(VARS)
for (i in seq_along(COUNTRIES)) {
  ii <- ((i - 1L) * k + 1L):(i * k)
  for (v in seq_along(VARS)) {
    G0cf <- sys$G0
    # Remove one response-variable column of Lambda_i,0 across all foreign countries.
    for (j in seq_along(COUNTRIES)) if (j != i) {
      jj <- ((j - 1L) * k + 1L):(j * k)
      G0cf[ii, jj[v]] <- 0
    }
    z <- score_system(G0cf, sys$G)
    rows[[length(rows) + 1L]] <- data.frame(
      scenario = "zero_current_foreign_star_column",
      country = COUNTRIES[i], variable = VARS[v],
      base_rcond = base["rcond"], base_kappa = base["kappa"],
      base_spectral_radius = base["spectral_radius"],
      rcond = z["rcond"], kappa = z["kappa"],
      spectral_radius = z["spectral_radius"],
      delta_rcond = z["rcond"] - base["rcond"],
      delta_kappa = z["kappa"] - base["kappa"],
      delta_spectral_radius = z["spectral_radius"] - base["spectral_radius"],
      stable = z["stable"]
    )
  }
}

out <- do.call(rbind, rows)
out$improves_stability <- out$delta_spectral_radius < 0
out <- out[order(out$delta_spectral_radius, -out$delta_rcond, na.last = TRUE), ]
out$rank_stability_improvement <- seq_len(nrow(out))

write.csv(out, file.path(OUT_DIR, "01_p1_current_foreign_star_block_decomposition.csv"), row.names = FALSE)
write.csv(head(out[out$scenario != "baseline", ], 20),
          file.path(OUT_DIR, "02_top20_stabilizing_blocks.csv"), row.names = FALSE)
writeLines(c("P=1 GLOBAL INSTABILITY DECOMPOSITION", "",
  "Each counterfactual zeros one response-variable column of a country's contemporaneous foreign-star block.",
  "A negative delta_spectral_radius means the block made the baseline system less stable.",
  "This is a diagnostic only; review the ranking before changing the model."),
  file.path(OUT_DIR, "README_p1_global_instability_decomposition.txt"))
print(head(out, 20), row.names = FALSE)
