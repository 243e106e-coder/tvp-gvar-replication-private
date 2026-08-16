#!/usr/bin/env Rscript

rm(list = ls())
set.seed(20260816)

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

BVAR <- cmpfun(BVAR)

sfInit(parallel = TRUE, cpus = CPU)
on.exit(try(sfStop(), silent = TRUE), add = TRUE)
sfExport(list = list(
  "mlag", "BVAR", "datahandling", "xglobal", "gW", "Daten", "cN",
  "bvartvpm", "saves", "burns", "thin", "ext.inst", "shrink.parm",
  "gpr_lag_mode", "tvpgvar_wex_lag"
))

predDens <- sfLapply(seq_along(cN), function(i) {
  BVAR(i, gW = gW, bigx = xglobal, Daten = Daten, cN = cN,
       nsave = saves, nburn = burns, thin_chain = thin,
       ext.inst = ext.inst, parms = shrink.parm)
})
sfStop()

save(predDens, Data.setup, file = "results/predDens_gpr_structural.rda")

# Confirm that the estimated regressor structure matches the requested GPR lag mode.
reg_diag <- do.call(rbind, lapply(seq_along(predDens), function(i) {
  rn <- dimnames(predDens[[i]]$ALPHA)[[2]]
  data.frame(
    country = cN[[i]],
    Wex = sum(rn == "Wex"),
    Wexlag1 = sum(rn == "Wexlag1"),
    Ylag1 = sum(rn == "Ylag1"),
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

if (any(reg_diag$Wex != expected_wex) ||
    any(reg_diag$Wexlag1 != expected_lag) ||
    any(reg_diag$Ylag1 != expected_ylag)) {
  print(reg_diag)
  stop("Estimated ALPHA regressor structure does not match requested GPR lag specification.")
}
reg_diag$gpr_lag_mode <- gpr_lag_mode
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

Tirf <- nrow(xglobal) - 1L
irf_dates <- Data.setup$quarters[-1]

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
  }
  cat("IRF/stability posterior draw", irep, "of", n_irf_draws, "complete\n")
}

if (any(!is.finite(stability_rho))) {
  stop("Non-finite spectral radius detected in stability audit.")
}

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
  paste0("saves=", saves, "; burns=", burns, "; thin=", thin),
  paste0("retained IRF draws=", n_irf_draws),
  paste0("IRF horizon=", nhor, " quarters"),
  paste0("IRF dates: ", paste(plot_dates, collapse = ", ")),
  paste0("Stability rule: spectral radius < 1; near-unit warning threshold=",
         near_unit_threshold),
  paste0("Overall stable share=", sprintf("%.1f%%", 100 * overall_stable_share),
         "; near-unit share=", sprintf("%.1f%%", 100 * overall_near_share),
         "; unstable share=", sprintf("%.1f%%", 100 * overall_unstable_share)),
  "Selected-date stability:",
  selected_stability_lines,
  "All-draw and stable-only IRF summaries are both saved; near-unit stable draws are flagged but not discarded.",
  "Interpret y as GDP-growth response; use cumulative GDP IRF for the implied GDP-level path."
)
writeLines(run_summary, "results/run_structural_summary.txt")
cat(paste(run_summary, collapse = "\n"), "\n")
