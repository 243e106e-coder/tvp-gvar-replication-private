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

source("R/BVAR_ttvp.r")
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
ext.inst <- FALSE
shrink.parm <- list(B_1 = 3, B_2 = 0.03, kappa0 = -0.1/20)

if (!is.finite(saves) || saves < 10) stop("TVPGVAR_SAVES is too small.")
if (!is.finite(burns) || burns < 10) stop("TVPGVAR_BURNS is too small.")
if (!is.finite(thin) || thin <= 0 || thin > 1) stop("TVPGVAR_THIN must lie in (0,1].")
if (!is.finite(shock_pct) || shock_pct <= 0) stop("TVPGVAR_GPR_SHOCK_PCT must be positive.")

BVAR <- cmpfun(BVAR)

sfInit(parallel = TRUE, cpus = CPU)
on.exit(try(sfStop(), silent = TRUE), add = TRUE)
sfExport(list = list(
  "mlag", "BVAR", "datahandling", "xglobal", "gW", "Daten", "cN",
  "bvartvpm", "saves", "burns", "thin", "ext.inst", "shrink.parm"
))

predDens <- sfLapply(seq_along(cN), function(i) {
  BVAR(i, gW = gW, bigx = xglobal, Daten = Daten, cN = cN,
       nsave = saves, nburn = burns, thin_chain = thin,
       ext.inst = ext.inst, parms = shrink.parm)
})
sfStop()

save(predDens, Data.setup, file = "results/predDens_gpr_structural.rda")

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
IRF_post <- array(
  NA_real_,
  c(Tirf, ncol(xglobal), nhor + 1L, n_irf_draws),
  dimnames = list(Data.setup$quarters[-1], colnames(xglobal), 0:nhor, NULL)
)

for (irep in seq_len(n_irf_draws)) {
  A.i <- lapply(A.list, function(z) z[, , , irep])
  S.i <- lapply(Sigma.posterior, function(z) z[, , , irep])
  for (tt in seq_len(Tirf)) {
    IRF_post[tt, , , irep] <- get_gpr_struct_irfa_t(
      tt = tt, draw_i = A.i, Sig_draw_i = S.i,
      x = t(xglobal), globalG = globalG, countries = cN,
      horz = nhor, normalization = "pct", shock_pct = shock_pct
    )$IRF_post
  }
  cat("IRF posterior draw", irep, "of", n_irf_draws, "complete\n")
}

save(IRF_post, file = "results/irf_gpr_structural.rda")

# Posterior summaries: 68% and 90% credible intervals.
med <- apply(IRF_post, c(1, 2, 3), median, na.rm = TRUE)
lo16 <- apply(IRF_post, c(1, 2, 3), quantile, probs = 0.16, na.rm = TRUE)
hi84 <- apply(IRF_post, c(1, 2, 3), quantile, probs = 0.84, na.rm = TRUE)
lo05 <- apply(IRF_post, c(1, 2, 3), quantile, probs = 0.05, na.rm = TRUE)
hi95 <- apply(IRF_post, c(1, 2, 3), quantile, probs = 0.95, na.rm = TRUE)

plot_dates <- c("2003Q1", "2008Q3", "2020Q1", "2022Q1")
plot_dates <- plot_dates[plot_dates %in% dimnames(IRF_post)[[1]]]
vars <- setdiff(colnames(xglobal), "US_gpr")

dd <- do.call(rbind, lapply(plot_dates, function(d) {
  ti <- match(d, dimnames(IRF_post)[[1]])
  do.call(rbind, lapply(vars, function(v) data.frame(
    date = d, variable = v, horizon = 0:nhor,
    median = med[ti, v, ], low68 = lo16[ti, v, ], high68 = hi84[ti, v, ],
    low90 = lo05[ti, v, ], high90 = hi95[ti, v, ],
    stringsAsFactors = FALSE
  )))
}))
write.csv(dd, "results/gpr_structural_irf_summary.csv", row.names = FALSE)

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
cumY <- IRF_post[, yidx, , , drop = FALSE]
for (tt in seq_len(dim(cumY)[1])) {
  for (vv in seq_len(dim(cumY)[2])) {
    for (ddraw in seq_len(dim(cumY)[4])) {
      cumY[tt, vv, , ddraw] <- cumsum(cumY[tt, vv, , ddraw])
    }
  }
}

cum_med <- apply(cumY, c(1, 2, 3), median, na.rm = TRUE)
cum_lo16 <- apply(cumY, c(1, 2, 3), quantile, probs = 0.16, na.rm = TRUE)
cum_hi84 <- apply(cumY, c(1, 2, 3), quantile, probs = 0.84, na.rm = TRUE)

y_dim_names <- dimnames(IRF_post)[[2]][yidx]
cum_summary <- do.call(rbind, lapply(plot_dates, function(d) {
  ti <- match(d, dimnames(IRF_post)[[1]])
  do.call(rbind, lapply(seq_along(y_dim_names), function(j) data.frame(
    date = d, variable = y_dim_names[j], horizon = 0:nhor,
    median = cum_med[ti, j, ], low68 = cum_lo16[ti, j, ], high68 = cum_hi84[ti, j, ],
    stringsAsFactors = FALSE
  )))
}))
write.csv(cum_summary, "results/gpr_structural_cumulative_gdp.csv", row.names = FALSE)

# A compact sign diagnostic for the current concern about y-direction.
sign_diag <- dd[grepl("_y$", dd$variable) & dd$horizon %in% c(0, 1, 4, 8, 12), ]
sign_diag$sign <- ifelse(sign_diag$median > 0, "positive", ifelse(sign_diag$median < 0, "negative", "zero"))
sign_diag$credible68_excludes_zero <- with(sign_diag, low68 > 0 | high68 < 0)
write.csv(sign_diag, "results/gdp_sign_diagnostic.csv", row.names = FALSE)

# Separate figures by historical date so different regimes are not visually mixed.
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
        title = paste0(d, ": responses to a +", shock_pct, "% structural global GPR shock (", vv, ")"),
        x = "Quarters after shock", y = "Response"
      )
    ggsave(sprintf("results/GPR_STRUCT_%s_%s.png", d, vv), p, width = 10, height = 13, dpi = 220)
  }

  zc <- cum_summary[cum_summary$date == d, ]
  p2 <- ggplot(zc, aes(horizon, median)) +
    geom_hline(yintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_ribbon(aes(ymin = low68, ymax = high68), alpha = 0.16) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~variable, scales = "free_y", ncol = 2) +
    theme_minimal(base_size = 10) +
    labs(
      title = paste0(d, ": cumulative GDP response to a +", shock_pct, "% structural global GPR shock"),
      x = "Quarters after shock", y = "Cumulative response"
    )
  ggsave(sprintf("results/GPR_STRUCT_%s_CUM_GDP.png", d), p2, width = 10, height = 13, dpi = 220)
}

run_summary <- c(
  "Structural global GPR TVP-GVAR run",
  paste0("US recursive order: ", paste(expected_us, collapse = " -> ")),
  paste0("Shock normalization: +", shock_pct, "% global GPR on impact"),
  paste0("saves=", saves, "; burns=", burns, "; thin=", thin),
  paste0("retained IRF draws=", n_irf_draws),
  paste0("IRF horizon=", nhor, " quarters"),
  paste0("IRF dates: ", paste(plot_dates, collapse = ", ")),
  "Interpret y as GDP-growth response; use cumulative GDP IRF for the implied GDP-level path."
)
writeLines(run_summary, "results/run_structural_summary.txt")
cat(paste(run_summary, collapse = "\n"), "\n")
