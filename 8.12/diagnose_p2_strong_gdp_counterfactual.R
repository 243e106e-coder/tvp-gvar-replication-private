#!/usr/bin/env Rscript

# Counterfactual timing diagnostic for the existing P2_strong posterior.
# IMPORTANT: This script DOES NOT re-estimate MCMC.
# It keeps the original posterior draws and original stable-draw selection,
# then changes only selected contemporaneous GDP channels in the assembled
# GVAR impact matrix G to identify what generates positive GDP impact IRFs.
#
# Modes:
#   1) baseline_pure_GPR
#      - Original contemporaneous G matrix.
#      - Pure US GPR reduced-form innovation (no automatic US GDP residual shock).
#   2) no_direct_GPRt_in_nonUS_GDP
#      - Set current GPR_t coefficient in each non-US GDP equation to zero.
#   3) no_current_foreignYt_in_all_GDP
#      - Set current foreign GDP (y*_t) coefficient in every GDP equation to zero.
#   4) no_direct_GPRt_and_no_current_foreignYt_in_GDP
#      - Apply both restrictions simultaneously.
#
# This is a diagnostic/counterfactual exercise, not the final re-estimated model.

suppressPackageStartupMessages(library(ggplot2))

src <- Sys.getenv("TVPGVAR_CF_SOURCE_DIR", "source-artifact")
out <- Sys.getenv("TVPGVAR_CF_OUT", "gdp-counterfactual-diagnostic")
dates <- trimws(strsplit(
  Sys.getenv("TVPGVAR_CF_DATES",
             "2003Q1,2008Q3,2014Q3,2020Q1,2022Q1,2023Q4"),
  ",", fixed = TRUE
)[[1]])
shock_pct <- as.numeric(Sys.getenv("TVPGVAR_CF_SHOCK_PCT", "10"))
horizon <- as.integer(Sys.getenv("TVPGVAR_CF_HORIZON", "12"))
near_unit <- as.numeric(Sys.getenv("TVPGVAR_NEAR_UNIT_THRESHOLD", "0.98"))

if (!is.finite(shock_pct) || shock_pct <= 0) stop("shock_pct must be positive.")
if (!is.finite(horizon) || horizon < 0L) stop("horizon must be >= 0.")

dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "plots"), recursive = TRUE, showWarnings = FALSE)

find_one <- function(name) {
  z <- list.files(
    src,
    pattern = paste0("^", gsub("\\.", "\\\\.", name), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(z) != 1L) stop("Expected exactly one ", name, "; found ", length(z))
  z[[1]]
}

readc <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
}

load(find_one("predDens_gpr_structural.rda"))
parms <- readc(find_one("stability_profile_parameters.csv"))
stab_source <- readc(find_one("stability_summary_selected_dates.csv"))

stopifnot(exists("predDens"), exists("Data.setup"))
if (nrow(parms) != 1L || parms$profile[[1]] != "P2_strong") {
  stop("Wrong artifact: expected P2_strong.")
}

p <- as.integer(parms$p[[1]])
q <- as.integer(parms$q[[1]])
if (p != 2L || q != 1L) {
  stop("Expected GVAR(2,1), found p=", p, ", q=", q)
}

xglobal <- Data.setup$bigx
countries <- as.character(Data.setup$countries)
quarters <- as.character(Data.setup$quarters)
Wlist <- lapply(predDens, function(z) z$W)
Alist <- lapply(predDens, function(z) z$ALPHA)
Slist <- lapply(predDens, function(z) z$SIGMApost)

us_i <- match("US", countries)
gpr_i <- match("US_gpr", colnames(xglobal))
if (is.na(us_i) || is.na(gpr_i)) stop("US / US_gpr not found.")

ndraw <- dim(Alist[[1]])[4]
irf_dates <- quarters[-seq_len(p)]
date_idx <- match(dates, irf_dates)
if (anyNA(date_idx)) {
  stop("Some requested dates are not in the posterior path: ",
       paste(dates[is.na(date_idx)], collapse = ", "))
}

modes <- c(
  "baseline_pure_GPR",
  "no_direct_GPRt_in_nonUS_GDP",
  "no_current_foreignYt_in_all_GDP",
  "no_direct_GPRt_and_no_current_foreignYt_in_GDP"
)

split_alpha <- function(A) {
  tags <- dimnames(A)[[2]]
  get_block <- function(tag) {
    ii <- which(tags == tag)
    if (!length(ii)) stop("Missing coefficient block: ", tag)
    A[, ii, , drop = FALSE]
  }

  lag_tags <- grep("^Ylag[0-9]+$", unique(tags), value = TRUE)
  lag_n <- as.integer(sub("^Ylag", "", lag_tags))
  pp <- max(lag_n)

  list(
    L0 = get_block("Wex"),
    Theta = lapply(seq_len(pp), function(k) get_block(paste0("Ylag", k))),
    Llag = lapply(seq_len(pp), function(k) {
      tag <- paste0("Wexlag", k)
      if (any(tags == tag)) get_block(tag) else NULL
    })
  )
}

align_lag <- function(B, Wi, k, cc) {
  B <- as.matrix(B)
  need <- nrow(Wi) - k
  if (ncol(B) == need) return(B)

  # Compatibility with the existing structural runner:
  # non-US lagged foreign block can omit lagged global GPR.
  if (cc != "US" && need == 6L && ncol(B) == 5L) {
    return(cbind(B, rep(0, nrow(B))))
  }

  stop("Foreign-lag dimension mismatch for ", cc,
       ": need ", need, ", got ", ncol(B))
}

rho_companion <- function(F) {
  K <- nrow(F[[1]])
  pp <- length(F)
  if (pp == 1L) {
    return(max(Mod(eigen(F[[1]], only.values = TRUE)$values)))
  }
  top <- do.call(cbind, F)
  lower <- cbind(
    diag(K * (pp - 1L)),
    matrix(0, K * (pp - 1L), K)
  )
  max(Mod(eigen(rbind(top, lower), only.values = TRUE)$values))
}

stab_cov <- function(S) {
  S <- (as.matrix(S) + t(as.matrix(S))) / 2
  ev <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  mn <- min(ev)
  if (mn <= 1e-10) {
    S <- S + diag(abs(mn) + 1e-10, nrow(S))
  }
  S
}

# Build the global contemporaneous matrix G and lag matrices F.
# Counterfactual restrictions are imposed ONLY on lambda0 before G is formed.
build_state <- function(tt, Adraw, Sdraw, mode = "baseline_pure_GPR") {
  if (!mode %in% modes) stop("Unknown mode: ", mode)

  Grow <- vector("list", length(countries))
  Hcountry <- vector("list", length(countries))
  S <- vector("list", length(countries))

  changed_direct_gpr <- 0L
  changed_foreign_y <- 0L

  for (i in seq_along(countries)) {
    cc <- countries[[i]]
    aa <- split_alpha(Adraw[[i]])
    Wi <- as.matrix(Wlist[[i]])
    k <- dim(aa$L0)[3]

    lambda0 <- t(aa$L0[tt, , ])

    # Current Wex ordering used by the existing model:
    #   US:     foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq
    #   non-US: foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq, global_gpr
    # Domestic equation order starts with y.
    if (mode %in% c(
      "no_direct_GPRt_in_nonUS_GDP",
      "no_direct_GPRt_and_no_current_foreignYt_in_GDP"
    )) {
      if (cc != "US") {
        if (ncol(lambda0) != 6L) {
          stop("Expected 6 current Wex columns for non-US country ", cc,
               "; found ", ncol(lambda0))
        }
        if (lambda0[1, 6] != 0) changed_direct_gpr <- changed_direct_gpr + 1L
        lambda0[1, 6] <- 0
      }
    }

    if (mode %in% c(
      "no_current_foreignYt_in_all_GDP",
      "no_direct_GPRt_and_no_current_foreignYt_in_GDP"
    )) {
      if (ncol(lambda0) < 1L) stop("No current foreign-y column for ", cc)
      if (lambda0[1, 1] != 0) changed_foreign_y <- changed_foreign_y + 1L
      lambda0[1, 1] <- 0
    }

    A0 <- cbind(diag(k), -lambda0)
    Grow[[i]] <- A0 %*% Wi

    Hcountry[[i]] <- lapply(seq_along(aa$Theta), function(kL) {
      theta <- t(aa$Theta[[kL]][tt, , ])
      if (is.null(aa$Llag[[kL]])) {
        lagfull <- matrix(0, nrow = k, ncol = nrow(Wi) - k)
      } else {
        lagfull <- align_lag(t(aa$Llag[[kL]][tt, , ]), Wi, k, cc)
      }
      cbind(theta, lagfull) %*% Wi
    })

    S[[i]] <- Sdraw[[i]][tt, , ]
  }

  G <- do.call(rbind, Grow)
  H <- lapply(seq_len(length(Hcountry[[1]])), function(kL) {
    do.call(rbind, lapply(Hcountry, function(z) z[[kL]]))
  })
  F <- lapply(H, function(h) solve(G, h))

  list(
    G = G,
    F = F,
    S = S,
    rho = rho_companion(F),
    changed_direct_gpr = changed_direct_gpr,
    changed_foreign_y = changed_foreign_y
  )
}

# Pure GPR reduced-form innovation in the US block.
# US structural ordering in the current model is:
#   GPR, y, dp, r, de, deq
make_pure_gpr_shock <- function(Sus) {
  Sus <- stab_cov(Sus)
  L <- t(chol(Sus))
  pure <- rep(0, nrow(L))
  pure[1] <- L[1, 1]
  pure
}

impact_from_us_shock <- function(G, Sblocks, usshock) {
  u <- unlist(lapply(seq_along(Sblocks), function(i) {
    if (i == us_i) usshock else rep(0, nrow(Sblocks[[i]]))
  }), use.names = FALSE)

  imp <- as.numeric(solve(G, u))
  names(imp) <- colnames(xglobal)

  target <- log1p(shock_pct / 100)
  if (!is.finite(imp[gpr_i]) || abs(imp[gpr_i]) < 1e-12) {
    stop("Cannot scale shock: contemporaneous GPR response is zero/non-finite.")
  }

  imp * target / imp[gpr_i]
}

propagate <- function(F, impact, H = 12L) {
  K <- length(impact)
  ans <- matrix(
    0,
    K,
    H + 1L,
    dimnames = list(names(impact), as.character(0:H))
  )
  ans[, 1] <- impact

  if (H >= 1L) {
    for (h in seq_len(H)) {
      for (lag in seq_len(min(length(F), h))) {
        ans[, h + 1] <- ans[, h + 1] +
          F[[lag]] %*% ans[, h - lag + 1]
      }
    }
  }

  ans
}

rows_irf <- list()
rows_stab <- list()
rows_rho <- list()
ni <- ns <- nr <- 0L

for (dd in seq_along(dates)) {
  d <- dates[[dd]]
  tt <- date_idx[[dd]]
  stable_n <- 0L
  near_n <- 0L

  message("Date: ", d)

  for (draw in seq_len(ndraw)) {
    Adraw <- lapply(Alist, function(z) z[, , , draw])
    Sdraw <- lapply(Slist, function(z) z[, , , draw])

    # IMPORTANT: stable-draw selection is based on the ORIGINAL model only.
    st0 <- build_state(tt, Adraw, Sdraw, "baseline_pure_GPR")
    if (!is.finite(st0$rho) || st0$rho >= 1) next

    stable_n <- stable_n + 1L
    if (st0$rho >= near_unit) near_n <- near_n + 1L

    pure_shock <- make_pure_gpr_shock(st0$S[[us_i]])

    for (mode in modes) {
      st <- if (mode == "baseline_pure_GPR") {
        st0
      } else {
        build_state(tt, Adraw, Sdraw, mode)
      }

      nr <- nr + 1L
      rows_rho[[nr]] <- data.frame(
        date = d,
        draw = draw,
        mode = mode,
        original_rho = st0$rho,
        counterfactual_rho = st$rho,
        counterfactual_stable = is.finite(st$rho) && st$rho < 1,
        changed_direct_gpr_equations = st$changed_direct_gpr,
        changed_foreign_y_equations = st$changed_foreign_y,
        stringsAsFactors = FALSE
      )

      # Keep the draw in the comparison even if the counterfactual F is unstable.
      # h=0 remains interpretable from G. For h>0, unstable counterfactual paths are
      # still exported and explicitly flagged by counterfactual_stable.
      imp <- impact_from_us_shock(st$G, st$S, pure_shock)
      path <- propagate(st$F, imp, H = horizon)

      for (cc in countries) {
        yn <- paste0(cc, "_y")
        for (h in 0:horizon) {
          ni <- ni + 1L
          rows_irf[[ni]] <- data.frame(
            date = d,
            draw = draw,
            country = cc,
            mode = mode,
            horizon = h,
            original_rho = st0$rho,
            counterfactual_rho = st$rho,
            counterfactual_stable = is.finite(st$rho) && st$rho < 1,
            y_response = path[yn, h + 1],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }

  ns <- ns + 1L
  rows_stab[[ns]] <- data.frame(
    date = d,
    draws = ndraw,
    stable_draws_recomputed = stable_n,
    near_unit_draws_recomputed = near_n,
    stringsAsFactors = FALSE
  )
}

irf_raw <- do.call(rbind, rows_irf)
rho_raw <- do.call(rbind, rows_rho)
stab_re <- do.call(rbind, rows_stab)

# 1) Reproduce original stable counts before interpreting the counterfactual.
stab_source$date <- as.character(stab_source$date)
check <- merge(
  stab_re,
  stab_source[, c("date", "stable_draws", "near_unit_draws")],
  by = "date",
  all.x = TRUE,
  sort = FALSE
)
check$stable_match <- check$stable_draws_recomputed == as.integer(check$stable_draws)
check$near_unit_match <- check$near_unit_draws_recomputed == as.integer(check$near_unit_draws)
write.csv(
  check,
  file.path(out, "01_stability_reconstruction_check.csv"),
  row.names = FALSE
)
if (any(!check$stable_match)) {
  stop("Recomputed original stable-draw counts do not match the source run.")
}

qv <- function(x, prob) {
  as.numeric(quantile(x, prob, names = FALSE, type = 8, na.rm = TRUE))
}

summ <- function(df, value, groups) {
  key <- interaction(df[groups], drop = TRUE, lex.order = TRUE)
  z <- lapply(split(df, key), function(a) {
    x <- a[[value]]
    x <- x[is.finite(x)]
    first <- a[1, groups, drop = FALSE]
    cbind(
      first,
      data.frame(
        draws = length(x),
        median = if (length(x)) median(x) else NA_real_,
        low68 = if (length(x)) qv(x, 0.16) else NA_real_,
        high68 = if (length(x)) qv(x, 0.84) else NA_real_,
        low90 = if (length(x)) qv(x, 0.05) else NA_real_,
        high90 = if (length(x)) qv(x, 0.95) else NA_real_,
        positive_share = if (length(x)) mean(x > 0) else NA_real_,
        stringsAsFactors = FALSE
      )
    )
  })
  ans <- do.call(rbind, z)
  rownames(ans) <- NULL
  ans
}

# 2) Counterfactual stability implications, keeping original stable draws fixed.
rho_summary <- do.call(rbind, lapply(
  split(rho_raw, interaction(rho_raw$date, rho_raw$mode, drop = TRUE)),
  function(a) {
    data.frame(
      date = a$date[[1]],
      mode = a$mode[[1]],
      original_stable_draws = nrow(a),
      counterfactual_stable_draws = sum(a$counterfactual_stable),
      counterfactual_stable_share_within_original_stable = mean(a$counterfactual_stable),
      median_counterfactual_rho = median(a$counterfactual_rho[is.finite(a$counterfactual_rho)]),
      p95_counterfactual_rho = qv(a$counterfactual_rho[is.finite(a$counterfactual_rho)], 0.95),
      stringsAsFactors = FALSE
    )
  }
))
rownames(rho_summary) <- NULL
write.csv(
  rho_summary,
  file.path(out, "02_counterfactual_stability_summary.csv"),
  row.names = FALSE
)

# 3) GDP IRF summaries using the SAME original stable draws in every mode.
gdp <- summ(
  irf_raw,
  "y_response",
  c("date", "country", "mode", "horizon")
)
write.csv(
  gdp,
  file.path(out, "03_GDP_counterfactual_IRF_summary.csv"),
  row.names = FALSE
)

# 4) h=0 comparison in long form.
gdp_h0 <- gdp[gdp$horizon == 0, ]
write.csv(
  gdp_h0,
  file.path(out, "04_GDP_h0_counterfactual_comparison.csv"),
  row.names = FALSE
)

# 5) Wide h=0 table: baseline and the three counterfactuals side-by-side.
wide_piece <- function(mode, prefix) {
  z <- gdp_h0[gdp_h0$mode == mode,
              c("date", "country", "median", "positive_share")]
  names(z)[3:4] <- paste0(prefix, c("_median", "_positive_share"))
  z
}

wide <- wide_piece("baseline_pure_GPR", "baseline")
wide <- merge(
  wide,
  wide_piece("no_direct_GPRt_in_nonUS_GDP", "no_direct_gpr"),
  by = c("date", "country"), all = TRUE
)
wide <- merge(
  wide,
  wide_piece("no_current_foreignYt_in_all_GDP", "no_foreign_y"),
  by = c("date", "country"), all = TRUE
)
wide <- merge(
  wide,
  wide_piece("no_direct_GPRt_and_no_current_foreignYt_in_GDP", "both_off"),
  by = c("date", "country"), all = TRUE
)

wide$no_direct_gpr_minus_baseline <- wide$no_direct_gpr_median - wide$baseline_median
wide$no_foreign_y_minus_baseline <- wide$no_foreign_y_median - wide$baseline_median
wide$both_off_minus_baseline <- wide$both_off_median - wide$baseline_median
wide$baseline_positive <- wide$baseline_median > 0
wide$no_direct_gpr_positive <- wide$no_direct_gpr_median > 0
wide$no_foreign_y_positive <- wide$no_foreign_y_median > 0
wide$both_off_positive <- wide$both_off_median > 0
wide$flip_after_no_direct_gpr <- wide$baseline_positive & !wide$no_direct_gpr_positive
wide$flip_after_no_foreign_y <- wide$baseline_positive & !wide$no_foreign_y_positive
wide$flip_after_both_off <- wide$baseline_positive & !wide$both_off_positive

write.csv(
  wide,
  file.path(out, "05_GDP_h0_channel_decomposition_wide.csv"),
  row.names = FALSE
)

# 6) Date-level sign diagnostic: how many country medians remain positive?
sign_summary <- do.call(rbind, lapply(
  split(gdp_h0, interaction(gdp_h0$date, gdp_h0$mode, drop = TRUE)),
  function(a) {
    data.frame(
      date = a$date[[1]],
      mode = a$mode[[1]],
      countries = nrow(a),
      positive_median_countries = sum(a$median > 0, na.rm = TRUE),
      zero_or_negative_median_countries = sum(a$median <= 0, na.rm = TRUE),
      positive_share_of_countries = mean(a$median > 0, na.rm = TRUE),
      median_abs_h0 = median(abs(a$median), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
rownames(sign_summary) <- NULL
write.csv(
  sign_summary,
  file.path(out, "06_GDP_h0_positive_country_share.csv"),
  row.names = FALSE
)

# 7) Country-level average impact across the requested dates.
country_summary <- do.call(rbind, lapply(
  split(gdp_h0, interaction(gdp_h0$country, gdp_h0$mode, drop = TRUE)),
  function(a) {
    data.frame(
      country = a$country[[1]],
      mode = a$mode[[1]],
      dates = nrow(a),
      mean_h0_median = mean(a$median, na.rm = TRUE),
      median_h0_median = median(a$median, na.rm = TRUE),
      positive_dates = sum(a$median > 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }
))
rownames(country_summary) <- NULL
write.csv(
  country_summary,
  file.path(out, "07_GDP_h0_country_summary.csv"),
  row.names = FALSE
)

# 8) Compact decision table for the exact question: which channel kills the positive h=0 sign?
decision <- aggregate(
  cbind(
    flip_after_no_direct_gpr,
    flip_after_no_foreign_y,
    flip_after_both_off
  ) ~ date,
  data = wide,
  FUN = sum
)
names(decision)[2:4] <- c(
  "countries_flip_nonpositive_when_direct_GPR_removed",
  "countries_flip_nonpositive_when_foreignY_removed",
  "countries_flip_nonpositive_when_both_removed"
)
write.csv(
  decision,
  file.path(out, "08_channel_sign_flip_decision_table.csv"),
  row.names = FALSE
)

# Plots ---------------------------------------------------------------
plot_df <- gdp_h0
plot_df$mode <- factor(plot_df$mode, levels = modes)
plot_df$country <- factor(plot_df$country, levels = rev(countries))
plot_df$date <- factor(plot_df$date, levels = dates)

p_heat <- ggplot(plot_df, aes(x = date, y = country, fill = median)) +
  geom_tile() +
  facet_wrap(~ mode, ncol = 1) +
  labs(
    title = "GDP h=0 response: contemporaneous-channel counterfactual",
    subtitle = paste0("Existing P2_strong posterior; +", shock_pct,
                      "% GPR; original stable draws held fixed"),
    x = NULL, y = NULL, fill = "Median h=0"
  ) +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  file.path(out, "plots", "GDP_h0_counterfactual_heatmap.png"),
  p_heat,
  width = 11,
  height = 16,
  dpi = 160
)

# Plot median IRFs for a focused set of countries; use all if fewer are found.
focus <- intersect(c("US", "CN", "JP", "EA", "KR", "TR"), countries)
if (!length(focus)) focus <- countries
irf_plot_df <- gdp[gdp$country %in% focus, ]
irf_plot_df$mode <- factor(irf_plot_df$mode, levels = modes)
irf_plot_df$date <- factor(irf_plot_df$date, levels = dates)

p_irf <- ggplot(
  irf_plot_df,
  aes(x = horizon, y = median, linetype = mode, group = mode)
) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(linewidth = 0.6) +
  facet_grid(country ~ date, scales = "free_y") +
  labs(
    title = "GDP median IRF under contemporaneous-channel counterfactuals",
    subtitle = paste0("Existing P2_strong posterior; +", shock_pct,
                      "% GPR; original stable draws held fixed"),
    x = "Horizon", y = "GDP response", linetype = "Mode"
  ) +
  theme_minimal(base_size = 9)

ggsave(
  file.path(out, "plots", "GDP_counterfactual_IRF_focus.png"),
  p_irf,
  width = 18,
  height = 15,
  dpi = 150
)

# Human-readable summary ------------------------------------------------
summary_lines <- character(0)
summary_lines <- append(summary_lines, "P2_strong GDP contemporaneous-channel counterfactual diagnostic")
summary_lines <- append(summary_lines, "NO MCMC re-estimation was performed.")
summary_lines <- append(summary_lines, paste0("Source profile: p=", p, ", q=", q))
summary_lines <- append(summary_lines, paste0("Shock: +", shock_pct, "% global GPR"))
summary_lines <- append(summary_lines, paste0("Horizon: 0-", horizon))
summary_lines <- append(summary_lines, paste0("Dates: ", paste(dates, collapse = ", ")))
summary_lines <- append(summary_lines, "")
summary_lines <- append(summary_lines, "Stable draws are selected from the ORIGINAL P2_strong system and held fixed across counterfactual modes.")
summary_lines <- append(summary_lines, "Mode 1: baseline pure-GPR residual shock.")
summary_lines <- append(summary_lines, "Mode 2: set current GPR_t -> GDP coefficient to zero in every non-US GDP equation.")
summary_lines <- append(summary_lines, "Mode 3: set current foreign GDP y*_t -> GDP coefficient to zero in every country's GDP equation.")
summary_lines <- append(summary_lines, "Mode 4: impose both restrictions.")
summary_lines <- append(summary_lines, "")
summary_lines <- append(summary_lines, "Interpret 08_channel_sign_flip_decision_table.csv first.")
summary_lines <- append(summary_lines, "If mode 2 removes most positive h=0 GDP signs, the direct contemporaneous GPR channel is the main source.")
summary_lines <- append(summary_lines, "If mode 2 does little but mode 3/mode 4 removes them, contemporaneous cross-country GDP feedback is the main source.")
summary_lines <- append(summary_lines, "If signs remain positive even in mode 4, inspect other contemporaneous channels (prices/rates/FX/equity) before changing the estimated model.")

writeLines(
  summary_lines,
  con = file.path(out, "README_counterfactual.txt")
)
cat(paste(summary_lines, collapse = "\n"), "\n")

message("Done. Outputs written to: ", out)
