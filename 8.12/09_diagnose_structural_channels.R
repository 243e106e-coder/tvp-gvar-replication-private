#!/usr/bin/env Rscript

# =============================================================================
# TVP-GVAR STRUCTURAL CHANNEL ACCOUNTING (NO MCMC)
#
# Goal
# ----
# Diagnose why the saved GPR structural IRF produces:
#   (i) mostly positive cumulative equity responses, and
#   (ii) an unusually large cumulative ZA REER response.
#
# The script uses the EXACT saved posterior from the successful GDP-loglevel
# dominant-unit [GPR,VIX] run. It DOES NOT re-estimate the model.
#
# For each country equation and each saved stable posterior draw, the identity is
#
#   y_i(h)
#     = Lambda0_i * Wex_i(h)
#       + Lambda1_i * Wex_i(h-1)
#       + Theta1_i * y_i(h-1)
#       + Theta2_i * y_i(h-2)
#
# because only the dominant GL block receives the structural innovation.
#
# Therefore the response of a target country variable can be accounted for
# exactly by:
#   - GPR channel
#   - VIX channel
#   - foreign-equity channel
#   - foreign-REER channel
#   - other foreign-macro channel
#   - own-lag dynamics
#
# IMPORTANT INTERPRETATION
# ------------------------
# This is an exact EQUATION-LEVEL ACCOUNTING decomposition of the saved IRF.
# It is NOT a causal mediation decomposition, because contemporaneous Wex values
# are jointly determined inside the GVAR system.
#
# For deq (equity return/change) and de (REER change), the script also cumulatively
# sums each channel DRAW BY DRAW so the channel contributions add exactly to the
# cumulative level-effect IRF used in the preceding post-processing step.
#
# Base R only. No external packages required.
# =============================================================================

options(stringsAsFactors = FALSE)

posterior_file <- Sys.getenv(
  "TVPGVAR_POSTERIOR_RDA",
  "prior_artifact/results/predDens_dominant_gpr_vix.rda"
)
irf_file <- Sys.getenv(
  "TVPGVAR_IRF_RDA",
  "prior_artifact/results/irf_dominant_gpr_vix.rda"
)
out_dir <- Sys.getenv(
  "TVPGVAR_CHANNEL_OUT",
  "results/channel_decomposition"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

for (f in c(posterior_file, irf_file)) {
  if (!file.exists(f)) stop("Missing required saved result: ", f)
}

loaded_post <- load(posterior_file)
if (!all(c("predDens", "Data.setup") %in% loaded_post)) {
  stop(
    "Posterior RDA must contain predDens and Data.setup. Found: ",
    paste(loaded_post, collapse = ", ")
  )
}

loaded_irf <- load(irf_file)
if (!all(c("IRF_post", "stable_mask") %in% loaded_irf)) {
  stop(
    "IRF RDA must contain IRF_post and stable_mask. Found: ",
    paste(loaded_irf, collapse = ", ")
  )
}

# -----------------------------------------------------------------------------
# Minimal coefficient helpers, kept self-contained so this diagnostic does not
# depend on the current version of the IRF reconstruction file.
# -----------------------------------------------------------------------------
split_alpha <- function(A) {
  if (length(dim(A)) != 3L) stop("Single-draw ALPHA must be [time,coef,equation].")
  tags <- dimnames(A)[[2]]
  eq <- dimnames(A)[[3]]
  if (is.null(tags) || is.null(eq)) stop("ALPHA array lacks coefficient/equation dimnames.")

  block <- function(tag, optional = FALSE) {
    ii <- which(tags == tag)
    if (!length(ii)) {
      if (!optional) stop("Missing coefficient block: ", tag)
      return(array(0, dim = c(dim(A)[1], 0L, dim(A)[3])))
    }
    A[, ii, , drop = FALSE]
  }

  lag_tags <- grep("^Ylag[0-9]+$", unique(tags), value = TRUE)
  if (!length(lag_tags)) stop("No domestic lag blocks in ALPHA.")
  p <- max(as.integer(sub("^Ylag", "", lag_tags)))

  list(
    Lambda0 = block("Wex", optional = TRUE),
    Lambda = lapply(seq_len(p), function(j) block(paste0("Wexlag", j), optional = TRUE)),
    Theta = lapply(seq_len(p), function(j) block(paste0("Ylag", j))),
    p = p,
    eq = eq
  )
}

slice_coef <- function(arr, tt, nrow_out) {
  if (dim(arr)[2] == 0L) return(matrix(0, nrow = nrow_out, ncol = 0L))
  z <- arr[tt, , , drop = FALSE]
  z <- matrix(z, nrow = dim(arr)[2], ncol = nrow_out)
  t(z)
}

qsum <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(c(
      low90 = NA, low68 = NA, median = NA, high68 = NA, high90 = NA,
      mean = NA, positive_share = NA, negative_share = NA, n = 0
    ))
  }
  qq <- quantile(x, c(.05, .16, .50, .84, .95), names = FALSE)
  c(
    low90 = qq[1], low68 = qq[2], median = qq[3],
    high68 = qq[4], high90 = qq[5],
    mean = mean(x),
    positive_share = mean(x > 0),
    negative_share = mean(x < 0),
    n = length(x)
  )
}

# -----------------------------------------------------------------------------
# Validate saved architecture.
# -----------------------------------------------------------------------------
xglobal <- as.matrix(Data.setup$bigx)
storage.mode(xglobal) <- "double"
units <- Data.setup$countries
country_units <- Data.setup$country_units
globalG <- Data.setup$gW

if (is.null(colnames(xglobal))) stop("Data.setup$bigx lacks column names.")
if (!identical(tail(colnames(xglobal), 2L), c("GL_gpr", "GL_vix"))) {
  stop("Expected dominant columns GL_gpr -> GL_vix at the end of bigx.")
}
if (!all(country_units %in% units) || !"GL" %in% units) {
  stop("Saved Data.setup unit labels are inconsistent.")
}
if (length(predDens) != length(units)) stop("predDens/unit count mismatch.")

irf_dates <- dimnames(IRF_post)[[1]]
irf_vars <- dimnames(IRF_post)[[2]]
irf_h <- dimnames(IRF_post)[[3]]
if (is.null(irf_dates) || is.null(irf_vars)) stop("IRF_post lacks date/variable dimnames.")
if (!identical(irf_vars, colnames(xglobal))) {
  stop("IRF variable order differs from saved Data.setup$bigx.")
}

h_num <- suppressWarnings(as.integer(irf_h))
if (anyNA(h_num)) h_num <- seq_len(dim(IRF_post)[3]) - 1L

n_draws <- dim(IRF_post)[4]
if (!all(dim(stable_mask) == c(length(irf_dates), n_draws))) {
  stop("stable_mask dimensions do not match IRF_post.")
}
if (any(vapply(predDens, function(z) dim(z$ALPHA)[4], integer(1)) != n_draws)) {
  stop("Posterior ALPHA draw count differs from IRF_post draw count.")
}

selected_dates <- c("2003Q1", "2008Q3", "2014Q3", "2020Q1", "2022Q1", "2023Q4")
selected_dates <- intersect(selected_dates, irf_dates)
if (!length(selected_dates)) stop("None of the requested diagnostic dates exists in IRF_post.")

# Targets: every country's equity equation + ZA REER equation.
target_specs <- rbind(
  data.frame(
    country = country_units,
    target_variable = paste0(country_units, "_deq"),
    target_concept = "equity",
    stringsAsFactors = FALSE
  ),
  data.frame(
    country = "ZA",
    target_variable = "ZA_de",
    target_concept = "reer",
    stringsAsFactors = FALSE
  )
)

# Wex order is fixed by the saved dominant mapping.
expected_wex <- c(
  "foreign_y", "foreign_dp", "foreign_r", "foreign_de", "foreign_deq",
  "global_gpr", "global_vix"
)

# -----------------------------------------------------------------------------
# Main draw-level accounting.
# -----------------------------------------------------------------------------
draw_rows <- list()
coef_rows <- list()
rr <- 0L
cr <- 0L

for (ddate in selected_dates) {
  tt <- match(ddate, irf_dates)
  stable_draws <- which(stable_mask[tt, ])
  if (!length(stable_draws)) {
    warning("No stable draws at ", ddate, "; skipping.")
    next
  }

  cat(
    "Date ", ddate, ": ", length(stable_draws), "/", n_draws,
    " stable draws\n", sep = ""
  )

  for (draw in stable_draws) {
    # Exact saved full-system response for this date and posterior draw:
    # variables x horizons.
    Xresp <- IRF_post[tt, , , draw, drop = FALSE]
    Xresp <- matrix(
      Xresp,
      nrow = length(irf_vars),
      ncol = length(h_num),
      dimnames = list(irf_vars, as.character(h_num))
    )

    # Loop over 14 country units once. Each country provides its equity target,
    # and ZA additionally provides its REER target.
    for (cc in country_units) {
      unit_i <- match(cc, units)
      pd <- predDens[[unit_i]]
      A4 <- pd$ALPHA
      A_draw <- A4[, , , draw, drop = FALSE]
      A_draw <- array(
        A_draw,
        dim = dim(A4)[1:3],
        dimnames = dimnames(A4)[1:3]
      )
      V <- split_alpha(A_draw)

      Wi <- as.matrix(globalG[[cc]])
      k <- length(V$eq)
      if (k != 5L) stop("Country ", cc, " does not have 5 equations.")
      if (nrow(Wi) != 12L) stop("Country ", cc, " mapping does not have 12 rows.")

      wex_names <- rownames(Wi)[(k + 1L):nrow(Wi)]
      if (!identical(wex_names, expected_wex)) {
        stop(
          "Wex row order mismatch for ", cc, ". Found: ",
          paste(wex_names, collapse = ", ")
        )
      }

      lam0 <- slice_coef(V$Lambda0, tt, k)
      if (!all(dim(lam0) == c(k, 7L))) stop("Lambda0 shape mismatch for ", cc)

      lam1 <- slice_coef(V$Lambda[[1L]], tt, k)
      if (!all(dim(lam1) == c(k, 7L))) stop("Lambda1 shape mismatch for ", cc)

      theta1 <- slice_coef(V$Theta[[1L]], tt, k)
      if (!all(dim(theta1) == c(k, k))) stop("Theta1 shape mismatch for ", cc)

      theta2 <- if (V$p >= 2L) {
        z <- slice_coef(V$Theta[[2L]], tt, k)
        if (!all(dim(z) == c(k, k))) stop("Theta2 shape mismatch for ", cc)
        z
      } else {
        matrix(0, k, k)
      }

      # Mapping from global response to own and Wex response, for all horizons.
      own_select <- Wi[seq_len(k), , drop = FALSE]
      wex_select <- Wi[(k + 1L):nrow(Wi), , drop = FALSE]
      own_resp <- own_select %*% Xresp       # 5 x H
      wex_resp <- wex_select %*% Xresp       # 7 x H

      targets_this_cc <- subset(target_specs, country == cc)
      for (ts in seq_len(nrow(targets_this_cc))) {
        target_var <- targets_this_cc$target_variable[ts]
        target_concept <- targets_this_cc$target_concept[ts]

        eq <- match(target_var, V$eq)
        if (is.na(eq)) {
          # Some saved estimators may use suffix-only equation labels.
          suffix <- sub(paste0("^", cc, "_"), "", target_var)
          eq <- match(suffix, V$eq)
        }
        if (is.na(eq)) {
          stop(
            "Could not locate target equation ", target_var,
            " in ALPHA equations: ", paste(V$eq, collapse = ", ")
          )
        }

        target_global <- match(target_var, irf_vars)
        if (is.na(target_global)) stop("Target not found in IRF_post: ", target_var)

        # Save coefficient draws for the three focal channels plus REER.
        cr <- cr + 1L
        coef_rows[[cr]] <- data.frame(
          date = ddate,
          draw = draw,
          country = cc,
          target_variable = target_var,
          target_concept = target_concept,
          lambda0_foreign_de = lam0[eq, match("foreign_de", expected_wex)],
          lambda0_foreign_deq = lam0[eq, match("foreign_deq", expected_wex)],
          lambda0_global_gpr = lam0[eq, match("global_gpr", expected_wex)],
          lambda0_global_vix = lam0[eq, match("global_vix", expected_wex)],
          lambda1_foreign_de = lam1[eq, match("foreign_de", expected_wex)],
          lambda1_foreign_deq = lam1[eq, match("foreign_deq", expected_wex)],
          lambda1_global_gpr = lam1[eq, match("global_gpr", expected_wex)],
          lambda1_global_vix = lam1[eq, match("global_vix", expected_wex)],
          stringsAsFactors = FALSE
        )

        for (hh_i in seq_along(h_num)) {
          h <- h_num[hh_i]

          # Current Wex contributions Lambda0 * Wex(h)
          curr <- lam0[eq, ] * wex_resp[, hh_i]
          names(curr) <- expected_wex

          # Lag-1 Wex contributions Lambda1 * Wex(h-1)
          lag1 <- setNames(rep(0, length(expected_wex)), expected_wex)
          if (hh_i >= 2L) {
            lag1 <- lam1[eq, ] * wex_resp[, hh_i - 1L]
            names(lag1) <- expected_wex
          }

          # Own lag contributions.  They are grouped because the main question
          # concerns GPR/VIX/foreign-market channels.
          own_lag1 <- setNames(rep(0, k), rownames(Wi)[seq_len(k)])
          own_lag2 <- setNames(rep(0, k), rownames(Wi)[seq_len(k)])
          if (hh_i >= 2L) {
            own_lag1 <- theta1[eq, ] * own_resp[, hh_i - 1L]
            names(own_lag1) <- rownames(Wi)[seq_len(k)]
          }
          if (V$p >= 2L && hh_i >= 3L) {
            own_lag2 <- theta2[eq, ] * own_resp[, hh_i - 2L]
            names(own_lag2) <- rownames(Wi)[seq_len(k)]
          }

          # Focal grouped channels.
          gpr_channel <- curr[["global_gpr"]] + lag1[["global_gpr"]]
          vix_channel <- curr[["global_vix"]] + lag1[["global_vix"]]
          foreign_equity_channel <- curr[["foreign_deq"]] + lag1[["foreign_deq"]]
          foreign_reer_channel <- curr[["foreign_de"]] + lag1[["foreign_de"]]
          other_foreign_macro_channel <- sum(
            curr[c("foreign_y", "foreign_dp", "foreign_r")] +
              lag1[c("foreign_y", "foreign_dp", "foreign_r")]
          )
          own_dynamics_channel <- sum(own_lag1) + sum(own_lag2)

          actual <- Xresp[target_global, hh_i]
          accounted <- (
            gpr_channel +
            vix_channel +
            foreign_equity_channel +
            foreign_reer_channel +
            other_foreign_macro_channel +
            own_dynamics_channel
          )
          accounting_error <- actual - accounted

          rr <- rr + 1L
          draw_rows[[rr]] <- data.frame(
            date = ddate,
            draw = draw,
            country = cc,
            target_variable = target_var,
            target_concept = target_concept,
            horizon = h,
            actual_response = actual,

            gpr_channel = gpr_channel,
            vix_channel = vix_channel,
            foreign_equity_channel = foreign_equity_channel,
            foreign_reer_channel = foreign_reer_channel,
            other_foreign_macro_channel = other_foreign_macro_channel,
            own_dynamics_channel = own_dynamics_channel,

            gpr_current = curr[["global_gpr"]],
            gpr_lag1 = lag1[["global_gpr"]],
            vix_current = curr[["global_vix"]],
            vix_lag1 = lag1[["global_vix"]],
            foreign_equity_current = curr[["foreign_deq"]],
            foreign_equity_lag1 = lag1[["foreign_deq"]],
            foreign_reer_current = curr[["foreign_de"]],
            foreign_reer_lag1 = lag1[["foreign_de"]],

            accounted_response = accounted,
            accounting_error = accounting_error,
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

if (!length(draw_rows)) stop("No channel-decomposition rows were produced.")

draw_tab <- do.call(rbind, draw_rows)
coef_tab <- do.call(rbind, coef_rows)

# Hard numerical audit.
max_err <- max(abs(draw_tab$accounting_error), na.rm = TRUE)
if (!is.finite(max_err)) stop("Non-finite accounting error.")
if (max_err > 1e-7) {
  stop(
    "Channel accounting does not reproduce the saved IRF. max |error| = ",
    signif(max_err, 8),
    ". Do not interpret channel results until the mapping is fixed."
  )
}

# -----------------------------------------------------------------------------
# Draw-by-draw cumulative contributions for the change variables de/deq.
# -----------------------------------------------------------------------------
channel_cols <- c(
  "actual_response",
  "gpr_channel", "vix_channel",
  "foreign_equity_channel", "foreign_reer_channel",
  "other_foreign_macro_channel", "own_dynamics_channel",
  "gpr_current", "gpr_lag1",
  "vix_current", "vix_lag1",
  "foreign_equity_current", "foreign_equity_lag1",
  "foreign_reer_current", "foreign_reer_lag1",
  "accounted_response", "accounting_error"
)

group_key <- interaction(
  draw_tab$date, draw_tab$draw, draw_tab$target_variable,
  drop = TRUE, lex.order = TRUE
)
cum_tab <- draw_tab
for (idx in split(seq_len(nrow(draw_tab)), group_key)) {
  idx <- idx[order(draw_tab$horizon[idx])]
  for (nm in channel_cols) {
    cum_tab[idx, nm] <- cumsum(draw_tab[idx, nm])
  }
}

# Keep a manageable draw-level audit file.
write.csv(
  draw_tab,
  file.path(out_dir, "channel_decomposition_direct_draws_stable.csv"),
  row.names = FALSE
)
write.csv(
  cum_tab,
  file.path(out_dir, "channel_decomposition_cumulative_draws_stable.csv"),
  row.names = FALSE
)
write.csv(
  coef_tab,
  file.path(out_dir, "focal_channel_coefficients_draws_stable.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Posterior summaries.
# -----------------------------------------------------------------------------
focal_channels <- c(
  "actual_response",
  "gpr_channel",
  "vix_channel",
  "foreign_equity_channel",
  "foreign_reer_channel",
  "other_foreign_macro_channel",
  "own_dynamics_channel"
)

summarise_channels <- function(tab, representation) {
  key <- interaction(
    tab$date, tab$country, tab$target_variable, tab$target_concept, tab$horizon,
    drop = TRUE, lex.order = TRUE
  )
  rows <- list()
  rr <- 0L

  for (idx in split(seq_len(nrow(tab)), key)) {
    d <- tab[idx, , drop = FALSE]
    for (ch in focal_channels) {
      qs <- qsum(d[[ch]])
      rr <- rr + 1L
      rows[[rr]] <- data.frame(
        date = d$date[1],
        country = d$country[1],
        target_variable = d$target_variable[1],
        target_concept = d$target_concept[1],
        representation = representation,
        horizon = d$horizon[1],
        channel = ch,
        median = qs[["median"]],
        mean = qs[["mean"]],
        low68 = qs[["low68"]],
        high68 = qs[["high68"]],
        low90 = qs[["low90"]],
        high90 = qs[["high90"]],
        positive_share = qs[["positive_share"]],
        negative_share = qs[["negative_share"]],
        draws = as.integer(qs[["n"]]),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

sum_direct <- summarise_channels(draw_tab, "direct_change_IRF")
sum_cum <- summarise_channels(cum_tab, "cumulative_level_effect")
sum_channels <- rbind(sum_direct, sum_cum)

write.csv(
  sum_channels,
  file.path(out_dir, "channel_decomposition_posterior_summary.csv"),
  row.names = FALSE
)

# Coefficient posterior summaries.
coef_names <- c(
  "lambda0_foreign_de", "lambda0_foreign_deq",
  "lambda0_global_gpr", "lambda0_global_vix",
  "lambda1_foreign_de", "lambda1_foreign_deq",
  "lambda1_global_gpr", "lambda1_global_vix"
)
coef_key <- interaction(
  coef_tab$date, coef_tab$country, coef_tab$target_variable, coef_tab$target_concept,
  drop = TRUE, lex.order = TRUE
)

coef_sum_rows <- list()
rr <- 0L
for (idx in split(seq_len(nrow(coef_tab)), coef_key)) {
  d <- coef_tab[idx, , drop = FALSE]
  for (nm in coef_names) {
    qs <- qsum(d[[nm]])
    rr <- rr + 1L
    coef_sum_rows[[rr]] <- data.frame(
      date = d$date[1],
      country = d$country[1],
      target_variable = d$target_variable[1],
      target_concept = d$target_concept[1],
      coefficient = nm,
      median = qs[["median"]],
      mean = qs[["mean"]],
      low68 = qs[["low68"]],
      high68 = qs[["high68"]],
      low90 = qs[["low90"]],
      high90 = qs[["high90"]],
      positive_share = qs[["positive_share"]],
      negative_share = qs[["negative_share"]],
      draws = as.integer(qs[["n"]]),
      stringsAsFactors = FALSE
    )
  }
}
coef_summary <- do.call(rbind, coef_sum_rows)
write.csv(
  coef_summary,
  file.path(out_dir, "focal_channel_coefficients_posterior_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Compact sign-count diagnostics for equity.
# -----------------------------------------------------------------------------
equity_sum <- subset(
  sum_channels,
  target_concept == "equity" &
    representation == "cumulative_level_effect" &
    horizon %in% c(0, 1, 4, 8, 12)
)

sign_count_rows <- list()
rr <- 0L
for (idx in split(
  seq_len(nrow(equity_sum)),
  interaction(equity_sum$date, equity_sum$horizon, equity_sum$channel, drop = TRUE)
)) {
  d <- equity_sum[idx, , drop = FALSE]
  rr <- rr + 1L
  sign_count_rows[[rr]] <- data.frame(
    date = d$date[1],
    horizon = d$horizon[1],
    channel = d$channel[1],
    positive_median_countries = sum(d$median > 0, na.rm = TRUE),
    negative_median_countries = sum(d$median < 0, na.rm = TRUE),
    significant_positive_68 = sum(d$low68 > 0, na.rm = TRUE),
    significant_negative_68 = sum(d$high68 < 0, na.rm = TRUE),
    economies = nrow(d),
    stringsAsFactors = FALSE
  )
}
equity_sign_counts <- do.call(rbind, sign_count_rows)
write.csv(
  equity_sign_counts,
  file.path(out_dir, "equity_channel_sign_counts.csv"),
  row.names = FALSE
)

# Coefficient sign counts across 14 equity equations.
eq_coef <- subset(coef_summary, target_concept == "equity")
coef_count_rows <- list()
rr <- 0L
for (idx in split(
  seq_len(nrow(eq_coef)),
  interaction(eq_coef$date, eq_coef$coefficient, drop = TRUE)
)) {
  d <- eq_coef[idx, , drop = FALSE]
  rr <- rr + 1L
  coef_count_rows[[rr]] <- data.frame(
    date = d$date[1],
    coefficient = d$coefficient[1],
    positive_median_countries = sum(d$median > 0, na.rm = TRUE),
    negative_median_countries = sum(d$median < 0, na.rm = TRUE),
    positive_share_ge_0_68 = sum(d$positive_share >= .68, na.rm = TRUE),
    negative_share_ge_0_68 = sum(d$negative_share >= .68, na.rm = TRUE),
    economies = nrow(d),
    stringsAsFactors = FALSE
  )
}
coef_sign_counts <- do.call(rbind, coef_count_rows)
write.csv(
  coef_sign_counts,
  file.path(out_dir, "equity_coefficient_sign_counts.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Focus tables: 2022Q1 equity and ZA REER.
# -----------------------------------------------------------------------------
focus_h <- c(0, 1, 4, 8, 12)

focus_equity_2022 <- subset(
  sum_channels,
  date == "2022Q1" &
    target_concept == "equity" &
    representation == "cumulative_level_effect" &
    horizon %in% focus_h
)
write.csv(
  focus_equity_2022,
  file.path(out_dir, "FOCUS_2022Q1_equity_cumulative_channels.csv"),
  row.names = FALSE
)

focus_za <- subset(
  sum_channels,
  country == "ZA" &
    target_variable == "ZA_de" &
    representation == "cumulative_level_effect" &
    horizon %in% focus_h
)
write.csv(
  focus_za,
  file.path(out_dir, "FOCUS_ZA_REER_cumulative_channels.csv"),
  row.names = FALSE
)

focus_coef_2022 <- subset(
  coef_summary,
  date == "2022Q1" & target_concept == "equity"
)
write.csv(
  focus_coef_2022,
  file.path(out_dir, "FOCUS_2022Q1_equity_coefficients.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# Human-readable automated diagnosis.
# -----------------------------------------------------------------------------
get_count <- function(date, coefficient) {
  d <- subset(
    coef_sign_counts,
    date == date & coefficient == coefficient
  )
  if (!nrow(d)) return(c(pos = NA, neg = NA))
  c(pos = d$positive_median_countries[1], neg = d$negative_median_countries[1])
}

# Avoid non-standard evaluation ambiguity in the helper above.
count_at <- function(dd, cf) {
  d <- coef_sign_counts[
    coef_sign_counts$date == dd & coef_sign_counts$coefficient == cf,
    , drop = FALSE
  ]
  if (!nrow(d)) return(c(pos = NA_integer_, neg = NA_integer_))
  c(
    pos = as.integer(d$positive_median_countries[1]),
    neg = as.integer(d$negative_median_countries[1])
  )
}

diagnostic_date <- if ("2022Q1" %in% selected_dates) "2022Q1" else tail(selected_dates, 1)
vix0 <- count_at(diagnostic_date, "lambda0_global_vix")
gpr0 <- count_at(diagnostic_date, "lambda0_global_gpr")
feq0 <- count_at(diagnostic_date, "lambda0_foreign_deq")

# Channel medians at h=12 (or maximum available horizon).
h_focus <- if (12L %in% h_num) 12L else max(h_num)
ch_focus <- subset(
  equity_sign_counts,
  date == diagnostic_date & horizon == h_focus
)

line_channel <- function(ch) {
  d <- ch_focus[ch_focus$channel == ch, , drop = FALSE]
  if (!nrow(d)) return(paste0(ch, ": unavailable"))
  paste0(
    ch, ": +median in ", d$positive_median_countries[1], "/",
    d$economies[1], "; -median in ", d$negative_median_countries[1], "/",
    d$economies[1]
  )
}

report <- c(
  "TVP-GVAR STRUCTURAL CHANNEL ACCOUNTING — NO MCMC",
  "================================================",
  paste0("Posterior source: ", posterior_file),
  paste0("IRF source: ", irf_file),
  paste0("Dates: ", paste(selected_dates, collapse = ", ")),
  paste0("Max exact accounting error: ", format(max_err, scientific = TRUE, digits = 4)),
  "",
  "METHOD",
  "Country-equation IRFs are decomposed exactly into contemporaneous Wex, lag-1 Wex, and own-lag contributions.",
  "GPR/VIX/foreign-equity/foreign-REER contributions are separated by coefficient block.",
  "For deq and de, cumulative channel effects are computed draw-by-draw before posterior quantiles.",
  "This is equation-level structural accounting, not causal mediation.",
  "",
  paste0("DIAGNOSTIC DATE: ", diagnostic_date),
  paste0(
    "Equity lambda0_global_vix median sign: positive in ", vix0[["pos"]],
    "/14, negative in ", vix0[["neg"]], "/14."
  ),
  paste0(
    "Equity lambda0_global_gpr median sign: positive in ", gpr0[["pos"]],
    "/14, negative in ", gpr0[["neg"]], "/14."
  ),
  paste0(
    "Equity lambda0_foreign_deq median sign: positive in ", feq0[["pos"]],
    "/14, negative in ", feq0[["neg"]], "/14."
  ),
  "",
  paste0("CUMULATIVE EQUITY CHANNEL SIGNS AT h=", h_focus),
  line_channel("actual_response"),
  line_channel("gpr_channel"),
  line_channel("vix_channel"),
  line_channel("foreign_equity_channel"),
  line_channel("foreign_reer_channel"),
  line_channel("other_foreign_macro_channel"),
  line_channel("own_dynamics_channel"),
  "",
  "HOW TO READ THE RESULT",
  "1) If lambda0_global_vix is mostly positive, the conditional contemporaneous VIX loading itself is a prime source of the positive equity IRF.",
  "2) If lambda0_global_vix is mostly negative but the VIX channel contribution is positive, the sign reversal comes from the endogenous VIX response / system propagation rather than the raw coefficient sign.",
  "3) If GPR and VIX channels are negative but foreign-equity or own-dynamics channels are positive and dominant, the positive cumulative equity response is generated by cross-country/dynamic feedback.",
  "4) For ZA_de, compare foreign_REER, GPR, VIX and own-dynamics cumulative channels. A dominant own-dynamics channel points to ZA local coefficients rather than raw REER scaling.",
  "",
  "Do not flip raw equity signs or rescale ZA REER based on this diagnostic alone."
)

writeLines(report, file.path(out_dir, "README_channel_diagnosis.txt"))
cat(paste(report, collapse = "\n"), "\n")

cat("\nChannel decomposition completed successfully.\n")
cat("Outputs: ", out_dir, "\n", sep = "")
