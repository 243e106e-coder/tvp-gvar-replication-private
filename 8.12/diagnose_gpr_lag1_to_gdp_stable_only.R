#!/usr/bin/env Rscript

# =============================================================================
# Diagnose direct GPR_{t-1} -> GDP coefficients in the P2 strong
# GDP-lagged-GPR TVP-GVAR posterior.
#
# IMPORTANT:
#   * NO MCMC re-estimation.
#   * Uses the exact stable_mask from irf_gpr_structural.rda.
#   * Non-US GDP: extracts the 6th Wexlag1 coefficient (global GPR_{t-1}).
#   * US GDP: extracts US_gpr from the Ylag1 domestic-lag block.
#   * Optionally merges h=1 GDP IRFs if gdp_sign_diagnostic_stable_only.csv
#     is available in the source artifact.
#
# Expected source artifact files:
#   predDens_gpr_structural.rda
#   irf_gpr_structural.rda
# Optional:
#   gdp_lagged_gpr_test_parameters.csv
#   gdp_sign_diagnostic_stable_only.csv
#
# Usage example:
#   TVPGVAR_GDP_LAG_COEF_SOURCE_DIR="source-artifact" \
#   TVPGVAR_GDP_LAG_COEF_OUT="gdp-gpr-lag-coef-diagnostic" \
#   Rscript --vanilla 8.12/diagnose_gpr_lag1_to_gdp_stable_only.R
# =============================================================================

src <- Sys.getenv("TVPGVAR_GDP_LAG_COEF_SOURCE_DIR", ".")
out <- Sys.getenv(
  "TVPGVAR_GDP_LAG_COEF_OUT",
  "gdp-gpr-lag-coef-diagnostic"
)
dates <- trimws(strsplit(
  Sys.getenv(
    "TVPGVAR_GDP_LAG_COEF_DATES",
    "2003Q1,2008Q3,2020Q1,2022Q1"
  ),
  ",", fixed = TRUE
)[[1]])

dir.create(out, recursive = TRUE, showWarnings = FALSE)

find_one <- function(name, required = TRUE) {
  z <- list.files(
    src,
    pattern = paste0("^", gsub("\\.", "\\\\.", name), "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (!length(z)) {
    if (required) stop("Required file not found under ", src, ": ", name)
    return(NA_character_)
  }
  if (length(z) > 1L) {
    stop(
      "Expected exactly one ", name, " under ", src,
      "; found ", length(z), ":\n", paste(z, collapse = "\n")
    )
  }
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

qv <- function(x, prob) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, prob, names = FALSE, type = 8))
}

summarize_vec <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(data.frame(
      draws_used = 0L,
      median = NA_real_,
      low68 = NA_real_,
      high68 = NA_real_,
      low90 = NA_real_,
      high90 = NA_real_,
      positive_share = NA_real_,
      sign = NA_character_,
      credible68_excludes_zero = NA,
      credible90_excludes_zero = NA,
      stringsAsFactors = FALSE
    ))
  }

  med <- median(x)
  lo68 <- qv(x, 0.16)
  hi68 <- qv(x, 0.84)
  lo90 <- qv(x, 0.05)
  hi90 <- qv(x, 0.95)

  data.frame(
    draws_used = length(x),
    median = med,
    low68 = lo68,
    high68 = hi68,
    low90 = lo90,
    high90 = hi90,
    positive_share = mean(x > 0),
    sign = if (med > 0) "positive" else if (med < 0) "negative" else "zero",
    credible68_excludes_zero = (lo68 > 0 || hi68 < 0),
    credible90_excludes_zero = (lo90 > 0 || hi90 < 0),
    stringsAsFactors = FALSE
  )
}

# -----------------------------------------------------------------------------
# 1. Load posterior and the exact stability mask from the original run.
# -----------------------------------------------------------------------------
pred_path <- find_one("predDens_gpr_structural.rda")
irf_path <- find_one("irf_gpr_structural.rda")

load(pred_path)
if (!exists("predDens") || !exists("Data.setup")) {
  stop("predDens_gpr_structural.rda must contain predDens and Data.setup.")
}

irf_env <- new.env(parent = emptyenv())
load(irf_path, envir = irf_env)
if (!exists("stable_mask", envir = irf_env, inherits = FALSE) ||
    !exists("stability_rho", envir = irf_env, inherits = FALSE)) {
  stop("irf_gpr_structural.rda must contain stable_mask and stability_rho.")
}
stable_mask <- get("stable_mask", envir = irf_env)
stability_rho <- get("stability_rho", envir = irf_env)

countries <- as.character(Data.setup$countries)
quarters <- as.character(Data.setup$quarters)

if (!"US" %in% countries) stop("US is missing from Data.setup$countries.")

A.list <- lapply(predDens, function(z) z$ALPHA)
n_draws <- dim(A.list[[1]])[4]
if (is.null(n_draws) || n_draws < 1L) stop("No posterior draws found in ALPHA.")
if (any(vapply(A.list, function(A) dim(A)[4], integer(1)) != n_draws)) {
  stop("Countries have inconsistent posterior draw counts.")
}

# Infer domestic lag order from ALPHA tags.
all_tags <- unique(unlist(lapply(A.list, function(A) dimnames(A)[[2]])))
ylag_tags <- grep("^Ylag[0-9]+$", all_tags, value = TRUE)
if (!length(ylag_tags)) stop("No Ylag blocks found in ALPHA.")
p <- max(as.integer(sub("^Ylag", "", ylag_tags)))

irf_dates <- quarters[-seq_len(p)]
Tirf <- length(irf_dates)

if (nrow(stable_mask) != Tirf || ncol(stable_mask) != n_draws) {
  stop(
    "stable_mask dimension mismatch: got ",
    paste(dim(stable_mask), collapse = "x"),
    ", expected ", Tirf, "x", n_draws, "."
  )
}
if (!all(dim(stability_rho) == dim(stable_mask))) {
  stop("stability_rho and stable_mask dimensions differ.")
}

date_idx <- match(dates, irf_dates)
if (anyNA(date_idx)) {
  stop(
    "Requested dates not found in posterior path: ",
    paste(dates[is.na(date_idx)], collapse = ", ")
  )
}

# Optional provenance/spec check.
param_path <- find_one("gdp_lagged_gpr_test_parameters.csv", required = FALSE)
if (!is.na(param_path)) {
  prm <- readc(param_path)
  if ("gpr_lag_mode" %in% names(prm) &&
      nrow(prm) >= 1L &&
      !identical(as.character(prm$gpr_lag_mode[[1]]), "gdp_lagged_only")) {
    stop(
      "Wrong artifact: expected gpr_lag_mode=gdp_lagged_only, found ",
      prm$gpr_lag_mode[[1]]
    )
  }
  if ("p" %in% names(prm) &&
      nrow(prm) >= 1L &&
      as.integer(prm$p[[1]]) != p) {
    stop(
      "Lag-order mismatch between parameter file and ALPHA: ",
      prm$p[[1]], " vs ", p
    )
  }
}

# -----------------------------------------------------------------------------
# 2. Extract the direct lagged-GPR coefficient in each GDP equation.
# -----------------------------------------------------------------------------
raw_rows <- list()
zero_rows <- list()
us_lag2_rows <- list()
nr <- nz <- nu2 <- 0L

for (dd in seq_along(dates)) {
  d <- dates[[dd]]
  tt <- date_idx[[dd]]
  stable_draws <- which(stable_mask[tt, ])

  if (!length(stable_draws)) {
    warning("No stable posterior draws at ", d, "; skipping.")
    next
  }

  for (i in seq_along(countries)) {
    cc <- countries[[i]]
    A <- A.list[[i]]

    reg_names <- dimnames(A)[[2]]
    eq_names <- dimnames(A)[[3]]

    if (is.null(reg_names) || is.null(eq_names)) {
      stop("ALPHA lacks regressor/equation names for ", cc, ".")
    }

    y_name <- paste0(cc, "_y")
    y_eq <- match(y_name, eq_names)
    if (is.na(y_eq)) stop("GDP equation not found for ", cc, ": ", y_name)

    if (cc != "US") {
      # In gdp_lagged_only mode, non-US Wexlag1 order is:
      # foreign_y, foreign_dp, foreign_r, foreign_de, foreign_deq, global_gpr.
      lag1_idx <- which(reg_names == "Wexlag1")
      if (length(lag1_idx) != 6L) {
        stop(
          cc, ": expected 6 Wexlag1 regressors in gdp_lagged_only mode; found ",
          length(lag1_idx), "."
        )
      }
      gpr_lag1_reg <- tail(lag1_idx, 1L)
      beta1 <- as.numeric(A[tt, gpr_lag1_reg, y_eq, stable_draws])

      # Exact-zero check for the omitted contemporaneous global GPR coefficient.
      wex_idx <- which(reg_names == "Wex")
      if (length(wex_idx) != 6L) {
        stop(cc, ": expected 6 current Wex regressors; found ", length(wex_idx), ".")
      }
      current_gpr_reg <- tail(wex_idx, 1L)
      beta0 <- as.numeric(A[tt, current_gpr_reg, y_eq, stable_draws])

      for (jj in seq_along(stable_draws)) {
        dr <- stable_draws[[jj]]

        nr <- nr + 1L
        raw_rows[[nr]] <- data.frame(
          date = d,
          country = cc,
          draw = dr,
          global_rho = stability_rho[tt, dr],
          coefficient_source = "Wexlag1_global_gpr",
          beta_GPR_lag1_to_GDP = beta1[[jj]],
          stringsAsFactors = FALSE
        )

        nz <- nz + 1L
        zero_rows[[nz]] <- data.frame(
          date = d,
          country = cc,
          draw = dr,
          current_GPR_to_GDP_padded_zero = beta0[[jj]],
          stringsAsFactors = FALSE
        )
      }

    } else {
      # US endogenous ordering is expected to be:
      # US_gpr, US_y, US_dp, US_r, US_de, US_deq.
      gpr_pos <- match("US_gpr", eq_names)
      if (is.na(gpr_pos)) stop("US_gpr not found in US equation names.")

      ylag1_idx <- which(reg_names == "Ylag1")
      if (length(ylag1_idx) != length(eq_names)) {
        stop(
          "US: Ylag1 width ", length(ylag1_idx),
          " does not match endogenous dimension ", length(eq_names), "."
        )
      }
      gpr_lag1_reg <- ylag1_idx[[gpr_pos]]
      beta1 <- as.numeric(A[tt, gpr_lag1_reg, y_eq, stable_draws])

      for (jj in seq_along(stable_draws)) {
        dr <- stable_draws[[jj]]
        nr <- nr + 1L
        raw_rows[[nr]] <- data.frame(
          date = d,
          country = cc,
          draw = dr,
          global_rho = stability_rho[tt, dr],
          coefficient_source = "Ylag1_US_gpr",
          beta_GPR_lag1_to_GDP = beta1[[jj]],
          stringsAsFactors = FALSE
        )
      }

      # US alone has a direct endogenous GPR lag-2 term under p=2.
      if (p >= 2L) {
        ylag2_idx <- which(reg_names == "Ylag2")
        if (length(ylag2_idx) != length(eq_names)) {
          stop(
            "US: Ylag2 width ", length(ylag2_idx),
            " does not match endogenous dimension ", length(eq_names), "."
          )
        }
        gpr_lag2_reg <- ylag2_idx[[gpr_pos]]
        beta2 <- as.numeric(A[tt, gpr_lag2_reg, y_eq, stable_draws])

        for (jj in seq_along(stable_draws)) {
          dr <- stable_draws[[jj]]
          nu2 <- nu2 + 1L
          us_lag2_rows[[nu2]] <- data.frame(
            date = d,
            country = "US",
            draw = dr,
            global_rho = stability_rho[tt, dr],
            beta_GPR_lag2_to_GDP = beta2[[jj]],
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
}

if (!length(raw_rows)) stop("No coefficient rows were extracted.")

raw <- do.call(rbind, raw_rows)
zero_raw <- if (length(zero_rows)) do.call(rbind, zero_rows) else data.frame()
us_lag2_raw <- if (length(us_lag2_rows)) do.call(rbind, us_lag2_rows) else data.frame()

# -----------------------------------------------------------------------------
# 3. Stable-only posterior summaries.
# -----------------------------------------------------------------------------
split_key <- interaction(
  raw$date, raw$country,
  drop = TRUE, lex.order = TRUE
)

summary_rows <- lapply(split(raw, split_key), function(z) {
  s <- summarize_vec(z$beta_GPR_lag1_to_GDP)
  cbind(
    z[1, c("date", "country", "coefficient_source"), drop = FALSE],
    s
  )
})
coef_summary <- do.call(rbind, summary_rows)
rownames(coef_summary) <- NULL

# Preserve requested date order and country order.
coef_summary$date <- factor(coef_summary$date, levels = dates)
coef_summary$country <- factor(coef_summary$country, levels = countries)
coef_summary <- coef_summary[
  order(coef_summary$date, coef_summary$country),
]
coef_summary$date <- as.character(coef_summary$date)
coef_summary$country <- as.character(coef_summary$country)

# Cross-country date-level sign summary.
date_summary <- do.call(rbind, lapply(dates, function(d) {
  z <- coef_summary[coef_summary$date == d, , drop = FALSE]
  if (!nrow(z)) return(NULL)
  data.frame(
    date = d,
    countries = nrow(z),
    positive_median_countries = sum(z$median > 0, na.rm = TRUE),
    negative_median_countries = sum(z$median < 0, na.rm = TRUE),
    significant_positive_68 = sum(z$low68 > 0, na.rm = TRUE),
    significant_negative_68 = sum(z$high68 < 0, na.rm = TRUE),
    significant_positive_90 = sum(z$low90 > 0, na.rm = TRUE),
    significant_negative_90 = sum(z$high90 < 0, na.rm = TRUE),
    cross_country_median_beta = median(z$median, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

# Timing-restriction zero check for non-US GDP equations.
zero_check <- data.frame()
if (nrow(zero_raw)) {
  zk <- interaction(zero_raw$date, zero_raw$country, drop = TRUE, lex.order = TRUE)
  zero_check <- do.call(rbind, lapply(split(zero_raw, zk), function(z) {
    x <- z$current_GPR_to_GDP_padded_zero
    data.frame(
      date = z$date[[1]],
      country = z$country[[1]],
      draws_used = length(x),
      max_abs_current_GPR_to_GDP = max(abs(x)),
      all_exact_zero = all(x == 0),
      stringsAsFactors = FALSE
    )
  }))
  rownames(zero_check) <- NULL
}

# US lag-2 summary.
us_lag2_summary <- data.frame()
if (nrow(us_lag2_raw)) {
  uk <- split(us_lag2_raw, us_lag2_raw$date)
  us_lag2_summary <- do.call(rbind, lapply(uk, function(z) {
    cbind(
      data.frame(date = z$date[[1]], country = "US", stringsAsFactors = FALSE),
      summarize_vec(z$beta_GPR_lag2_to_GDP)
    )
  }))
  rownames(us_lag2_summary) <- NULL
}

# -----------------------------------------------------------------------------
# 4. Optional merge with h=1 GDP IRF summary from the same artifact.
# -----------------------------------------------------------------------------
irf_merge <- data.frame()
irf_corr <- data.frame()

gdp_diag_path <- find_one("gdp_sign_diagnostic_stable_only.csv", required = FALSE)
if (!is.na(gdp_diag_path)) {
  gdp <- readc(gdp_diag_path)

  need <- c("date", "variable", "horizon", "median", "low68", "high68", "low90", "high90")
  if (all(need %in% names(gdp))) {
    gdp$country <- sub("_y$", "", as.character(gdp$variable))
    h1 <- gdp[
      gdp$horizon == 1 &
        gdp$date %in% dates &
        gdp$country %in% countries &
        grepl("_y$", gdp$variable),
      c(
        "date", "country", "median", "low68", "high68",
        "low90", "high90", intersect(c("draws_used"), names(gdp))
      ),
      drop = FALSE
    ]

    names(h1)[names(h1) == "median"] <- "GDP_IRF_h1_median"
    names(h1)[names(h1) == "low68"] <- "GDP_IRF_h1_low68"
    names(h1)[names(h1) == "high68"] <- "GDP_IRF_h1_high68"
    names(h1)[names(h1) == "low90"] <- "GDP_IRF_h1_low90"
    names(h1)[names(h1) == "high90"] <- "GDP_IRF_h1_high90"
    if ("draws_used" %in% names(h1)) {
      names(h1)[names(h1) == "draws_used"] <- "GDP_IRF_h1_draws_used"
    }

    irf_merge <- merge(
      coef_summary,
      h1,
      by = c("date", "country"),
      all.x = TRUE,
      sort = FALSE
    )

    irf_corr <- do.call(rbind, lapply(dates, function(d) {
      z <- irf_merge[irf_merge$date == d, , drop = FALSE]
      ok <- is.finite(z$median) & is.finite(z$GDP_IRF_h1_median)
      if (sum(ok) < 3L) {
        return(data.frame(
          date = d, countries_used = sum(ok),
          pearson = NA_real_, spearman = NA_real_,
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        date = d,
        countries_used = sum(ok),
        pearson = cor(z$median[ok], z$GDP_IRF_h1_median[ok], method = "pearson"),
        spearman = cor(z$median[ok], z$GDP_IRF_h1_median[ok], method = "spearman"),
        stringsAsFactors = FALSE
      )
    }))
  } else {
    warning(
      "gdp_sign_diagnostic_stable_only.csv found, but required columns are missing. ",
      "Skipping beta-vs-h1 merge."
    )
  }
}

# -----------------------------------------------------------------------------
# 5. Write outputs.
# -----------------------------------------------------------------------------
write.csv(
  coef_summary,
  file.path(out, "01_gpr_lag1_to_gdp_summary_stable_only.csv"),
  row.names = FALSE
)
write.csv(
  raw,
  file.path(out, "02_gpr_lag1_to_gdp_draws_stable_only.csv"),
  row.names = FALSE
)
write.csv(
  date_summary,
  file.path(out, "03_gpr_lag1_to_gdp_date_sign_summary.csv"),
  row.names = FALSE
)

if (nrow(zero_check)) {
  write.csv(
    zero_check,
    file.path(out, "04_nonUS_current_GPR_to_GDP_zero_check.csv"),
    row.names = FALSE
  )
}

if (nrow(us_lag2_summary)) {
  write.csv(
    us_lag2_summary,
    file.path(out, "05_US_gpr_lag2_to_gdp_summary_stable_only.csv"),
    row.names = FALSE
  )
}

if (nrow(irf_merge)) {
  write.csv(
    irf_merge,
    file.path(out, "06_gpr_lag1_beta_vs_GDP_h1_IRF.csv"),
    row.names = FALSE
  )
}
if (nrow(irf_corr)) {
  write.csv(
    irf_corr,
    file.path(out, "07_gpr_lag1_beta_vs_GDP_h1_IRF_correlation_by_date.csv"),
    row.names = FALSE
  )
}

# Compact human-readable report.
report <- c(
  "Direct GPR lag-1 -> GDP coefficient diagnostic",
  paste0("Source directory: ", normalizePath(src, mustWork = FALSE)),
  paste0("Domestic lag order p: ", p),
  paste0("Posterior draws: ", n_draws),
  paste0("Dates: ", paste(dates, collapse = ", ")),
  "",
  "Date-level sign summary:",
  paste(capture.output(print(date_summary, row.names = FALSE)), collapse = "\n"),
  "",
  "Country/date coefficient summaries (sorted by date, descending median):"
)

for (d in dates) {
  z <- coef_summary[coef_summary$date == d, , drop = FALSE]
  z <- z[order(-z$median), ]
  report <- c(
    report,
    "",
    paste0("===== ", d, " ====="),
    paste(capture.output(print(
      z[, c(
        "country", "draws_used", "median", "low68", "high68",
        "low90", "high90", "positive_share",
        "credible68_excludes_zero", "credible90_excludes_zero"
      )],
      row.names = FALSE
    )), collapse = "\n")
  )
}

if (nrow(irf_corr)) {
  report <- c(
    report,
    "",
    "Cross-country association between lag-1 beta and h=1 GDP IRF:",
    paste(capture.output(print(irf_corr, row.names = FALSE)), collapse = "\n")
  )
}

writeLines(report, file.path(out, "README_diagnostic_summary.txt"))

cat("\nDiagnostic complete.\n")
cat("Output directory:", out, "\n\n")
print(date_summary, row.names = FALSE)

cat("\nKey interpretation:\n")
cat("* Mostly positive/significant beta_GPR_lag1_to_GDP -> direct lagged-GPR channel is a main source of positive h=1 GDP IRFs.\n")
cat("* Beta near zero/negative while h=1 IRF remains positive -> contemporaneous non-GDP propagation / structural shock identification is more important.\n")
cat("* Non-US zero-check must report all_exact_zero=TRUE under gdp_lagged_only.\n")
