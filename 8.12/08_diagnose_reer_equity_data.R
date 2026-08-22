#!/usr/bin/env Rscript

# =============================================================================
# TVP-GVAR DATA SANITY DIAGNOSTIC: REER + EQUITY + VIX
#
# Purpose
# -------
# Diagnose the two remaining suspicious patterns WITHOUT re-running MCMC:
#   1) ZA exchange-rate (REER) response is much larger than other economies.
#   2) Equity IRF medians are mostly positive after a GPR shock.
#
# This script checks the ACTUAL model input used by the successful GDP-loglevel
# dominant-unit run. In the current builder:
#   *_de  <- REER_DLOG
#   *_deq <- EQ_RETURN
#   GL_vix = log(quarterly mean VIX)
#
# It does NOT change data and does NOT estimate the TVP-GVAR.
# Base R only: no external R packages are required.
# =============================================================================

options(stringsAsFactors = FALSE)

input_file <- Sys.getenv(
  "TVPGVAR_SANITY_INPUT",
  "prior_artifact/data/model_input_dominant.csv"
)
out_dir <- Sys.getenv(
  "TVPGVAR_SANITY_OUT",
  "results/data_sanity_diagnostic"
)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) stop("Missing model input: ", input_file)

dat <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
required <- c(
  "Quarter", "GL_gpr", "GL_vix",
  unlist(lapply(countries, function(cc) c(paste0(cc, "_de"), paste0(cc, "_deq"))), use.names = FALSE)
)
missing_cols <- setdiff(required, names(dat))
if (length(missing_cols)) stop("Missing required columns: ", paste(missing_cols, collapse = ", "))

if (anyDuplicated(dat$Quarter)) stop("Duplicate quarters detected.")

quarter_index <- function(q) {
  q <- as.character(q)
  if (any(!grepl("^[12][0-9]{3}Q[1-4]$", q))) stop("Invalid quarter label detected.")
  as.integer(substr(q, 1, 4)) * 4L + as.integer(substr(q, 6, 6))
}
qidx <- quarter_index(dat$Quarter)
if (is.unsorted(qidx)) stop("Quarters are not sorted chronologically.")
if (any(diff(qidx) != 1L)) warning("Sample contains quarter gaps.")

num_required <- setdiff(required, "Quarter")
for (nm in num_required) dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
if (any(!is.finite(as.matrix(dat[, num_required, drop = FALSE])))) {
  stop("NA/NaN/Inf found in required model-input columns.")
}

# VIX is log(level), so the change relevant for an equity-return sign sanity
# check is Delta log(VIX), not the VIX level itself.
dat$d_GL_vix <- c(NA_real_, diff(dat$GL_vix))
dat$d_GL_gpr <- c(NA_real_, diff(dat$GL_gpr))

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 8L || sd(x[ok]) == 0 || sd(y[ok]) == 0) return(NA_real_)
  cor(x[ok], y[ok])
}

# -----------------------------------------------------------------------------
# 1) Cross-country scale diagnostic for REER changes and equity returns
# -----------------------------------------------------------------------------
stat_one <- function(x) {
  x <- x[is.finite(x)]
  c(
    n = length(x),
    mean = mean(x),
    sd = sd(x),
    median = median(x),
    min = min(x),
    q01 = unname(quantile(x, .01)),
    q05 = unname(quantile(x, .05)),
    q95 = unname(quantile(x, .95)),
    q99 = unname(quantile(x, .99)),
    max = max(x),
    max_abs = max(abs(x))
  )
}

scale_rows <- list()
rr <- 0L
for (concept in c("de", "deq")) {
  for (cc in countries) {
    nm <- paste0(cc, "_", concept)
    ss <- stat_one(dat[[nm]])
    rr <- rr + 1L
    scale_rows[[rr]] <- data.frame(
      country = cc,
      concept = concept,
      variable = nm,
      n = as.integer(ss[["n"]]),
      mean = ss[["mean"]],
      sd = ss[["sd"]],
      median = ss[["median"]],
      min = ss[["min"]],
      q01 = ss[["q01"]],
      q05 = ss[["q05"]],
      q95 = ss[["q95"]],
      q99 = ss[["q99"]],
      max = ss[["max"]],
      max_abs = ss[["max_abs"]],
      stringsAsFactors = FALSE
    )
  }
}
scale_tab <- do.call(rbind, scale_rows)

for (concept in c("de", "deq")) {
  ix <- scale_tab$concept == concept
  med_sd <- median(scale_tab$sd[ix], na.rm = TRUE)
  med_ma <- median(scale_tab$max_abs[ix], na.rm = TRUE)
  scale_tab$sd_ratio_to_cross_country_median[ix] <- scale_tab$sd[ix] / med_sd
  scale_tab$maxabs_ratio_to_cross_country_median[ix] <- scale_tab$max_abs[ix] / med_ma
}

scale_tab$scale_flag <- "OK"
scale_tab$scale_flag[
  scale_tab$sd_ratio_to_cross_country_median >= 5 |
  scale_tab$maxabs_ratio_to_cross_country_median >= 5
] <- "VERY_LARGE"
scale_tab$scale_flag[
  scale_tab$scale_flag == "OK" & (
    scale_tab$sd_ratio_to_cross_country_median >= 3 |
    scale_tab$maxabs_ratio_to_cross_country_median >= 3
  )
] <- "LARGE"
scale_tab$scale_flag[
  scale_tab$sd_ratio_to_cross_country_median <= .20 &
  scale_tab$maxabs_ratio_to_cross_country_median <= .20
] <- "VERY_SMALL"

write.csv(
  scale_tab,
  file.path(out_dir, "scale_summary_reer_equity.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 2) Largest absolute observations by country and variable
# -----------------------------------------------------------------------------
extreme_rows <- list()
rr <- 0L
for (concept in c("de", "deq")) {
  for (cc in countries) {
    nm <- paste0(cc, "_", concept)
    ord <- order(abs(dat[[nm]]), decreasing = TRUE)
    ord <- head(ord, 10L)
    for (rank in seq_along(ord)) {
      ii <- ord[[rank]]
      rr <- rr + 1L
      extreme_rows[[rr]] <- data.frame(
        country = cc,
        concept = concept,
        variable = nm,
        rank_abs = rank,
        Quarter = dat$Quarter[ii],
        value = dat[[nm]][ii],
        abs_value = abs(dat[[nm]][ii]),
        stringsAsFactors = FALSE
      )
    }
  }
}
extreme_tab <- do.call(rbind, extreme_rows)
write.csv(
  extreme_tab,
  file.path(out_dir, "top10_extreme_observations_reer_equity.csv"),
  row.names = FALSE
)

# Dedicated ZA REER audit.
za_de <- subset(extreme_tab, country == "ZA" & concept == "de")
write.csv(
  za_de,
  file.path(out_dir, "ZA_REER_top10_extremes.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 3) Equity sign sanity checks in well-known stress quarters
# -----------------------------------------------------------------------------
stress_quarters <- c(
  "2008Q3", "2008Q4",
  "2020Q1", "2020Q2",
  "2022Q1", "2022Q2"
)
stress_quarters <- stress_quarters[stress_quarters %in% dat$Quarter]

crash_rows <- list()
rr <- 0L
for (qq in stress_quarters) {
  ii <- match(qq, dat$Quarter)
  for (cc in countries) {
    nm <- paste0(cc, "_deq")
    rr <- rr + 1L
    crash_rows[[rr]] <- data.frame(
      Quarter = qq,
      country = cc,
      equity_return = dat[[nm]][ii],
      sign = ifelse(dat[[nm]][ii] > 0, "positive",
                    ifelse(dat[[nm]][ii] < 0, "negative", "zero")),
      stringsAsFactors = FALSE
    )
  }
}
crash_long <- do.call(rbind, crash_rows)
write.csv(
  crash_long,
  file.path(out_dir, "equity_stress_quarter_returns.csv"),
  row.names = FALSE
)

crash_summary <- do.call(rbind, lapply(split(crash_long, crash_long$Quarter), function(d) {
  data.frame(
    Quarter = d$Quarter[1],
    negative = sum(d$equity_return < 0),
    positive = sum(d$equity_return > 0),
    zero = sum(d$equity_return == 0),
    cross_country_median = median(d$equity_return),
    cross_country_mean = mean(d$equity_return),
    economies = nrow(d),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  crash_summary,
  file.path(out_dir, "equity_stress_quarter_sign_summary.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 4) Equity/VIX correlation diagnostic
# -----------------------------------------------------------------------------
cor_rows <- lapply(countries, function(cc) {
  x <- dat[[paste0(cc, "_deq")]]
  data.frame(
    country = cc,
    cor_equity_with_logVIX_level = safe_cor(x, dat$GL_vix),
    cor_equity_with_dlogVIX = safe_cor(x, dat$d_GL_vix),
    cor_equity_with_logGPR_level = safe_cor(x, dat$GL_gpr),
    cor_equity_with_dlogGPR = safe_cor(x, dat$d_GL_gpr),
    stringsAsFactors = FALSE
  )
})
cor_tab <- do.call(rbind, cor_rows)
write.csv(
  cor_tab,
  file.path(out_dir, "equity_vix_gpr_correlations.csv"),
  row.names = FALSE
)

# -----------------------------------------------------------------------------
# 5) Automatic diagnostic flags
# -----------------------------------------------------------------------------
za_scale <- subset(scale_tab, country == "ZA" & concept == "de")
equity_scale_bad <- subset(scale_tab, concept == "deq" & scale_flag != "OK")
reer_scale_bad <- subset(scale_tab, concept == "de" & scale_flag != "OK")

# Strong sign-convention evidence: the classic crash quarters should usually
# have broad negative equity returns, and equity returns should generally be
# negatively correlated with Delta log(VIX).
classic <- crash_summary[crash_summary$Quarter %in% c("2008Q4", "2020Q1"), , drop = FALSE]
classic_negative_ok <- nrow(classic) > 0 && all(classic$negative >= 10L)
median_cor_dvix <- median(cor_tab$cor_equity_with_dlogVIX, na.rm = TRUE)
cor_negative_ok <- is.finite(median_cor_dvix) && median_cor_dvix < 0

if (classic_negative_ok && cor_negative_ok) {
  equity_sign_conclusion <- paste0(
    "LIKELY_CORRECT: equity input sign convention looks economically normal. ",
    "Classic crash quarters are broadly negative and median cor(EQ_RETURN, Delta log VIX) = ",
    sprintf("%.3f", median_cor_dvix),
    ". If structural GPR equity IRFs remain mostly positive, focus next on identification/mapping rather than flipping the raw equity sign."
  )
} else if (cor_negative_ok) {
  equity_sign_conclusion <- paste0(
    "MIXED_BUT_NOT_REVERSED: median cor(EQ_RETURN, Delta log VIX) is negative (",
    sprintf("%.3f", median_cor_dvix),
    "), so a wholesale equity-sign reversal is not supported. Inspect stress-quarter values and source construction before changing transformations."
  )
} else {
  equity_sign_conclusion <- paste0(
    "REQUIRES_REVIEW: equity returns are not showing the expected broad negative relation with Delta log VIX. Median correlation = ",
    sprintf("%.3f", median_cor_dvix),
    ". Check the EQ_RETURN construction/sign before touching model identification."
  )
}

za_conclusion <- if (nrow(za_scale) == 1L && za_scale$scale_flag != "OK") {
  paste0(
    "FLAGGED: ZA_REER_DLOG is a cross-country scale/outlier concern. sd ratio = ",
    sprintf("%.2f", za_scale$sd_ratio_to_cross_country_median),
    "x; max-abs ratio = ",
    sprintf("%.2f", za_scale$maxabs_ratio_to_cross_country_median),
    "x. Inspect ZA_REER_top10_extremes.csv and the original REER source/transform."
  )
} else {
  paste0(
    "NOT_SCALE_FLAGGED: ZA_REER_DLOG is not >=3x the cross-country median by SD or max absolute observation. ",
    "Its unusual IRF may therefore come from model coefficients/dynamics rather than a simple input scale mismatch."
  )
}

# Compact machine-readable verdict table.
verdict <- data.frame(
  diagnostic = c(
    "equity_input_sign",
    "ZA_REER_scale",
    "equity_cross_country_scale_flags",
    "REER_cross_country_scale_flags"
  ),
  result = c(
    equity_sign_conclusion,
    za_conclusion,
    paste(nrow(equity_scale_bad), "equity series flagged LARGE/VERY_LARGE/VERY_SMALL"),
    paste(nrow(reer_scale_bad), "REER series flagged LARGE/VERY_LARGE/VERY_SMALL")
  ),
  stringsAsFactors = FALSE
)
write.csv(verdict, file.path(out_dir, "automatic_verdicts.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# 6) Human-readable report
# -----------------------------------------------------------------------------
fmt_quarters <- if (nrow(crash_summary)) {
  paste(apply(crash_summary, 1, function(z) {
    paste0(z[["Quarter"]], ": negative=", z[["negative"]],
           ", positive=", z[["positive"]],
           ", median=", sprintf("%.5f", as.numeric(z[["cross_country_median"]])))
  }), collapse = "\n")
} else "No requested stress quarters found."

report <- c(
  "TVP-GVAR REER / EQUITY SANITY DIAGNOSTIC",
  "========================================",
  paste0("Input: ", input_file),
  paste0("Sample: ", dat$Quarter[1], " - ", tail(dat$Quarter, 1)),
  paste0("Observations: ", nrow(dat)),
  "",
  "IMPORTANT",
  "- This is a data/post-processing diagnostic only. NO MCMC was run.",
  "- *_de is interpreted as REER_DLOG from the current project builder.",
  "- *_deq is interpreted as EQ_RETURN from the current project builder.",
  "- GL_vix is log quarterly-mean VIX; correlations also use Delta log(VIX).",
  "",
  "EQUITY SIGN VERDICT",
  equity_sign_conclusion,
  "",
  "ZA REER VERDICT",
  za_conclusion,
  "",
  "STRESS-QUARTER EQUITY SIGNS",
  fmt_quarters,
  "",
  paste0("Median country cor(EQ_RETURN, Delta log VIX): ", sprintf("%.4f", median_cor_dvix)),
  paste0("Equity cross-country scale flags: ", nrow(equity_scale_bad)),
  paste0("REER cross-country scale flags: ", nrow(reer_scale_bad)),
  "",
  "NEXT INTERPRETATION RULE",
  "If equity input signs look normal here but structural equity IRFs stay positive, do NOT flip EQ_RETURN just to obtain negative IRFs. The next diagnostic should inspect structural shock propagation / contemporaneous and lagged mapping.",
  "If ZA_REER_DLOG is strongly scale-flagged, inspect the original ZA REER series and transformation before re-estimating the model."
)
writeLines(report, file.path(out_dir, "README_diagnostic.txt"))
cat(paste(report, collapse = "\n"), "\n")

cat("\nFiles written to: ", out_dir, "\n", sep = "")
