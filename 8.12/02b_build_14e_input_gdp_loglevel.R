#!/usr/bin/env Rscript

# =============================================================================
# GDP LOG-LEVEL A/B EXPERIMENT
#
# Purpose
# -------
# Keep the existing validated macro/GPR/trade pipeline unchanged, but replace
# each country's GDP growth variable (GDP_DLOG -> *_y) with an indexed log-level
# series reconstructed from the SAME validated GDP_DLOG observations.
#
# This is deliberately an A/B transformation experiment:
#   baseline: y_t = Delta log(real GDP_t)
#   test:     y_t = log(indexed real GDP_t), index at first sample = 100
#
# If GDP_DLOG = log(Y_t)-log(Y_{t-1}), then for t >= 2:
#   log_index_t - log_index_{t-1} = GDP_DLOG_t
# up to the user-supplied scale factor.
#
# No GDP observations are invented. The first in-sample level is simply
# normalized to 100; this adds a country-specific constant that is absorbed by
# the intercept. For a formal paper version, raw real-GDP levels can later be
# used directly if desired.
# =============================================================================

options(stringsAsFactors = FALSE)

baseline_builder <- Sys.getenv(
  "TVPGVAR_BASELINE_BUILDER",
  "8.12/02_build_14e_input_structural.R"
)
if (!file.exists(baseline_builder)) stop("Missing baseline builder: ", baseline_builder)

# Build the exact current baseline first (macro + GPR + trade weights).
source(baseline_builder, local = FALSE)

input_file <- Sys.getenv("TVPGVAR_BASE_INPUT", "data/model_input.csv")
if (!file.exists(input_file)) stop("Baseline builder did not create: ", input_file)

countries <- c("AU","BR","CA","CH","CN","EA","UK","JP","KR","NO","SG","TR","US","ZA")
scale_factor <- as.numeric(Sys.getenv("TVPGVAR_GDP_DLOG_SCALE", "1"))
base_index <- as.numeric(Sys.getenv("TVPGVAR_GDP_BASE_INDEX", "100"))

if (!is.finite(scale_factor) || scale_factor <= 0) stop("TVPGVAR_GDP_DLOG_SCALE must be > 0.")
if (!is.finite(base_index) || base_index <= 0) stop("TVPGVAR_GDP_BASE_INDEX must be > 0.")

model <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!"Quarter" %in% names(model)) stop("model_input.csv is missing Quarter.")

# Save an untouched copy for exact A/B auditing.
dir.create("data", showWarnings = FALSE, recursive = TRUE)
write.csv(model, "data/model_input_before_gdp_loglevel.csv", row.names = FALSE)

# Detect an obvious scale mismatch before transforming.
all_dlog <- unlist(lapply(countries, function(cc) {
  nm <- paste0(cc, "_y")
  if (!nm %in% names(model)) stop("Missing GDP variable: ", nm)
  as.numeric(model[[nm]])
}), use.names = FALSE)
all_dlog <- all_dlog[is.finite(all_dlog)]

if (!length(all_dlog)) stop("No finite GDP_DLOG values found.")
if (scale_factor == 1 && quantile(abs(all_dlog), 0.99, na.rm = TRUE) > 0.75) {
  stop(
    "GDP_DLOG appears too large for raw log differences. ",
    "If the source is 100*Delta log(GDP), rerun with TVPGVAR_GDP_DLOG_SCALE=100."
  )
}

# Reconstruction audit.
audit <- list()
log_base <- log(base_index)

for (cc in countries) {
  nm <- paste0(cc, "_y")
  dlog <- as.numeric(model[[nm]]) / scale_factor
  if (any(!is.finite(dlog))) stop("Non-finite GDP_DLOG for ", cc)

  # First in-sample GDP level is normalized to index=100.
  # dlog[1] is the pre-sample -> first-sample change and cannot identify the
  # absolute first level; it is intentionally not used to move the anchor.
  level_log <- numeric(length(dlog))
  level_log[1] <- log_base
  if (length(dlog) >= 2L) {
    level_log[2:length(dlog)] <- log_base + cumsum(dlog[2:length(dlog)])
  }

  # Exact reconstruction check for all within-sample changes.
  err <- if (length(dlog) >= 2L) max(abs(diff(level_log) - dlog[-1])) else 0
  if (!is.finite(err) || err > 1e-12) {
    stop("GDP log-level reconstruction failed for ", cc, "; max error=", err)
  }

  model[[nm]] <- level_log
  audit[[cc]] <- data.frame(
    country = cc,
    first_quarter = model$Quarter[1],
    last_quarter = tail(model$Quarter, 1),
    base_index = base_index,
    dlog_scale_factor = scale_factor,
    first_log_level = level_log[1],
    last_log_level = tail(level_log, 1),
    implied_last_index = exp(tail(level_log, 1)),
    max_diff_reconstruction_error = err,
    stringsAsFactors = FALSE
  )
}

audit <- do.call(rbind, audit)
write.csv(audit, "results/gdp_loglevel_reconstruction_audit.csv", row.names = FALSE)

# Overwrite ONLY the model inputs consumed by downstream VIX/dominant scripts.
write.csv(model, "data/model_input.csv", row.names = FALSE)
write.csv(model, "data/model_input_structural.csv", row.names = FALSE)

# Append explicit provenance/method notes.
summary_file <- "results/build_structural_input_summary.txt"
extra <- c(
  "",
  "GDP LOG-LEVEL EXPERIMENT",
  "GDP transformation: indexed log level reconstructed from the validated GDP_DLOG series",
  paste0("GDP first-sample index normalization: ", base_index),
  paste0("GDP_DLOG scale divisor: ", scale_factor),
  "For t>=2, diff(log-level index) exactly equals GDP_DLOG / scale divisor.",
  "The first in-sample log level is an arbitrary normalization; the intercept absorbs the additive constant.",
  "This is an A/B specification test; all non-GDP macro series, GPR, VIX, weights, priors and identification remain unchanged."
)
cat(extra, sep = "\n", file = summary_file, append = TRUE)
cat(paste(extra, collapse = "\n"), "\n")

cat("GDP log-level input created successfully: ", input_file, "\n", sep = "")
